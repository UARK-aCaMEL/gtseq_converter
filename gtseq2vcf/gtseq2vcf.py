#!/usr/bin/env python3
import argparse
import logging
import subprocess
from pathlib import Path
import pandas as pd
import pysam


def setup_logger(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    fh = logging.FileHandler(log_file)
    fh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s"))
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s"))
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


class GTseqToVCF:
    def __init__(self, csv_file, prefix, drop_patterns=None, radseq_file=None):
        self.csv_file = Path(csv_file)
        self.prefix = prefix
        self.drop_patterns = drop_patterns or []
        self.radseq_file = Path(radseq_file) if radseq_file else None
        self.logger = setup_logger(Path(f"{self.prefix}.log"))

    def load_and_parse(self):
        self.logger.info(f"Loading GTseq CSV: {self.csv_file}")
        df = pd.read_csv(self.csv_file)

        # Drop known metadata columns if present
        meta_cols = [
            "Total Genotypes",
            "Total Missing",
            "%missing",
            "Raw Reads",
            "On-Target Reads",
            "%On-Target",
            "%GT",
            "IFI",
        ]
        df.drop(columns=[c for c in meta_cols if c in df.columns],
                inplace=True, errors='ignore')

        # Drop unwanted patterns
        for p in self.drop_patterns:
            df = df.loc[:, ~df.columns.str.contains(p, case=False)]
        df.columns = df.columns.str.strip()

        # parse loci into DataFrame
        sample_col = df.columns[0]
        loci = df.columns[1:]
        chrom, pos = zip(*(c.split("_", 1) for c in loci))
        df2 = df.set_index(sample_col).T
        df2.insert(0, "CHROM", chrom)
        df2.insert(1, "POS", pos)
        df2.index = [f"{c}_{p}" for c, p in zip(chrom, pos)]

        # resolve duplicate sample columns by most non-missing genotypes
        cols = list(df2.columns)
        col_positions = {}
        for idx, name in enumerate(cols):
            col_positions.setdefault(name, []).append(idx)
        keep_idxs = []
        for name, idxs in col_positions.items():
            if len(idxs) == 1:
                keep_idxs.append(idxs[0])
            else:
                self.logger.info(f"Found {len(idxs)} duplicate columns for sample '{name}': {idxs}")
                best_idx, best_count = idxs[0], -1
                for i in idxs:
                    count = df2.iloc[:, i].apply(
                        lambda g: any(ch not in {'-', 'N', '0', '00'} for ch in str(g))
                    ).sum()
                    self.logger.info(f"  Column {i} has {count} non-missing genotypes")
                    if count > best_count:
                        best_count, best_idx = count, i
                self.logger.info(f"Keeping column {best_idx} for '{name}' ({best_count} non-missing)")
                keep_idxs.append(best_idx)
        df2 = df2.iloc[:, keep_idxs]

        self.df2 = df2
        self.snps = df2.index.tolist()
        self.samples = [c for c in df2.columns if c not in ('CHROM', 'POS')]
        self.logger.info(f"Parsed {len(self.snps)} loci for {len(self.samples)} samples.")

        # determine REF/ALT from GTseq
        from collections import Counter
        refs, alts = {}, {}
        for snpid in self.snps:
            bases = [b for gt in df2.loc[snpid, self.samples] for b in gt if b not in {'-', 'N', '0', '00'}]
            cnt = Counter(bases)
            refs[snpid] = cnt.most_common(1)[0][0] if cnt else '.'
            alts[snpid] = cnt.most_common(2)[1][0] if len(cnt) > 1 else '.'

        # override REF/ALT with radseq if provided
        if self.radseq_file:
            self.logger.info(f"Parsing RADseq VCF for REF/ALT: {self.radseq_file}")
            vcf = pysam.VariantFile(str(self.radseq_file))
            rad_samples = len(vcf.header.samples)
            rad_loci = sum(1 for _ in vcf.fetch())
            self.logger.info(f"RADseq input: {rad_samples} samples, {rad_loci} loci")
            for rec in vcf.fetch():
                key = f"{rec.chrom}_{rec.pos}"
                if key in refs:
                    refs[key], alts[key] = rec.ref, rec.alts[0] if rec.alts else '.'
                    self.logger.info(f"Using RADseq REF/ALT at {key}: {refs[key]}/{alts[key]}")

        self.refs, self.alts = refs, alts

    def write_plain_vcf(self):
        chroms = sorted(set(self.df2['CHROM']))
        contig_lines = [f"##contig=<ID={c}>" for c in chroms]
        header = [
            "##fileformat=VCFv4.2",
            "##source=GTseqToVCF",
            "##reference=GenomeRef",
        ] + contig_lines + [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(self.samples)
        ]
        tmp = Path(f"{self.prefix}.tmp.vcf")
        self.logger.info(f"Writing temporary VCF: {tmp}")
        with open(tmp, 'w') as f:
            f.write("\n".join(header) + "\n")
            for snpid in self.snps:
                chrom, pos = snpid.split("_", 1)
                r, a = self.refs[snpid], self.alts[snpid]
                gt = ["/".join('0' if b==r else '1' if b==a else '.' for b in self.df2.loc[snpid, self.samples][s]) for s in self.samples]
                line = [chrom, pos, snpid, r, a, '.', 'PASS', '.', 'GT'] + gt
                f.write("\t".join(map(str, line)) + "\n")
        return tmp

    def sort_compress_index(self, vcf_in, vcf_out):
        self.logger.info(f"Sorting and compressing {vcf_in} -> {vcf_out}")
        subprocess.run(['bcftools', 'sort', '-Oz', '-o', str(vcf_out), str(vcf_in)], check=True)
        subprocess.run(['tabix', '-p', 'vcf', str(vcf_out)], check=True)
        vcf_in.unlink()

    def merge_with_radseq(self, tmp_vcf):
        final = Path(f"{self.prefix}.vcf.gz")
        self.logger.info(f"Merging with RADseq {self.radseq_file} -> {final}")

        # --- NEW: prune RADseq samples that duplicate GTseq samples (keep GTseq) ---
        # If there is any overlap in sample names, drop those from RADseq before merging.
        pruned_rad = Path(f"{self.prefix}.radseq_pruned.vcf.gz")
        if self.samples:
            exclude_list = ",".join(self.samples)  # names to exclude from RADseq
            self.logger.info(f"Pruning {self.radseq_file} to drop {len(self.samples)} GTseq samples if present.")
            # bcftools view -s ^S1,S2,... keeps all samples except listed
            subprocess.run([
                'bcftools', 'view', '-s', f'^{exclude_list}', '-Oz',
                '-o', str(pruned_rad), str(self.radseq_file)
            ], check=True)
            subprocess.run(['tabix', '-p', 'vcf', str(pruned_rad)], check=True)
            rad_for_merge = pruned_rad
        else:
            rad_for_merge = self.radseq_file

        subprocess.run([
            'bcftools', 'merge', '-m', 'all', '-O', 'z',
            str(rad_for_merge), str(tmp_vcf), '-o', str(final)
        ], check=True)
        subprocess.run(['tabix', '-p', 'vcf', str(final)], check=True)

        # Clean up pruned RAD if created
        if pruned_rad.exists():
            pruned_rad.unlink()
            pruned_tbi = pruned_rad.with_suffix(pruned_rad.suffix + '.tbi')
            if pruned_tbi.exists():
                pruned_tbi.unlink()

        tmp_vcf.unlink()
        return final


def get_args():
    p = argparse.ArgumentParser(description="GTseq→VCF (+ optional RADseq merge)")
    p.add_argument('--gtseq', required=True)
    p.add_argument('--prefix', required=True)
    p.add_argument('--drop', nargs='*', default=[] )
    p.add_argument('--radseq', help='BGZIPped+indexed VCF for RADseq')
    return p.parse_args()


def main():
    args = get_args()
    conv = GTseqToVCF(args.gtseq, args.prefix, args.drop, args.radseq)
    conv.load_and_parse()
    tmp = conv.write_plain_vcf()
    sorted_tmp = Path(f"{args.prefix}.tmp.vcf.gz")
    conv.sort_compress_index(tmp, sorted_tmp)
    # report GTseq-only file stats
    v_gt = pysam.VariantFile(str(sorted_tmp))
    gt_samples = len(v_gt.header.samples)
    gt_loci = sum(1 for _ in v_gt.fetch())
    conv.logger.info(f"GTseq-only VCF {sorted_tmp}: {gt_samples} samples, {gt_loci} loci")
    if args.radseq:
        final = conv.merge_with_radseq(sorted_tmp)
        v_fin = pysam.VariantFile(str(final))
        fin_samples = len(v_fin.header.samples)
        fin_loci = sum(1 for _ in v_fin.fetch())
        conv.logger.info(f"Final merged VCF {final}: {fin_samples} samples, {fin_loci} loci")
    else:
        conv.logger.info(f"Final GTseq VCF remains {sorted_tmp}")

if __name__ == '__main__':
    main()
