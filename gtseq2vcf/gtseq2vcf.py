#!/usr/bin/env python3
import argparse
import gzip
import logging
import shutil
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns
from scipy.stats import pearsonr

fontsize = 18
dpi = 300

# Adjust matplotlib settings globally.
sizes = {
    "axes.labelsize": fontsize,
    "axes.titlesize": fontsize,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "xtick.top": False,
    "ytick.right": False,
    "figure.titlesize": fontsize,
    "figure.labelsize": fontsize,
    "xtick.labelsize": fontsize,
    "ytick.labelsize": fontsize,
    "font.size": fontsize,
    "legend.fontsize": fontsize,
    "legend.title_fontsize": fontsize,
    "legend.frameon": False,
    "legend.markerscale": 2.0,
    "figure.dpi": dpi,
    "savefig.dpi": dpi,
}

sns.set_context("paper", rc=sizes)
plt.rcParams.update(sizes)
mpl.rcParams.update(sizes)


def setup_logger(log_file):
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.handlers = []
    Path(log_file).parents[0].mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(log_file)
    sh = logging.StreamHandler()
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(name)s - %(message)s")
    fh.setFormatter(fmt)
    sh.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(sh)


class GTseqToVCF:
    def __init__(self, filepath, output_dir, prefix, plot_filetype, str2drop=None, debug=False):
        self.filepath = filepath
        self.prefix = str(prefix)
        self.plt_ft = plot_filetype.lstrip(".").lower()
        self.debug = debug
        self.data = None
        self.header = []
        self.vcf_data = []
        self.snp_ids = []
        self.sample_ids = []
        self.header_key = "Sample"

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        logpath = self.output_dir / f"{self.prefix}_logfile.txt"
        setup_logger(logpath)
        self.logger = logging.getLogger(__name__)

        self.metadata_cols = ["Raw Reads", "On-Target Reads", "%On-Target", "%GT", "IFI"]
        self.str2drop = str2drop

    def load_data(self):
        self.logger.info("Loading GT-seq CSV file...")
        df = pd.read_csv(self.filepath, header=0)
        self.metadata = pd.concat([df.pop(c) for c in self.metadata_cols if c in df.columns], axis=1)
        df.columns = df.columns.str.strip()
        if self.str2drop:
            patterns = self.str2drop if isinstance(self.str2drop, list) else [self.str2drop]
            for p in patterns:
                df = df.loc[:, ~df.columns.str.lower().str.contains(p.lower())]
        self.data = df
        self.logger.info("Successfully loaded GT-seq CSV file!")

    def parse_sample_column(self):
        self.logger.info("Parsing loci and sample data from GT-seq file...")
        sample_col = self.data.columns[0]
        locus_cols = [c for c in self.data.columns if "_" in c and c != sample_col]
        if not locus_cols:
            raise ValueError("No locus columns found (expecting 'CHR_POS').")
        chrom, pos = zip(*(c.split("_", 1) for c in locus_cols))
        df = self.data[[sample_col] + locus_cols].set_index(sample_col).T
        df.insert(0, "CHROM", chrom)
        df.insert(1, "POS", pos)
        df.index = pd.MultiIndex.from_arrays([df["CHROM"], df["POS"]], names=("CHROM", "POS"))
        df = df.drop(columns=["CHROM", "POS"])
        if not df.index.get_level_values("POS").str.isnumeric().all():
            msg = "Non-numeric values in POS level"
            self.logger.error(msg)
            raise ValueError(msg)
        self.parsed = df.copy()
        self.logger.info("Successfully parsed loci and sample data!")

    def calculate_ref_alt_alleles(self):
        self.logger.info("Calculating reference and alternate alleles...")
        df = self.parsed
        def get_ref_alt(row):
            counts = Counter(base for genotype in row for base in genotype if base not in {"-", "N", "0"})
            if counts:
                common = counts.most_common()
                return common[0][0], (common[1][0] if len(common) > 1 else ".")
            else:
                return None, "."
        refs_alts = df.apply(get_ref_alt, axis=1)
        self.ref_dict = {idx: ra[0] for idx, ra in refs_alts.items()}
        self.alt_dict = {idx: ra[1] for idx, ra in refs_alts.items()}
        self.logger.info("REF/ALT dictionaries built.")

    def transpose_data(self):
        df = self.parsed.reset_index()
        df = df.rename(columns={"level_0": "CHROM", "level_1": "POS"})
        df["ID"] = df["CHROM"].astype(str) + "_" + df["POS"].astype(str)
        df = df.set_index("ID")
        df = df.loc[:, ~df.columns.duplicated()]
        self.data = df
        self.snp_ids = list(df.index)
        seen = set()
        unique = []
        for c in df.columns:
            if c not in ["CHROM", "POS"] and c not in seen:
                unique.append(c)
                seen.add(c)
        self.sample_ids = unique

    def create_vcf_header(self):
        self.logger.info("Creating VCF header...")
        samples = "\t".join(self.sample_ids)
        if not samples:
            raise ValueError("No samples found for VCF header.")
        self.header = [
            "##fileformat=VCFv4.2",
            "##source=GTseqToVCFConverter",
            "##reference=GenomeRef",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samples,
        ]

    def format_genotype_data(self):
        self.logger.info("Formatting genotype data...")
        geno_matrix = self.data[self.sample_ids].values
        refs = [self.ref_dict.get(idx, ".") for idx in self.snp_ids]
        alts = [self.alt_dict.get(idx, ".") for idx in self.snp_ids]
        def encode_row(row, ref, alt):
            return ["/".join("0" if b == ref else "1" if b == alt else "." for b in gt) for gt in row]
        formatted = [encode_row(row, r, a) for row, r, a in zip(geno_matrix, refs, alts)]
        for idx, genos in zip(self.snp_ids, formatted):
            chrom, pos = idx.split("_", 1)
            ref = self.ref_dict.get(idx, ".")
            alt = self.alt_dict.get(idx, ".")
            line = [chrom, pos, idx, ref, alt, ".", "PASS", ".", "GT"] + genos
            self.vcf_data.append("\t".join(map(str, line)))
        self.logger.info("Successfully formatted genotype data!")

    def write_vcf(self, filename):
        vcf_dir = self.output_dir / "vcfs"
        vcf_dir.mkdir(exist_ok=True, parents=True)
        path = vcf_dir / filename
        with open(path, "w") as fh:
            for h in self.header:
                fh.write(h + "\n")
            for l in self.vcf_data:
                fh.write(l + "\n")

    def is_gzip_file(self, filepath):
        with open(filepath, "rb") as f:
            return f.read(2) == b"\x1f\x8b"

    def read_vcf(self, filepath):
        try:
            self.ensure_tabix_indexed(filepath)
        except OSError:
            if self.is_gzip_file(filepath):
                newfp = filepath.with_suffix(filepath.suffix + ".ungzip")
                with gzip.open(filepath, "rb") as fin, open(newfp, "wb") as fout:
                    shutil.copyfileobj(fin, fout)
                filepath = newfp
            else:
                raise
        if ".gz" not in "".join(filepath.suffixes):
            filepath = filepath.with_suffix(filepath.suffix + ".gz")
        return pysam.VariantFile(str(filepath)), filepath

    def sort_vcf_file(self, filepath):
        sorted_fp = filepath.with_suffix(".sorted.vcf.gz")
        with pysam.VariantFile(str(filepath)) as vin, pysam.BGZFile(str(sorted_fp), "w") as vout:
            vout.write(str(vin.header).encode())
            for r in sorted(vin.fetch(), key=lambda x: (x.contig, x.start)):
                vout.write(str(r).encode())
        return sorted_fp

    def ensure_tabix_indexed(self, filepath):
        tbi = filepath.with_suffix(filepath.suffix + ".tbi")
        if not tbi.exists():
            pysam.tabix_index(str(filepath), preset="vcf", force=True)
        return filepath

    def preprocess_vcf_files(self, input_vcf_fn, is_final_output=False):
        if not is_final_output:
            input_vcf = self.output_dir / "vcfs" / input_vcf_fn
        else:
            input_vcf = input_vcf_fn
        sorted_fp = self.sort_vcf_file(input_vcf)
        return self.ensure_tabix_indexed(sorted_fp)

    def subset_vcf_by_locus_ids(self, vcf_path, output_filename="subset.vcf"):
        v1, _ = self.read_vcf(Path(vcf_path))
        outpth = self.output_dir / "vcfs" / output_filename
        vout = pysam.VariantFile(str(outpth), "w", header=v1.header)
        for r in v1.fetch():
            lid = f"{r.chrom}_{r.pos}"
            if lid in self.snp_ids:
                try:
                    vout.write(r)
                except OSError:
                    pass
        vout.close()

    def merge_vcf_files_with_precedence(self, v1_path, v2_path, output_fn, prefer_v1=False):
        v1, _ = self.read_vcf(Path(v1_path))
        v2, _ = self.read_vcf(Path(v2_path))
        combined = sorted(set(v1.header.samples) | set(v2.header.samples))
        hdr = v1.header.copy()
        for s in combined:
            if s not in hdr.samples:
                hdr.add_sample(s)
        outpth = self.output_dir / "vcfs" / output_fn
        vout = pysam.VariantFile(str(outpth), "w", header=hdr)
        v2_recs = {f"{r.chrom}_{r.pos}": r for r in v2.fetch()}
        for r1 in v1.fetch():
            key = f"{r1.chrom}_{r1.pos}"
            r2 = v2_recs.get(key)
            flipped = False
            if r2 and self.check_allele_flip(r1, r2):
                flipped = True
                r2.ref, r2.alts = r2.alts[0], (r2.ref,)
            nr = vout.new_record()
            nr.chrom, nr.pos = r1.chrom, r1.pos
            nr.ref, nr.alts = r1.ref, r1.alts
            src = r2 if r2 and not prefer_v1 else r1
            for f in src.filter.keys():
                nr.filter.add(f)
            for k, v in src.info.items():
                nr.info[k] = v
            for s in combined:
                if s in v1.header.samples:
                    gt1 = r1.samples[s]["GT"]
                else:
                    gt1 = (None, None)
                if r2 and s in v2.header.samples:
                    gt2 = r2.samples[s]["GT"]
                else:
                    gt2 = (None, None)
                if flipped:
                    gt2 = tuple(1 - a if a is not None else a for a in gt2)
                chosen = gt1 if (not prefer_v1 and gt1 != (None, None)) else gt2
                nr.samples[s]["GT"] = chosen
            vout.write(nr)
        vout.close()
        return self.preprocess_vcf_files(output_fn)

    def replace_missing_genotypes(self, input_vcf_path, output_vcf_path):
        opener = gzip.open if ".gz" in "".join(Path(input_vcf_path).suffixes) else open
        with opener(input_vcf_path, "rt") as orig, open(output_vcf_path, "w") as out:
            for line in orig:
                if line.startswith("#"):
                    out.write(line)
                else:
                    cols = line.rstrip("\n").split("\t")
                    fmt_idx = cols[8].split(":").index("GT") if "GT" in cols[8] else None
                    if fmt_idx is not None:
                        for i in range(9, len(cols)):
                            parts = cols[i].split(":")
                            if parts[fmt_idx] == ".":
                                parts[fmt_idx] = "./."
                            cols[i] = ":".join(parts)
                    out.write("\t".join(cols) + "\n")

    def check_allele_flip(self, r1, r2):
        a1 = (r1.ref,) + tuple(r1.alts[:1] if r1.alts else (None,))
        a2 = (r2.ref,) + tuple(r2.alts[:1] if r2.alts else (None,))
        return a1 == tuple(reversed(a2))

    def find_shared_loci(self, v1_path, v2_path):
        v1, _ = self.read_vcf(Path(v1_path))
        v2, _ = self.read_vcf(Path(v2_path))
        loci1 = {f"{r.chrom}-{r.pos}" for r in v1.fetch()}
        loci2 = {f"{r.chrom}-{r.pos}" for r in v2.fetch()}
        return [x.split("-") for x in loci1 & loci2]

    def are_genotypes_concordant(self, gt1, gt2, flipped):
        if flipped:
            gt2 = tuple(1 - a if a is not None else a for a in gt2)
        g1 = sorted(str(a) if a is not None else "." for a in gt1)
        g2 = sorted(str(a) if a is not None else "." for a in gt2)
        return g1 == g2

    def calculate_overall_concordance(self, vcf1_path, vcf2_path, output_dir):
        self.logger.info("Calculating sample and locus concordance...")

        vcf1, vcf1_path = self.read_vcf(vcf1_path)
        vcf2, vcf2_path = self.read_vcf(vcf2_path)

        shared_loci = self.find_shared_loci(vcf1_path, vcf2_path)

        sample_concordance = defaultdict(
            lambda: {"concordant": 0, "discordant": 0, "missing": 0}
        )

        # Get common samples between the two VCF files
        common_samples = set(vcf1.header.samples) & set(vcf2.header.samples)

        for chrom, pos in shared_loci:
            pos = int(pos)
            record1 = next(vcf1.fetch(chrom, pos - 1, pos), None)
            record2 = next(vcf2.fetch(chrom, pos - 1, pos), None)

            if record1 is None or record2 is None:
                continue

            # Check and correct for allele flips
            ref_flipped = self.check_allele_flip(record1, record2)

            if ref_flipped and record2.alts is not None:
                # Swap ref and alt in record2 for comparison purposes
                record2.ref, record2.alts = record2.alts[0], (record2.ref,)

            for sample in common_samples:
                gt1 = record1.samples[sample]["GT"]
                gt2 = record2.samples[sample]["GT"]

                if None in gt1 or None in gt2:
                    sample_concordance[sample]["missing"] += 1
                elif self.are_genotypes_concordant(gt1, gt2, ref_flipped):
                    sample_concordance[sample]["concordant"] += 1
                else:
                    sample_concordance[sample]["discordant"] += 1

                    if self.debug:
                        self.logger.debug(
                            f"Discordant sample: {sample}, GT1: {gt1}, GT2: {gt2}, flipped: {ref_flipped}"
                        )

        sample_summary = pd.DataFrame(sample_concordance).T

        tabdir = output_dir / "tables"
        tabdir.mkdir(exist_ok=True, parents=True)

        plotdir = output_dir / "plots"
        plotdir.mkdir(exist_ok=True, parents=True)

        sample_summary["total"] = (
            sample_summary["concordant"] + sample_summary["discordant"]
        )
        sample_summary["concordant_pct"] = (
            sample_summary["concordant"] / sample_summary["total"]
        ) * 100
        sample_summary["discordant_pct"] = (
            sample_summary["discordant"] / sample_summary["total"]
        ) * 100
        sample_summary["missing_pct"] = (
            sample_summary["missing"]
            / (sample_summary["total"] + sample_summary["missing"])
        ) * 100

        ssumfn = tabdir / f"{self.prefix}_sample_concordance_summary.csv"
        sample_summary.to_csv(ssumfn, header=True, index=False)

        concordant_mean = sample_summary["concordant_pct"].mean()
        concordant_median = sample_summary["concordant_pct"].median()
        concordant_std = sample_summary["concordant_pct"].std()

        discordant_mean = sample_summary["discordant_pct"].mean()
        discordant_median = sample_summary["discordant_pct"].median()
        discordant_std = sample_summary["discordant_pct"].std()

        missing_mean = sample_summary["missing_pct"].mean()
        missing_median = sample_summary["missing_pct"].median()
        missing_std = sample_summary["missing_pct"].std()

        summary_stats = {
            "Concordant Mean": concordant_mean,
            "Concordant Median": concordant_median,
            "Concordant StdDev": concordant_std,
            "Discordant Mean": discordant_mean,
            "Discordant Median": discordant_median,
            "Discordant StdDev": discordant_std,
            "Missing Mean": missing_mean,
            "Missing Median": missing_median,
            "Missing StdDev": missing_std,
        }

        ssfn = tabdir / f"{self.prefix}_concordance_summary_stats.csv"
        pd.DataFrame(summary_stats, index=[0]).to_csv(ssfn, index=False)

        plt.figure(figsize=(10, 6))

        dfmelt = sample_summary.melt(
            value_vars=["concordant_pct", "discordant_pct", "missing_pct"],
            var_name="Variable",
            value_name="Value",
            id_vars=["total", "concordant", "discordant", "missing"],
        )

        # Remove any with all missing data.
        dfmelt = dfmelt.dropna(subset=["Value"], how="any")

        # Make the tick labels clearer.
        dfmelt["Variable"] = dfmelt["Variable"].map(
            {
                "concordant_pct": "Concordant",
                "missing_pct": "Missing",
                "discordant_pct": "Discordant",
            }
        )

        # Plot violin plots with boxplots overlaid.
        sns.violinplot(
            data=dfmelt,
            x="Variable",
            y="Value",
            hue="Variable",
            hue_order=["Concordant", "Discordant", "Missing"],
            palette="turbo",
            density_norm="count",
            inner=None,
            linewidth=0,
            saturation=0.4,
        )

        sns.boxplot(
            data=dfmelt,
            x="Variable",
            y="Value",
            hue="Variable",
            palette="turbo",
            boxprops={"zorder": 2},
            hue_order=["Concordant", "Discordant", "Missing"],
            width=0.2,
        )
        plt.ylim([0, 110])
        plt.title("Per-sample Concordance, Discordance, and Missing Data")
        plt.xlabel("Genotype Agreement")
        plt.ylabel("Percentage")

        fn = f"{self.prefix}_concordance_discordance_boxplot.{self.plt_ft}"
        pltfn = plotdir / fn
        plt.savefig(pltfn, bbox_inches="tight", facecolor="white")
        plt.close()

        plt.figure(figsize=(10, 6))

        sns.regplot(
            data=sample_summary, x="discordant_pct", y="missing_pct", color="darkorchid"
        )

        dftmp = sample_summary[["discordant_pct", "missing_pct"]]
        dftmp = dftmp.dropna(axis=0, how="any")

        corr_coef, pval = pearsonr(x=dftmp["discordant_pct"], y=dftmp["missing_pct"])

        pval = "P < 0.0001" if pval < 0.0001 else f"P = {round(pval, 2)}"
        rval = round(corr_coef, 2)
        plottext = f"R = {rval} ({pval})"

        plt.annotate(plottext, xy=(0.85, 0.95), xycoords="axes fraction")
        plt.ylim([0, 110])
        plt.title("Regression between Discordance and Missing Data Percent")
        plt.xlabel("Discordance (%)")
        plt.ylabel("Sample-wise Missing Data (%)")

        pltfn = plotdir / f"{self.prefix}_discordance_missing_regression.{self.plt_ft}"
        plt.savefig(pltfn, bbox_inches="tight", facecolor="white")
        plt.close()

        self.logger.info("Successfully calculated concordances!")
        self.logger.info(f"Saved final plot to: {str(pltfn)}")


def get_args():
    parser = argparse.ArgumentParser(
        description="Convert GTseq CSV → VCF, optionally merge with ddRADseq VCF and plot discordance.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--gtseq", type=str, required=True, help="Path to the input GTseq CSV file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory for output files.")
    parser.add_argument("--prefix", type=str, required=True, help="Prefix for output files.")
    parser.add_argument("--radseq", type=str, required=False, default=None, help="Path to the input ddRADseq VCF file.")
    parser.add_argument("--plot_filetype", type=str, default="pdf", help="Filetype for output plots.")
    parser.add_argument("--str2drop", type=str, nargs="*", default=[], help="Strings to drop from GTseq columns.")
    return parser.parse_args()

def main():
    args = get_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    converter = GTseqToVCF(
        filepath=args.gtseq,
        output_dir=output_dir,
        prefix=args.prefix,
        plot_filetype=args.plot_filetype,
        str2drop=args.str2drop,
    )

    converter.load_data()
    converter.parse_sample_column()
    converter.calculate_ref_alt_alleles()
    converter.transpose_data()
    converter.create_vcf_header()
    converter.format_genotype_data()
    gtseq_vcf = f"{converter.prefix}_gtseq.vcf"
    converter.write_vcf(gtseq_vcf)
    gtseq_vcf_path = converter.preprocess_vcf_files(Path(gtseq_vcf))

    if args.radseq:
        subset_vcf = f"{converter.prefix}_radseq_subset.vcf"
        converter.subset_vcf_by_locus_ids(args.radseq, output_filename=subset_vcf)
        subset_path = converter.preprocess_vcf_files(Path(subset_vcf))
        merged_tmp = f"{converter.prefix}_merged.tmp.vcf"
        merged_path = converter.merge_vcf_files_with_precedence(
            subset_path, gtseq_vcf_path, merged_tmp, prefer_v1=False
        )
        final_name = merged_path.stem.replace(".tmp","") + merged_path.suffix
        merged_dir = output_dir / "merged_vcf"
        merged_dir.mkdir(exist_ok=True)
        final_vcf = merged_dir / final_name
        converter.replace_missing_genotypes(merged_path, final_vcf)
        final_vcf_path = converter.preprocess_vcf_files(final_vcf, is_final_output=True)
    else:
        converter.logger.info("No --radseq provided; skipping merge/concordance.")
        final_vcf_path = gtseq_vcf_path

    for tmp in output_dir.glob("*.tmp.vcf"):
        tmp.unlink()

    if not Path(final_vcf_path).exists():
        converter.logger.error(f"Final VCF not found: {final_vcf_path}")
        raise FileNotFoundError(f"{final_vcf_path} missing")

    converter.logger.info(f"Final VCF → {final_vcf_path}")
    converter.logger.info("All done.")

if __name__ == "__main__":
    main()
