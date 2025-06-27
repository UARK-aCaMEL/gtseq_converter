# GTSeq2VCF

A command-line tool to convert GTseq CSV data into VCF format and optionally merge with an existing ddRADseq VCF. Only panel loci are retained, and duplicate samples are reduced to the one with the fewest missing genotypes.

## Installation

Install from the project root:

```bash
pip install .
```

## Usage

Print the help message:

```bash
python gtseq2vcf/gtseq2vcf.py -h
```

### Required Arguments

* `--gtseq`: Path to the GT-seq CSV file.
* `--output_dir`: Directory where output files will be saved.
* `--prefix`: Prefix string for all output filenames.

### Optional Arguments

* `--radseq`: Path to a ddRADseq VCF file to merge with the GTseq-derived VCF.
* `--str2drop`: One or more string patterns; any locus names containing these patterns will be excluded before conversion. Provide multiple patterns separated by spaces.

### Example

Exclude non-standard loci and merge with a RADseq VCF:

```bash
python gtseq2vcf/gtseq2vcf.py \
  --gtseq data/my_gtseq_data.csv \
  --radseq data/my_radseq.vcf.gz \
  --output_dir results \
  --prefix analysis1 \
  --str2drop Ovi_ _PRNP_
```

## Functionality

1. **Convert GTseq CSV → VCF**
   Infers REF/ALT alleles, encodes genotypes, and writes a VCFv4.2 file.

2. **Merge with RADseq (optional)**
   Subsets the RADseq VCF to panel loci, merges samples, and resolves duplicates by keeping the sample with the fewest missing genotypes.

## Output Structure

* `<output_dir>/vcfs` – generated and intermediate VCF files
* `<output_dir>/merged_vcf` – final merged VCF when `--radseq` is used

Use the `<prefix>` to namespace multiple analyses within the same output directory.
