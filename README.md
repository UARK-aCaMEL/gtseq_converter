```markdown
# GTSeq2VCF

A command-line tool to convert GTseq CSV data into VCF format and optionally merge with an existing ddRADseq VCF. Only the loci present in the GT-seq panel are retained, and duplicate samples that overlap are reduced to the duplicate with the fewest missing genotypes.

## Installation

Install from the project root using pip:

```

pip install .

```

This will install all required Python dependencies along with the `gtseq2vcf` package.

## Usage

Print the help message:

```

python gtseq2vcf/gtseq2vcf.py -h

```

### Required Arguments

- `--gtseq` : Path to the GT-seq CSV file (as returned by the core facility).
- `--output_dir` : Directory where output files will be saved.
- `--prefix` : Prefix string for all output filenames.

### Optional Arguments

- `--radseq` : Path to an existing ddRADseq VCF file to merge with the GTseq-derived VCF.
- `--str2drop` : One or more string patterns; any locus names containing these patterns will be excluded before conversion. Provide multiple patterns separated by spaces without quotes.

### Example

Exclude non-standard loci by pattern and merge with a ddRADseq VCF:

```

python gtseq2vcf/gtseq2vcf.py&#x20;
\--gtseq data/my\_gtseq\_data.csv&#x20;
\--radseq data/my\_radseq.vcf.gz&#x20;
\--output\_dir results&#x20;
\--prefix analysis1&#x20;
\--str2drop Ovi\_ *PRNP*

```

## Functionality

1. **GTseq CSV → VCF**  
   Converts GT-seq genotype CSV into a standard VCFv4.2 file, inferring reference and alternate alleles and encoding genotypes appropriately.

2. **Optional RADseq Merge**  
   If a RADseq VCF is provided, retains only panel loci, subsets the RADseq VCF, and merges samples. In cases of duplicate sample IDs, the sample with the fewest missing genotypes is kept.

3. **Output Structure**  
   - `<output_dir>/vcfs`: Contains generated and intermediate VCF files.  
   - `<output_dir>/merged_vcf`: Contains the final merged VCF if `--radseq` is used.  

Use the prefix to organize multiple analyses within the same output directory.
```
