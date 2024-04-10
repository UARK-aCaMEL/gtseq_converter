# GTSeq2VCF

A tool to convert GTseq data to VCF format and merge with ddRADseq VCF file.  

## Installation

Navigate to the project's parent directory, and enter the command:

```pip install .```

## Usage

You can see the usage by typing:

```python gtseq2vcf/gtseq2vcf.py -h```

### Excluding Loci with Patterns

To exclude GT-Seq loci using a pattern, you just need to specify the ```--str2drop``` argument when running the script, followed by any number of patterns to exclude, separated by a space (no quotes surrounding). For example:  

```
python gtseq2vcf/gtseq2vcf.py --gtseq data/my_gtseq_data.csv --radseq data/my_radseq_data.vcf.gz --output_dir analysis --str2drop Ovi_ _PRNP_
```

## What does the script do?

The script takes two VCF files as input:

1. VCF file derived from double digest RAD sequencing (ddRADseq) or RAD sequencing (RADseq). The VCF file is expected to have been generated with ipyrad, but other formats *might* work.  

2. CSV (comma-delimited) file derived from GT-Seq panel.

The script then merges the two VCF files on the GT-Seq panel's loci and combines the samples from the two VCF files. If there are overlapping samples (duplicates), it keeps only the sample with the least missing genotypes.

You need to specify an output directory for the merged output files to be written to, and also it will make a nice visualization of discordance for the REF and ALT alleles between the two VCF files at the merged loci.  

