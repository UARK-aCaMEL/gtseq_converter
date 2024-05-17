# GTSeq2VCF

A tool to convert GTseq CSV data to VCF format and merge with an existing ddRADseq VCF file. Only the loci present in the GT-seq panel are retained, and duplicate samples that overlap are reduced to the duplicate with the least missing genotypes.  

Some nice visualizations and tables are generated. The tables are standard CSV files and can be easily re-loaded into Python or R if desired.

## Installation

Navigate to the project's parent directory, and enter the command:

```pip install .```

Eventually, we will upload this package to PyPi so that it can be installed easier, but for now, just install as above. It will install all necessary Python dependencies as well as the package.  

## Usage

You can see the usage by typing:

```python gtseq2vcf/gtseq2vcf.py -h```

### Command-line Script Arguments

There are several required options:

1. --gtseq - Path to the GT-seq CSV file, as returned from the core facilty.  
2. --radseq - Path to the VCF file to merge the GT-seq data with. Should contain only SNP data.
3. --output_dir - Path to the output directory where outputs will be saved.  
4. --prefix - String denoting the prefix to use for all output files.  

There are also two optional argument:  

+ --str2drop - Patterns of loci to remove from GT-seq data before merging with VCF file. This allows us to remove loci that do not contain numeric position information. The pattern can be any portion of a string contained within the locus names (header line of the GT-seq CSV file), and any loci names containing the pattern are removed before the merge. Finally, if you want to use multiple patterns, just separate the patterns by a space on the command-line. No need to use quotation marks surrounding the strings.  

+ --plot_filetype - Filetype to use for plots. Defaults to "pdf", but you can also use "png" or "jpg" instead if you would like to. The plots are then saved to the specified format. Technically you can use any file type supported by matplotlib, but to-date, we have only tested the "pdf", "png", and "jpg" options.  

### Excluding Loci with Patterns

To exclude GT-Seq loci using a pattern, you just need to specify the ```--str2drop``` argument when running the script, followed by any number of patterns to exclude, separated by a space (no quotes surrounding). For example:  

```
python gtseq2vcf/gtseq2vcf.py --gtseq data/my_gtseq_data.csv --radseq data/my_radseq_data.vcf.gz --output_dir analysis --prefix test --str2drop Ovi_ _PRNP_
```

**NOTE**: If you get an error that says the "POS" index contains non-numeric values, you probably need to specify loci to remove that do not conform to the standard pattern in the GT-seq CSV file header using the ```--str2drop``` option. Here is an example of how to use it:  

```--str2drop locus_pattern1 locus_pattern2 ...```   

In the ```--str2drop``` option, patterns contained with the offending locus names can be used, and the strings do not have to be the entire locus name. If there is more than one locus to remove by pattern, then you can separate them with a space on the command-line, as in the above code block example.  

## What does the script do?

The script takes two VCF files as input:  

1. VCF file derived from double digest RAD sequencing (ddRADseq) or RAD sequencing (RADseq). The VCF file is expected to have been generated with ipyrad, but other formats *might* work.   

2. CSV (comma-delimited) file derived from GT-Seq panel.  

The script then merges the two VCF files on the GT-Seq panel's loci and combines the samples from the two VCF files. If there are overlapping samples (duplicates), it keeps only the sample with the least missing genotypes.   

You need to specify an output directory for the merged output files to be written to, and also it will make a nice visualization of per-sample concordance, discordance, and missing data percentages. It also saves all the data to CSV files, which can be found in the ```<output_dir>/tables``` directory.    

Visualizations can be found in the ```<output_dir>/plots``` directory, and the final, merged output can be found in the ```<output_dir>/merged_vcf``` directory. Intermediate VCF files can be found in ```<output_dir>/vcfs```.  

Finally, you must specify a prefix with the ```--prefix <my_prefix>``` option. This allows you to save outputs from multiple analyses into one output directory.  

