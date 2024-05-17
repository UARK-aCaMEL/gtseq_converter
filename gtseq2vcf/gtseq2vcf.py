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
    logger = logging.getLogger()  # Root logger
    logger.setLevel(logging.INFO)

    # Clear existing handlers
    logger.handlers = []

    Path(log_file).parents[0].mkdir(parents=True, exist_ok=True)

    # Create handlers
    file_handler = logging.FileHandler(log_file)
    stream_handler = logging.StreamHandler()

    # Create formatters and add it to handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    )
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)


class GTseqToVCF:
    def __init__(
        self, filepath, output_dir, prefix, plot_filetype, str2drop=None, debug=False
    ):
        self.filepath = filepath
        self.prefix = prefix
        self.plt_ft = plot_filetype
        self.debug = debug
        self.data = None
        self.header = []
        self.vcf_data = []
        self.snp_ids = None
        self.sample_ids = None
        self.header_key = "Sample"

        # Do some handling and preprocessing to prevent errors downstream.
        if self.plt_ft.startswith("."):
            self.plt_ft = self.plt_ft.lstrip(".")
        self.plt_ft = self.plt_ft.lower()

        if not isinstance(self.prefix, str):
            self.prefix = str(self.prefix)

        # Set up logfile.
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        logpath = output_dir / f"{self.prefix}_logfile.txt"

        setup_logger(logpath)
        self.logger = logging.getLogger(__name__)

        self.metadata_cols = [
            "Raw Reads",
            "On-Target Reads",
            "%On-Target",
            "%GT",
            "IFI",
        ]

        self.str2drop = str2drop

    def load_data(self):
        self.logger.info("Loading GT-seq CSV file...")

        # Load the data
        self.data = pd.read_csv(self.filepath, header=0)

        # Remove metadata columns and store in self.metadata
        self.metadata = pd.concat(
            [self.data.pop(x) for x in self.metadata_cols if x in self.metadata_cols],
            axis=1,
        )

        self.data.columns = self.data.columns.str.strip()

        if self.str2drop is not None:
            if isinstance(self.str2drop, list):
                for to_drop in self.str2drop:
                    self.data = self.data.loc[
                        :,
                        ~self.data.columns.str.lower().str.startswith(to_drop.lower()),
                    ]
                    self.data = self.data.loc[
                        :, ~self.data.columns.str.lower().str.contains(to_drop.lower())
                    ]
            elif isinstance(self.str2drop, str):
                self.data = self.data.loc[
                    :,
                    ~self.data.columns.str.lower().str.startswith(
                        self.str2drop.lower()
                    ),
                ]
                self.data = self.data.loc[
                    :,
                    ~self.data.columns.str.lower().str.contains(self.str2drop.lower()),
                ]

            else:
                msg = f"'str2drop' must be a list of strings, but got: {type(self.str2drop)}"
                self.logger.error(msg)
                raise TypeError(msg)

        self.logger.info("Successfully loaded GT-seq CSV file!")

    def parse_sample_column(self):
        """
        Parses the column headers into separate 'CHROM' and 'POS' columns, and reorders the DataFrame to have 'CHROM' and 'POS' as the first two columns.
        Validates the 'POS' to be strictly numeric and the genotype columns to contain valid allele pairs. Assumes 'Sample' is the index of the DataFrame.

        Returns:
            None
        """
        self.logger.info("Parsing loci and sample data from GT-seq file..")

        # Extract 'CHROM' and 'POS' from the column headers\
        chrom_pos_split = [col.split("_", 1) for col in self.data.columns]
        chrom = [x[0] for x in chrom_pos_split if "Sample" not in x]
        pos = [x[1] for x in chrom_pos_split if "Sample" not in x]

        self.chrom = chrom
        self.pos = pos

        # Now, set the first column as index and transpose the DataFrame
        self.data = self.data.set_index(self.data.columns[0]).T
        self.data.insert(0, "CHROM", chrom)
        self.data.insert(1, "POS", pos)

        # After transposition, create a MultiIndex from the first two rows (which contain 'CHROM' and 'POS')
        self.data.index = pd.MultiIndex.from_arrays(
            [self.data.iloc[:, 0], self.data.iloc[:, 1]], names=("CHROM", "POS")
        )

        # Drop the now redundant first two rows
        self.data = self.data.drop(["CHROM", "POS"], axis=1)

        # Validate the 'POS' to ensure all values are numeric
        if not self.data.index.get_level_values("POS").str.isnumeric().all():
            msg = "Non-numeric values found in 'POS' index level. Try using the '--str2drop' option to remove offending loci."
            self.logger.error(msg)
            raise ValueError(msg)

        # Define valid allele pairs
        valid_alleles = {
            "AA",
            "TT",
            "CC",
            "GG",
            "AT",
            "TA",
            "AC",
            "CA",
            "AG",
            "GA",
            "TC",
            "CT",
            "TG",
            "GT",
            "CG",
            "GC",
            "--",
            "NN",
            "00",
        }

        # Validate the genotype columns to ensure they contain valid allele
        # pairs
        for sample in self.data:
            if (
                not self.data[sample]
                .apply(lambda x: x in valid_alleles or x == "00" or x == "0")
                .all()
            ):
                msg = f"Invalid alleles found in sample {sample}: {self.data[sample].tolist()}"
                self.logger.error(msg)
                raise ValueError(msg)

        self.data = self.data.T

        if not self.data.empty:
            self.logger.info("Successfully parsed loci and sample data!")
        else:
            msg = "GT-seq data parsing was unsuccessful. Yielded empty pandas DataFrame. Please check the file format."
            self.logger.error(msg)
            raise pd.errors.EmptyDataError(msg)

    def calculate_ref_alt_alleles(self):
        """
        Calculates the reference (REF) and alternate (ALT) alleles for each SNP across all samples. Assumes the DataFrame has 'CHROM' and 'POS' as a MultiIndex for the columns and the rest of the columns are SNP genotypes.
        """

        self.logger.info("Calculating reference and alternate alleles...")

        # Define a helper function to calculate REF and ALT alleles for a
        # series of genotypes
        def get_ref_alt(genotypes):
            # Flatten the genotype pairs, split by any non-allele character,
            # and count occurrences, excluding missing data

            allele_counts = Counter(
                allele
                for genotype in genotypes
                for allele in genotype
                if allele not in {"-", "N", "0", "00"}
            )

            if allele_counts:
                # Get the most common alleles
                common_alleles = allele_counts.most_common()

                ref_allele = common_alleles[0][0] if common_alleles else None
                alt_allele = common_alleles[1][0] if len(common_alleles) > 1 else "."
            else:
                # No valid alleles present
                ref_allele, alt_allele = None, "."

            return ref_allele, alt_allele

        # Apply the helper function to each SNP column to calculate REF and ALT alleles
        ref_alt_alleles = self.data.apply(get_ref_alt, axis=0)

        refs = ref_alt_alleles.iloc[0]
        alts = ref_alt_alleles.iloc[1]

        self.data = self.data.T

        # Assign the REF and ALT alleles to the corresponding MultiIndex levels
        self.data["REF"] = refs
        self.data["ALT"] = alts

        self.logger.info("Successfully assigned REF and ALT alleles!")

    def join_multiindex(self, df, separator="_"):
        """
        Joins the values of a Pandas MultiIndex into one string per row.

        Args:
        df: A pandas DataFrame with a MultiIndex.separator: A string separator used to join MultiIndex values.

        Returns:
        A pandas Series with the joined strings for each row.
        """
        # Ensure the DataFrame has a MultiIndex
        if not isinstance(df.index, pd.MultiIndex):
            raise ValueError("The DataFrame does not have a MultiIndex.")

        # Join the MultiIndex levels with the specified separator
        return df.index.to_series().apply(lambda x: separator.join(map(str, x)))

    def transpose_data(self):
        # Transpose the dataframe, so loci are rows and samples are columns
        self.data = self.data.T.reset_index()

        # Rename columns to reflect that the first column is now 'sample'
        self.data = self.data.rename(columns={"index": self.header_key})

        # Extract SNP IDs from the column names
        self.snp_ids = self.data.columns[1:].tolist()
        self.snp_ids = ["_".join(x) for x in self.snp_ids]

        # Extract sample IDs which are now the column names, starting from the first sample column
        self.sample_ids = self.data[self.header_key].tolist()
        self.sample_ids = [x for x in self.sample_ids if x not in ["REF", "ALT"]]

    def create_vcf_header(self):
        # VCF header lines start with '##' and the column header line starts with '#'

        self.logger.info("Creating VCF header...")

        self.sample_ids = [
            x.split("_")[-1] for x in self.sample_ids if x.startswith("GTseek")
        ]

        if not self.sample_ids:
            msg = "Error parsing sampleIDs from GT-seq file."
            self.logger.error(msg)
            raise ValueError(msg)

        sample_ids = "\t".join(self.sample_ids)

        self.header = [
            "##fileformat=VCFv4.2",
            "##source=GTseqToVCFConverter",
            "##reference=GenomeRef",
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + sample_ids.strip("\n"),
        ]

        self.logger.info("Successfully created VCF header!")

    def format_genotype_data(self):
        """Formats the genotype data for each SNP according to VCF specifications."""

        self.logger.info("Formatting genotype data...")

        self.data = self.data.T
        self.data = self.data.iloc[1:, :]
        self.data["ID"] = self.join_multiindex(self.data)
        self.data = self.data.reset_index()

        cols = ["CHROM", "POS"] + self.sample_ids + ["REF", "ALT", "ID"]
        self.data.columns = cols
        self.data = self.data.set_index("ID")

        # Ensure SNP IDs are sorted properly for consistent output
        self.snp_ids = sorted(
            self.snp_ids, key=lambda x: (x.split("_")[0], int(x.split("_")[1]))
        )

        # Replace genotype strings with VCF format genotypes
        for snp in self.snp_ids:
            ref_allele = self.data.at[snp, "REF"]
            alt_allele = self.data.at[snp, "ALT"]

            for sample in self.sample_ids:
                genotype = self.data.at[snp, sample]

                # Initialize a list to hold the formatted genotype indices
                formatted_genotype = []

                # Iterate over each character in the genotype string
                for allele in genotype:
                    if allele == ref_allele:
                        formatted_genotype.append("0")
                    elif allele == alt_allele or allele in "ATCG":
                        formatted_genotype.append("1")
                    else:
                        formatted_genotype.append(".")  # For missing alleles

                # Join the allele indices with '/' to form the VCF genotype
                vcf_genotype = "/".join(formatted_genotype)
                self.data.at[snp, sample] = vcf_genotype

        self.data = self.data.copy()
        self.data["ID"] = self.data.index.copy()

        # Format the data for VCF output
        for snp in self.snp_ids:
            chrom = self.data.at[snp, "CHROM"]
            pos = self.data.at[snp, "POS"]
            ref = self.data.at[snp, "REF"]
            alt = self.data.at[snp, "ALT"]
            loc_id = self.data.at[snp, "ID"]

            # Construct the VCF line for the SNP
            vcf_line = [chrom, pos, loc_id, ref, alt, ".", "PASS", ".", "GT"]

            # Append the genotype data for each sample
            genotypes = [self.data.at[snp, sample] for sample in self.sample_ids]
            vcf_line.extend(genotypes)

            # Join the line into a string and append to the VCF data
            self.vcf_data.append("\t".join(vcf_line))

        self.logger.info("Successfully formatted genotype data!")

    def check_allele_flip(self, record1, record2):
        """Determines if there is an allele flip between two VCF records using allele comparison.

        Args:
            record1 (pysam.VariantRecord): The first VCF record.
            record2 (pysam.VariantRecord): The second VCF record.

        Returns:
            bool: True if there is an allele flip, False otherwise.
        """
        if record1.alts is not None:
            if len(record1.alts) > 1:
                alts = (record1.alts[0],)
            else:
                alts = record1.alts
            alleles1 = (record1.ref,) + alts
        else:
            alleles1 = (record1.ref,) + (None,)
        if record2.alts is not None:
            if len(record2.alts) > 1:
                alts = (record2.alts[0],)
            else:
                alts = record2.alts
            alleles2 = (record2.ref,) + alts
        else:
            alleles2 = (record2.ref,) + (None,)

        # Check if alleles1 is a reverse of alleles2
        ref_flipped = alleles1 == tuple(reversed(alleles2))
        return ref_flipped

    def write_vcf(self, output_filename):
        vcf_dir = self.output_dir / "vcfs"
        vcf_dir.mkdir(exist_ok=True, parents=True)
        pth = vcf_dir / output_filename

        with open(pth, "w") as vcf_file:
            # Write the header lines
            for line in self.header:
                vcf_file.write(line + "\n")
            # Write the SNP data lines
            for line in self.vcf_data:
                vcf_file.write(line + "\n")

    def subset_by_locus_ids(self):
        """
        Subsets the data by the locus IDs which are the joined CHROM and POS columns.

        Raises:
            AttributeError: If snp_ids have not been set before this method is called.

        Returns:
            pandas.DataFrame: A subset of the original DataFrame containing only the rows with index matching the locus IDs in self.snp_ids.
        """
        self.logger.info("Subsetting loci to GT-seq panel...")

        # Ensure that the snp_ids attribute has been populated
        if self.snp_ids is None:
            msg = "snp_ids must be set before calling subset_by_locus_ids."
            self.logger.error(msg)
            raise AttributeError(msg)

        # Create a MultiIndex from the CHROM and POS columns if it's not already done
        if not isinstance(self.data.index, pd.MultiIndex):
            self.data = self.data.set_index(["CHROM", "POS"])

        # Find the intersection of the data's index with the snp_ids
        locus_ids = set(self.join_multiindex(self.data, separator="_"))
        subset_ids = locus_ids.intersection(self.snp_ids)

        # Subset the DataFrame based on the matching locus IDs
        subset_data = self.data.loc[subset_data.index.isin(subset_ids)]

        if subset_data.empty:
            msg = "Encountered an issue subsetting the VCF data to the GT-seq panel loci. Subsetted DataFrame was empty. Are the loci present in both the VCF and GT-seq files?"
            self.logger.error(msg)
            raise pd.errors.EmptyDataError(msg)

        self.logger.info("Successfully subset VCF loci to GT-seq panel!")

        return subset_data

    def subset_vcf_by_locus_ids(
        self, vcf_path, output_filename="phase6_gtseq_subset.vcf"
    ):
        """
        Subsets the VCF file by the locus IDs which are the joined CHROM and POS columns.

        Args:
            vcf_path (str): Path to the input VCF file.
            output_path (str): Directory path where the subsetted VCF file will be written.
            output_filename (str, optional): Name of the subsetted VCF file.Defaults to "phase6_gtseq_subset.vcf".

        Returns:
            None
        """
        vcf_in, vcf_path = self.read_vcf(vcf_path)

        outpth = self.output_dir / "vcfs" / output_filename.name

        # Create a new VCF file for the output
        vcf_out = pysam.VariantFile(outpth, "w", header=vcf_in.header)

        # Read through the VCF file and write only the records with locus IDs in self.snp_ids
        for record in vcf_in:
            locus_id = f"{record.chrom}_{record.pos}"
            if locus_id in self.snp_ids:
                try:
                    vcf_out.write(record)
                except OSError:
                    pass

        # Close the VCF files
        vcf_in.close()
        vcf_out.close()

    def is_gzip_file(self, filepath):
        with open(filepath, "rb") as test_f:
            # Test if is gzip file (not bgzip).
            return test_f.read(2) == b"\x1f\x8b"

    def read_vcf(self, filepath):
        """
        Open a VCF file using pysam and return the VariantFile object.

        Args:
            filepath (str): The file path to the VCF file.

        Returns:
            pysam.VariantFile: VariantFile object for the VCF.
        """
        # Make sure that VCF file is tabix indexed.
        try:
            self.ensure_tabix_indexed(filepath)
        except OSError as e:
            if self.is_gzip_file(filepath):
                self.logger.warning(
                    f"Input VCF file {str(filepath)} is GZipped, but BGZip or uncompressed format was expected. Attempting to gunzip."
                )
                shutil.copyfile(
                    filepath, filepath.with_suffix(str(filepath.suffix) + ".gzip")
                )

                filepath_orig = filepath
                filepath = filepath.with_name(f"unzip_{str(filepath.name) + '.gzip'}")

                with gzip.open(filepath_orig, "r") as f_in, open(
                    filepath, "wb"
                ) as f_out:
                    shutil.copyfileobj(f_in, f_out)
                warnings.warn(f"Copied {str(filepath)} to unzip_{str(filepath)}")

            elif "building of index" in str(e):
                try:
                    self.logger.warning(
                        "VCF file is unsorted. Attempting to sort the file and retry tabix indexing."
                    )
                    sorted_filepath = self.sort_vcf_file(filepath)
                    self.ensure_tabix_indexed(sorted_filepath)
                    filepath = sorted_filepath
                except OSError as e2:
                    msg = f"Unknown error encountered when indexing VCF file: {str(e2)}"
                    self.logger.error(msg)
                    raise OSError(msg)
            else:
                msg = f"Unknown error encountered when indexing VCF file: {str(e)}"
                self.logger.error(msg)
                raise OSError(msg)

        if not ".gz" in filepath.suffixes:
            filepath = filepath.with_suffix(filepath.suffix + ".gz")

        return pysam.VariantFile(filepath), filepath

    def sort_vcf_file(self, filepath):
        """
        Sort a VCF file.

        Args:
            filepath (str): The file path to the VCF file.

        Returns:
            Path: The file path to the sorted VCF file.
        """

        self.logger.warning(
            f"Saved sorted VCF file to: {str(filepath.with_suffix('.sorted.vcf.gz'))}"
        )
        sorted_filepath = filepath.with_suffix(".sorted.vcf.gz")

        with pysam.VariantFile(filepath, "r") as vcf_in, pysam.BGZFile(
            sorted_filepath, "w"
        ) as vcf_out:
            header = vcf_in.header
            vcf_out.write(str(header).encode())

            records = []
            for record in vcf_in:
                records.append(record)

            sorted_records = sorted(records, key=lambda r: (r.contig, r.start))

            for record in sorted_records:
                vcf_out.write(str(record).encode())

        return sorted_filepath

    def sort_and_write_vcf(self, input_path, output_fn):
        """
        Sorts a VCF file by chromosome and position and writes the sorted file to a new location.
        The input can be either a .vcf or .vcf.gz file, and the output will match the input format.

        Args:
            input_path (str): Path to the original unsorted VCF file.
            output_fn (str): Filename to write the sorted VCF file. Will be saved to self.output_dir.
        """
        vcf_in, input_path = self.read_vcf(input_path)
        vcf_dir = self.output_dir / "vcfs"
        vcf_dir.mkdir(exist_ok=True, parents=True)

        outpth = vcf_dir / output_fn.name

        # Create the output VariantFile object with the same header as the input
        with pysam.VariantFile(outpth, "w", header=vcf_in.header) as vcf_out:
            # Sort records by CHROM and POS before writing
            sorted_records = sorted(vcf_in, key=lambda x: (x.chrom, x.pos))
            for record in sorted_records:
                vcf_out.write(record)

    def ensure_tabix_indexed(self, filepath):
        """
        Checks if a .vcf file is indexed with Tabix; if not, indexes the file.
        If the file is compressed (.gz), it ensures the .tbi index exists. If not, it creates it.

        Notes:
            tabix_index() will automatically compress the VCF file and append ".gz" to the filename.

        Args:
            filepath (pathlib.Path): The path to the VCF file. Can be bgzipped or not.
        """
        if not Path(filepath.stem + ".tbi").exists():
            pysam.tabix_index(
                str(filepath), preset="vcf", force=True, keep_original=True
            )
            self.logger.info(f"Indexed VCF file: {filepath}")

        if ".gz" not in filepath.suffixes:
            return filepath.with_suffix(filepath.suffix + ".gz")
        return filepath

    def preprocess_vcf_files(self, input_vcf_fn, is_final_output=False):
        """Pre-processes input VCF files by sorting and indexing them.

        Returns:
            pathlib.Path: modified path with ".gz" extension (if file wasn't already compressed).
        """
        if not is_final_output:
            input_vcf = self.output_dir / "vcfs" / input_vcf_fn
        else:
            input_vcf = input_vcf_fn

        if not input_vcf.name.startswith("sorted_"):
            # Prepend prefix to sorted VCF filename.
            sorted_vcf = input_vcf.with_stem("sorted_" + input_vcf.stem)
        else:
            sorted_vcf = input_vcf

        # Read, sort, and index VCF file
        self.sort_and_write_vcf(input_path=input_vcf, output_fn=sorted_vcf)
        return self.ensure_tabix_indexed(filepath=sorted_vcf)

    def calculate_missing_percentage_for_loci(self, vcf_path):
        """
        Calculate the percentage of missing genotype data for each locus in a VCF file.

        Args:
            vcf_path (str): The file path to the VCF file.

        Returns:
            dict: A dictionary with locus identifiers (CHROM_POS) as keys and the percentage of missing data as values.
        """
        vcf, vcf_path = self.read_vcf(vcf_path)
        sample_count = len(vcf.header.samples)
        loci_missing_counts = {}

        for record in vcf:
            missing_count = sum(
                1
                for sample in record.samples
                if record.samples[sample]["GT"] == (None, None)
            )
            loci_missing_counts[f"{record.chrom}_{record.pos}"] = (
                missing_count / sample_count
            ) * 100

        return loci_missing_counts

    def merge_vcf_files_with_missing_data_handling(
        self, vcf1_path, vcf2_path, output_fn
    ):
        """
        Merge two VCF files by choosing for each locus the genotypes from the file with the least missing data. If a sample is only present in one file, its genotype is still included.

        Args:
            vcf1_path (pathlib.Path): Path to the first VCF file.
            vcf2_path (pathlib.Path): Path to the second VCF file.
            output_fn (pathlib.Path): Filename of the output VCF file. Will have ".gz" appended as extension if not already present, and will be sorted, compressed, and indexed.
        """
        self.logger.info("Merging VCF and GT-seq data...")

        vcf_dir = self.output_dir / "vcfs"
        vcf_dir.mkdir(exist_ok=True, parents=True)
        outpth = vcf_dir / output_fn

        # Calculate missing data percentages for each locus in both VCFs
        vcf1_missing = self.calculate_missing_percentage_for_loci(vcf1_path)
        vcf2_missing = self.calculate_missing_percentage_for_loci(vcf2_path)

        self.calculate_overall_concordance(vcf1_path, vcf2_path, self.output_dir)

        # Read in VCF files
        vcf1, vcf1_path = self.read_vcf(vcf1_path)
        vcf2, vcf2_path = self.read_vcf(vcf2_path)

        # Prepare output VCF with header from vcf1, including all samples from
        # both files
        combined_samples = sorted(set(vcf1.header.samples) | set(vcf2.header.samples))

        header = vcf1.header.copy()
        for sample in combined_samples:
            if sample not in vcf1.header.samples:
                header.add_sample(sample)
        output_vcf = pysam.VariantFile(outpth, "w", header=header)

        # Cache vcf2 records for efficiency
        vcf2_records = {f"{rec.chrom}_{rec.pos}": rec for rec in vcf2.fetch()}

        for record1 in vcf1.fetch():
            locus_key = f"{record1.chrom}_{record1.pos}"
            record2 = vcf2_records.get(locus_key, None)

            new_record = output_vcf.new_record()
            new_record.chrom = record1.chrom
            new_record.pos = record1.pos
            new_record.ref = record1.ref
            new_record.id = record1.id
            new_record.alts = record1.alts
            new_record.qual = record1.qual

            # Setting FILTER and INFO fields from the chosen record (prioritize record1 if equal or no record2)
            if record2 and vcf2_missing.get(locus_key, 0) < vcf1_missing.get(
                locus_key, float("inf")
            ):
                for filter_key in record2.filter.keys():
                    new_record.filter.add(filter_key)
                for key, value in record2.info.items():
                    new_record.info[key] = value
            else:
                for filter_key in record1.filter.keys():
                    new_record.filter.add(filter_key)
                for key, value in record1.info.items():
                    new_record.info[key] = value

            for sample in combined_samples:
                gt1 = (None, None)
                if sample in record1.header.samples:
                    gt1 = (
                        record1.samples[sample]["GT"]
                        if record1.samples[sample]["GT"] != (None, None)
                        else (None, None)
                    )

                gt2 = (None, None)
                if record2 and sample in record2.header.samples:
                    gt2 = (
                        record2.samples[sample]["GT"]
                        if record2.samples[sample]["GT"] != (None, None)
                        else (None, None)
                    )

                # Determine which genotype to use based on missing data
                if not record2 or (
                    gt1 != (None, None)
                    and vcf1_missing.get(locus_key, 0)
                    <= vcf2_missing.get(locus_key, float("inf"))
                ):
                    chosen_gt = gt1
                else:
                    chosen_gt = gt2 if gt2 != (None, None) else gt1

                # Apply the chosen genotype
                if chosen_gt != (
                    None,
                    None,
                ):  # Avoid setting explicitly missing genotypes if not necessary
                    new_record.samples[sample]["GT"] = chosen_gt
                else:
                    # Handle explicitly missing data if needed (e.g., set to
                    # './.' or leave as default missing)
                    pass

            # Write the new record after setting all sample genotypes
            output_vcf.write(new_record)

        output_vcf.close()

        self.logger.info("Successfully merged VCF and GT-seq data!")

        # Tabix index and sort output vcf file.
        return self.preprocess_vcf_files(outpth.name)

    def find_shared_loci(self, vcf1_path, vcf2_path):

        self.logger.info("Finding shared loci between VCF and GT-seq data...")

        vcf1, vcf1_path = self.read_vcf(vcf1_path)
        vcf2, vcf2_path = self.read_vcf(vcf2_path)

        # Extract loci from both files
        vcf1_loci = {f"{record.chrom}-{record.pos}" for record in vcf1.fetch()}
        vcf2_loci = {f"{record.chrom}-{record.pos}" for record in vcf2.fetch()}

        # Determine shared loci
        shared_loci = vcf1_loci.intersection(vcf2_loci)

        return [x.split("-") for x in shared_loci]

    def are_genotypes_concordant(self, gt1, gt2, ref_flipped):
        """
        Check if the genotypes are concordant, accounting for possible allele flips and unphased genotypes.
        """
        if ref_flipped:
            # Flip the alleles in gt2 if ref and alt alleles are flipped
            flipped_gt2 = []
            for allele in gt2:
                if allele == 0:
                    flipped_gt2.append(1)
                elif allele == 1:
                    flipped_gt2.append(0)
                else:
                    flipped_gt2.append(allele)
            gt2 = tuple(flipped_gt2)

        # Convert alleles to strings for comparison and handle None values
        gt1_str = [str(allele) if allele is not None else "." for allele in gt1]
        gt2_str = [str(allele) if allele is not None else "." for allele in gt2]

        # Sort alleles to handle unphased genotypes
        gt1_sorted = sorted(gt1_str)
        gt2_sorted = sorted(gt2_str)

        return gt1_sorted == gt2_sorted

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

    def replace_missing_genotypes(self, input_vcf_path, output_vcf_path):
        # Check if the input file is gzip-compressed based on its extension.
        if ".gz" in input_vcf_path.suffixes:
            file_open = gzip.open
        else:
            file_open = open

        # Use the determined open function to read the input file appropriately.
        with file_open(input_vcf_path, "rt") as original, open(
            output_vcf_path, "w"
        ) as updated:
            for line in original:
                if not line.startswith("#"):
                    # Process non-header lines.
                    entries = line.strip().split("\t")
                    format_idx = (
                        entries[8].split(":").index("GT")
                        if "GT" in entries[8]
                        else None
                    )

                    # Only proceed if the GT format is found.
                    if format_idx is not None:
                        for i in range(9, len(entries)):
                            sample_data = entries[i].split(":")
                            # Replace missing genotype values with './.'.
                            if sample_data[format_idx] == ".":
                                sample_data[format_idx] = "./."
                            entries[i] = ":".join(sample_data)
                    updated.write("\t".join(entries) + "\n")
                else:
                    # Directly write header lines.
                    updated.write(line)


def get_args():
    """
    Parses command line arguments for converting GTseq data to VCF format,merging with ddRADseq data, and generating discordance visualizations.

    Returns:
        Namespace: An argparse.Namespace object containing all the parsed command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="""This tool converts GTseq data into VCF format, merges it with ddRADseq VCF file, and handles missing data. It simplifies the bioinformatics pipeline for genomic researchers and provides a visualization of discordance for each of the shared loci.""",
        epilog="""Example usage:
        python script.py --gtseq mydata.csv --radseq radseqdata.vcf --output_dir ./processed
        
        This command processes 'mydata.csv', merges it with 'radseqdata.vcf', and outputs a merged VCF along with discordance visualization to the './processed' directory.
        
        """,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--gtseq",
        type=str,
        required=True,
        help="Path to the input GTseq CSV file. The file should contain genotyping data\n"
        "in CSV format, obtained from GTseq genotyping assays.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="The directory where the output files will be stored. This includes\n"
        "the merged VCF file and any discordance visualizations. The directory\n"
        "will be created if it does not exist.",
    )
    parser.add_argument(
        "--prefix", type=str, required=True, help="Prefix to use for output files."
    )
    parser.add_argument(
        "--radseq",
        type=str,
        required=True,
        help="Path to the input ddRADseq VCF file. This file should contain variant\n"
        "calling data from ddRADseq sequencing assays, formatted as a VCF file.",
    )
    parser.add_argument(
        "--plot_filetype",
        type=str,
        required=False,
        default="pdf",
        help="Filetype to use with output plots. Supported options: 'pdf', 'png', 'jpg'. Untested with other matplotlib options. Default: 'pdf'.",
    )
    parser.add_argument(
        "--str2drop",
        type=str,
        nargs="*",
        default=[],
        help="Optional. A list of strings to be excluded from the column names\n"
        "in the GTseq CSV file. Useful for removing unwanted data or metadata\n"
        "columns from the analysis.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    # Directory and file setup
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gtseq_output_vcf_fn = Path("gtseq_output.vcf")
    subset_vcf_fn = Path("radseq_subset.vcf")
    merged_vcf_fn = Path("merged_gtseq_radseq.tmp.vcf")

    # Initialize and process GTseq data
    converter = GTseqToVCF(
        filepath=args.gtseq,
        output_dir=output_dir,
        prefix=args.prefix,
        plot_filetype=args.plot_filetype,
        str2drop=args.str2drop,
    )

    converter.logger.info(f"Will save output files to: {str(output_dir)}")

    converter.load_data()
    converter.parse_sample_column()
    converter.calculate_ref_alt_alleles()
    converter.transpose_data()
    converter.create_vcf_header()
    converter.format_genotype_data()

    # Write the GTseq VCF
    converter.write_vcf(gtseq_output_vcf_fn)

    # Sort by CHROM and POS and tabix index file (if not already indexed)
    gtseq_output_vcf_path = converter.preprocess_vcf_files(gtseq_output_vcf_fn)

    # Subset RADseq VCF based on GTseq loci, write it, and index it
    converter.subset_vcf_by_locus_ids(Path(args.radseq), output_filename=subset_vcf_fn)

    subset_vcf_path = converter.preprocess_vcf_files(subset_vcf_fn)

    # Merge the GTseq VCF with the subsetted RADseq VCF
    merged_vcf_path = converter.merge_vcf_files_with_missing_data_handling(
        subset_vcf_path, gtseq_output_vcf_path, merged_vcf_fn
    )

    # Create the new filename by removing '.tmp' from the stem, then adding
    # back the rest of the desired suffix
    new_vcf_filename = merged_vcf_path.stem.replace(".tmp", "")

    # Use the .with_name() method to get the new path with corrected filename
    final_vcf_path = merged_vcf_path.with_name(new_vcf_filename)
    final_vcf_dir = converter.output_dir / "merged_vcf"
    final_vcf_dir.mkdir(exist_ok=True, parents=True)
    final_vcf_path = final_vcf_dir / final_vcf_path.name

    converter.replace_missing_genotypes(merged_vcf_path, final_vcf_path)
    converter.preprocess_vcf_files(final_vcf_path, is_final_output=True)

    # Clean up temporary / intermediate files.
    for pth in output_dir.iterdir():
        if ".tmp" in pth.suffixes:
            pth.unlink()

    if not final_vcf_path.exists() and not final_vcf_path.is_file():
        msg = f"Output VCF file {final_vcf_path} could not be found. The program may have encountered an unknown error."
        converter.logger.error(msg)
        raise FileNotFoundError(msg)

    converter.logger.info(f"The output VCF file is saved to: {final_vcf_path}")
    converter.logger.info("Script execution completed successfully!")
