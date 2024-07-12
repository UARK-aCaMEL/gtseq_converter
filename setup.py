from setuptools import find_packages, setup

setup(
    name="GTseq2VCF",
    version="1.0.1",
    author="Bradley T. Martin, Ph.D.",
    author_email="evobio721@gmail.com",
    description="A tool to convert GTseq data to VCF format and merge with an existing ddRADseq VCF file, keeping only the loci in the GT-seq panel.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/btmartin721/GTseq2VCF",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
        "seaborn",
        "pysam",
        "scipy",
    ],
    python_requires=">=3.11",
    entry_points={
        "console_scripts": [
            "gtseq2vcf=gtseq2vcf.gtseq2vcf:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="genomics,GTseq,GT-seq,VCF,panel,ddRAD,bioinformatics,merge,genetics,population genetics",
)
