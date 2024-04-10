from setuptools import find_packages, setup

setup(
    name="GTseq2VCF",
    version="0.1.0",
    author="Bradley T. Martin, Ph.D.",
    author_email="evobio721@gmail.com",
    description="A tool to convert GTseq data to VCF format and merge with ddRADseq VCF file.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/btmartin721/GTseq2VCF",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "matplotlib",
        "seaborn",
        "pysam",
    ],
    python_requires=">=3.11",
    entry_points={
        "console_scripts": [
            "gtseq2vcf=gtseq2vcf.py:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="genomics GTseq VCF ddRADseq bioinformatics merge",
)
