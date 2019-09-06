# ntJoin

Scaffold multiple assemblies using minimizer graphs

## Usage
```
Usage: ntJoin assemble prefix=<prefix> list_files='List of fasta files' list_weights='List of weights (ints) per assembly'

Note: ensure the lists of assemblies and weights are in the same order, and the files in the lists are separated by a space

Options:
prefix	Prefix of intermediate output files [out]
t	Number of threads [4]
k	K-mer size for minimizers [32]
w	Window size for minimizers [1000]
n	Minimum edge weight [2]
g	Minimum gap size [20]
m	Minimum percentage of increasing/decreasing minimizer positions to orient contig [90]"
mkt	If True, use Mann-Kendall Test to predict contig orientation (more computationally-intensive) [False]"
```
--------
## Installation Instructions

#### Installing ntJoin from the source code
```sh
git clone https://github.com/bcgsc/ntJoin.git
cd src
make
```

## Dependencies
* python3 ([pybedtools](https://daler.github.io/pybedtools/), [python-igraph](https://igraph.org/python/), [pymannkendall](https://pypi.org/project/pymannkendall/))
* [bedtools v2.21.0+](https://bedtools.readthedocs.io/en/latest/)
* [samtools](https://github.com/samtools/samtools)

Python dependencies can be installed with:
```sh
pip3 install -r requirements.txt
```
Bedtools and samtools can be installed using [Homebrew](https://brew.sh) on macOS or [Linuxbrew](http://linuxbrew.sh) on Linux:
```sh
brew install samtools bedtools
```
