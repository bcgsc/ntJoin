[![Release](https://img.shields.io/github/release/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/releases)
[![Issues](https://img.shields.io/github/issues/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/issues)
[![Conda](https://img.shields.io/conda/dn/bioconda/ntjoin?label=Conda)](https://anaconda.org/bioconda/ntjoin)

![Logo](https://github.com/bcgsc/ntJoin/blob/master/ntjoin-logo.png)

# ntJoin

Scaffolding draft assemblies using reference assemblies and minimizer graphs

## Description of the algorithm

ntJoin takes a target assembly and one or more 'reference' assembly as input, and uses information from the reference(s) to scaffold the target assembly. The 'reference' assemblies can be true reference assembly builds, or a different draft genome assemblies.

Instead of using costly alignments, ntJoin uses a more lightweight approach using minimizer graphs to yield a mapping between the input assemblies. 

**Main steps in the algorithm:**

1. Generate an ordered minimizer sketch for each contig of each input assembly
2. Filter the minimizers to only retain minimizers that are:
    * Unique within each assembly
    * Found in all assemblies (target + all references)
3. Build a minimizer graph
    * Nodes: minimizers
    * Edges: between minimizers that are adjacent in at least one of the assemblies. Edge weights are the sum of weights of the assemblies that support an edge.
4. Filter the graph based on the minimum edge weight (`n`)
5. For each node that is a branch node (degree > 2), filter the incident edges with an increasing edge threshold
6. Each linear path is converted to a list of oriented target assembly contig regions to scaffold together
7. Target assembly scaffolds are printed out


### Credits

Original concept: Rene Warren

Design and implementation: Lauren Coombe


### Citing ntJoin

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/stargazers) and for using, developing and promoting this free software!

If you use ntJoin in your research, please cite:

Lauren Coombe, Vladimir Nikolic, Justin Chu, Inanc Birol, Rene L Warren: ntJoin: Fast and lightweight assembly-guided scaffolding using minimizer graphs. _Bioinformatics_ (2020) doi: https://doi.org/10.1093/bioinformatics/btaa253.


## Usage

```
Usage: ntJoin assemble target=<target scaffolds> references='List of reference assemblies' reference_weights='List of weights per reference assembly'

Options:
target			Target assembly to be scaffolded in fasta format
references		List of reference files (separated by a space, in fasta format)
target_weight		Weight of target assembly [1]
reference_weights	List of weights of reference assemblies
prefix			Prefix of intermediate output files [out.k<k>.w<w>.n<n>]
t			Number of threads [4]
k			K-mer size for minimizers [32]
w			Window size for minimizers (bp) [1000]
n			Minimum graph edge weight [1]
g			Minimum gap size (bp) [20]
G           		Maximum gap size (bp) (0 if no maximum) [0]
m			Minimum percentage of increasing/decreasing minimizer positions to orient contig [90]
mkt			If True, use Mann-Kendall Test to predict contig orientation (computationally-intensive, overrides 'm') [False]
agp			If True, output AGP file describing output scaffolds [False]
no_cut			If True, will not cut contigs at putative misassemblies [False]

Notes: 
	- Ensure the lists of reference assemblies and weights are in the same order, and that both are space-separated
	- Ensure all assembly files are in the current working directory
```

Running `ntJoin help` prints the help documentation.

### Example

* Target assembly to scaffold: my_scaffolds.fa 
* Assembly to use as 'reference': assembly_ref1.fa
* Giving the target asssembly a weight of '1' and reference assembly a weight of '2'
* Using k=32, w=500
* **Ensure that all input assembly files are in or have soft-links to the current working directory**

```
ntJoin assemble target=my_scaffolds.fa target_weight=1 references='assembly_ref1.fa' reference_weights='2' k=32 w=500
```

### Output files

* Scaffolded targeted assembly (`<target assembly>.k<k>.w<w>.n<n>.all.scaffolds.fa`)
* Path file describing how target assembly was scaffolded (`<prefix>.path`)
* Unfiltered minimizer graph in dot format (`<prefix>.mx.dot`)
* If agp=True specified, AGP describing how target assembly was scaffolded (`<prefix>.agp`)

### Parameter considerations

* We recommend setting the reference weight(s) to be higher than the target weight

## Installation Instructions
#### Installing ntJoin using Brew
ntJoin can be installed using [Homebrew](https://brew.sh) on macOS or [Linuxbrew](http://linuxbrew.sh) on Linux:
```sh
brew install brewsci/bio/ntjoin
```

#### Installing ntJoin using Conda
```sh
conda install -c bioconda ntjoin
```

#### Installing ntJoin from the source code
```sh
git clone https://github.com/bcgsc/ntJoin.git
cd ntJoin/src
make
```

## Dependencies

* python3 ([pybedtools](https://daler.github.io/pybedtools/), [python-igraph](https://igraph.org/python/), [pymannkendall](https://pypi.org/project/pymannkendall/))
* [bedtools v2.29.2+](https://bedtools.readthedocs.io/en/latest/)
* [samtools](https://github.com/samtools/samtools)
* [zlib](https://www.zlib.net/)

Python dependencies can be installed with:
```sh
pip3 install -r requirements.txt
```

## Testing your Installation
See `tests/test_installation.sh` to test your ntJoin installation and see an example command.

## License

ntJoin Copyright (c) 2020 British Columbia Cancer Agency Branch.  All rights reserved.

ntJoin is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>

