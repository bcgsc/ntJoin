[![Release](https://img.shields.io/github/release/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/ntJoin/total?logo=github)](https://github.com/bcgsc/ntJoin/releases/download/v1.0.1/ntJoin-1.0.1.tar.gz)
[![Issues](https://img.shields.io/github/issues/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/issues)

![Logo](https://github.com/bcgsc/ntJoin/blob/master/ntjoin-logo.png)

# ntJoin

Scaffolding assemblies using reference assemblies and minimizer graphs

## Description of the algorithm
------------

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

### Citing ntJoin
------------

Thank you for your [![Stars](https://img.shields.io/github/stars/bcgsc/ntJoin.svg)](https://github.com/bcgsc/ntJoin/stargazers) and for using, developing and promoting this free software!

If you use ntJoin in your research, please cite:
TBD


## Usage
------------

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
m			Minimum percentage of increasing/decreasing minimizer positions to orient contig [90]
mkt			If True, use Mann-Kendall Test to predict contig orientation (computationally-intensive, overrides 'm') [False]

Note: ensure the lists of reference assemblies and weights are in the same order, and that both are space-separated
```

### Example
------------

* Target assembly to scaffold: my_scaffolds.fa 
* Two assemblies to use as 'references': assembly_ref1.fa, assembly_ref2.fa
* Giving the target asssembly a weight of '1' and each reference assembly a weight of '2'
* Using k=32, w=500

```
ntJoin assemble target=my_scaffolds.fa target_weight=1 references='assembly_ref1.fa assembly_ref2.fa' reference_weights='2 2' k=32 w=500
```

### Output files
------------

* Scaffolded targeted assembly (`<target assembly>.k<k>.w<w>.n<n>.all.scaffolds.fa`)
* Path file describing how target assembly was scaffolded (`<prefix>.mx.path`)
* Unfiltered minimizer graph in dot format (`<prefix>.mx.dot`)

## Installation Instructions
------------

#### Installing ntJoin from the source code
```sh
git clone https://github.com/bcgsc/ntJoin.git
cd src
make
```

## Dependencies
------------

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

## License
------------

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

