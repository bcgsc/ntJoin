## minimizer-assembly-scaffolder

Scaffolding multiple assemblies using minimizer graphs

#### Usage: 
```
Usage: minimizer_assembler-make assemble prefix=<prefix> list_files='List of fasta files' list_weights='List of weights (ints) per assembly'

Note: ensure the lists of assemblies and weights are in the same order, and that both are space-separated
Commands:
t	Number of threads [4]
k	K-mer size for minimizers [64]
w	Window size for minimizers [1000]
n	Min. edge weight [2]
```

#### Requirements:
* python3 ([pybedtools](https://daler.github.io/pybedtools/), [networkx](https://networkx.github.io/))
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [samtools](https://github.com/samtools/samtools)