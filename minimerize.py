#!/usr/bin/env pypy3
"""
Compute the minimizers of a nucleotide sequence.
Written by Shaun Jackman @sjackman
Edits by Lauren Coombe @lcoombe
"""

import argparse
import re
from hash import hash_kmer
from read_fasta import read_fasta

ACGT = re.compile("^[ACGT]+$")

def kmerize(k, seq):
    "Iterator over the kmers of a string."
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        if ACGT.match(kmer):
            yield (kmer, i)

def minimerize(k, w, seq):
    "Return the minimizers of a string."
    hashes = [hash_kmer(kmer) for kmer in kmerize(k, seq)]
    minimizers = []
    previous_minimizer = -1
    for i in range(0, len(hashes) - w + 1):
        minimizer, minimizer_i, minimizer_pos = min((x[0], j, x[1]) for (j, x) in enumerate(hashes[i : i + w]))
        minimizer_i += i
        if minimizer_i > previous_minimizer:
            previous_minimizer = minimizer_i
            minimizers.append((minimizer, minimizer_pos))
    return minimizers

def main():
    parser = argparse.ArgumentParser(description="Minimizerize sequences, with positions tracked")
    parser.add_argument("-k", help="k-mer size [32]", default=32, type=int)
    parser.add_argument("-w", help="window size [32]", default=32, type=int)
    parser.add_argument("fasta", help="Fasta file")
    args = parser.parse_args()

    with open(args.fasta) as fin:
        for name, seq, bx, _ in read_fasta(fin):
            print(name, sep="", end="\t")
            mxs_formatted = ["{0[0]}:{0[1]}".format(tup) for tup in minimerize(args.k, args.w, seq.upper())]
            print(*mxs_formatted)



if __name__ == "__main__":
    main()