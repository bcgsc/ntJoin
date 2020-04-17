#!/usr/bin/env python3

"""
Scans through KMC output of 2 files, prints out kmers in common between the two
"""

import argparse

def main():
    "Only print kmers in common between the KMC input files"
    parser = argparse.ArgumentParser(description="Find common kmers between KMC output")
    parser.add_argument("-t", help="Target KMC file", required=True)
    parser.add_argument("-r", help="Reference KMC file", required=True)
    args = parser.parse_args()

    ref_kmers = set()

    with open(args.r, 'r') as ref_kmers_file:
        for kmer in ref_kmers_file:
            kmer, count = kmer.strip().split("\t")
            ref_kmers.add(kmer)

    with open(args.t, 'r') as target_kmers_file:
        for kmer in target_kmers_file:
            kmer, count = kmer.strip().split("\t")
            if kmer in ref_kmers:
                print(kmer, count, sep="\t")


if __name__ == "__main__":
    main()
