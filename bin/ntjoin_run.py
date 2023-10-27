#!/usr/bin/env python3
"""
ntJoin: Scaffolding assemblies and analyzing synteny
using reference assemblies and minimizer graphs
Written by Lauren Coombe (@lcoombe)
"""
import argparse
import sys
from ntjoin_assemble import NtjoinScaffolder

def parse_arguments():
    "Parse ntJoin arguments"
    parser = argparse.ArgumentParser(
            description="ntJoin: Scaffolding genome assemblies using reference assemblies and minimizer graphs",
            epilog="Note: Script expects that each input minimizer TSV file has a matching fasta file.\n"
                   "Example: myscaffolds.fa.k32.w1000.tsv - myscaffolds.fa is the expected matching fasta",
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of references")
    parser.add_argument("-s", help="Target scaffolds minimizer TSV file", required=True)
    parser.add_argument("-l", help="Weight of target genome assembly [1]",
                        required=False, default=1, type=float)
    parser.add_argument("-r",
                        help="List of reference assembly weights (in quotes, separated by spaces, "
                                "in same order as minimizer TSV files)",
                        required=True, type=str)
    parser.add_argument("-p", help="Output prefix [out]", default="out",
                        type=str, required=False)
    parser.add_argument("-n", help="Minimum edge weight [1]", default=1, type=int)
    parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
    parser.add_argument("-g", help="Minimum gap size (bp)", required=False, default=20, type=int)
    parser.add_argument("-G", help="Maximum gap size (bp) (0 if no maximum threshold)", required=False,
                        default=0, type=int)
    parser.add_argument("--mkt", help="Use Mann-Kendall Test to orient contigs (slower, overrides m)",
                        action='store_true')
    parser.add_argument('-m', help="Require at least m %% of minimizer positions to be "
                                    "increasing/decreasing to assign contig orientation [90]\n "
                                    "Note: Only used with --mkt is NOT specified", default=90, type=int)
    parser.add_argument('-t', help="Number of threads for multiprocessing [1]", default=1, type=int)
    parser.add_argument("-v", "--version", action='version', version='ntJoin v1.1.1')
    parser.add_argument("--agp", help="Output AGP file describing scaffolds", action="store_true")
    parser.add_argument("--no_cut", help="Do not cut input contigs, place in most representative path",
                        action="store_true")
    parser.add_argument("--overlap", help="Attempt to detect and trim overlapping joined sequences",
                        action="store_true")
    parser.add_argument("--overlap_gap",
                        help="Length of gap introduced between overlapping, trimmed segments [20]",
                        type=int, default=20)
    parser.add_argument("--overlap_k", help="Kmer size used for overlap minimizer step [15]",
                        type=int, default=15)
    parser.add_argument("--overlap_w", help="Window size used for overlap minimizer step [10]",
                        type=int, default=10)
    parser.add_argument("--btllib_t", help="Number of threads for btllib wrapper functions "
                                            "(computing minimizers, reading fasta file) [4]",
                        type=int, default=4)


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()

def main():
    "Run ntJoin"
    args = parse_arguments()
    NtjoinScaffolder(args).main_scaffolder()


if __name__ == "__main__":
    main()
