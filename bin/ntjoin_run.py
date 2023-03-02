#!/usr/bin/env python3
"""
ntJoin: Scaffolding assemblies and analyzing synteny
using reference assemblies and minimizer graphs
Written by Lauren Coombe (@lcoombe)
"""
import argparse
import sys
from ntjoin_synteny import NtjoinSynteny
from ntjoin_assemble import NtjoinScaffolder

def parse_arguments():
    "Parse ntJoin arguments"
    parser = argparse.ArgumentParser(
        description="ntJoin: Genome analysis using reference assemblies and minimizer graphs",
        )
    parser.add_argument("-v", "--version", action='version', version='ntJoin v1.1.1')

    subparsers = parser.add_subparsers(dest="mode")
    scaffold_parser = subparsers.add_parser("scaffold",
    help="Scaffold the input target assembly using the supplied reference(s)",
                epilog="Note: Script expects that each input minimizer TSV file has a matching fasta file.\n"
                "Example: myscaffolds.fa.k32.w1000.tsv - myscaffolds.fa is the expected matching fasta",
        formatter_class=argparse.RawTextHelpFormatter)
    scaffold_parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of references")
    scaffold_parser.add_argument("-s", help="Target scaffolds minimizer TSV file", required=True)
    scaffold_parser.add_argument("-l", help="Weight of target genome assembly [1]",
                        required=False, default=1, type=float)
    scaffold_parser.add_argument("-r",
                        help="List of reference assembly weights (in quotes, separated by spaces, "
                                "in same order as minimizer TSV files)",
                        required=True, type=str)
    scaffold_parser.add_argument("-p", help="Output prefix [out]", default="out",
                        type=str, required=False)
    scaffold_parser.add_argument("-n", help="Minimum edge weight [1]", default=1, type=int)
    scaffold_parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
    scaffold_parser.add_argument("-g", help="Minimum gap size (bp)", required=False, default=20, type=int)
    scaffold_parser.add_argument("-G", help="Maximum gap size (bp) (0 if no maximum threshold)", required=False,
                        default=0, type=int)
    scaffold_parser.add_argument("--mkt", help="Use Mann-Kendall Test to orient contigs (slower, overrides m)",
                        action='store_true')
    scaffold_parser.add_argument('-m', help="Require at least m %% of minimizer positions to be "
                                    "increasing/decreasing to assign contig orientation [90]\n "
                                    "Note: Only used with --mkt is NOT specified", default=90, type=int)
    scaffold_parser.add_argument('-t', help="Number of threads for multiprocessing [1]", default=1, type=int)
    scaffold_parser.add_argument("-v", "--version", action='version', version='ntJoin v1.1.1')
    scaffold_parser.add_argument("--agp", help="Output AGP file describing scaffolds", action="store_true")
    scaffold_parser.add_argument("--no_cut", help="Do not cut input contigs, place in most representative path",
                        action="store_true")
    scaffold_parser.add_argument("--overlap", help="Attempt to detect and trim overlapping joined sequences",
                        action="store_true")
    scaffold_parser.add_argument("--overlap_gap",
                        help="Length of gap introduced between overlapping, trimmed segments [20]",
                        type=int, default=20)
    scaffold_parser.add_argument("--overlap_k", help="Kmer size used for overlap minimizer step [15]",
                        type=int, default=15)
    scaffold_parser.add_argument("--overlap_w", help="Window size used for overlap minimizer step [10]",
                        type=int, default=10)
    scaffold_parser.add_argument("--btllib_t", help="Number of threads for btllib wrapper functions "
                                            "(computing minimizers, reading fasta file) [4]",
                        type=int, default=4)

    synteny_parser = subparsers.add_parser("synteny", help="Extract syntenic blocks from input assemblies")
    synteny_parser.add_argument("FILES", nargs="+", help="Minimizer TSV files of input assemblies")
    synteny_parser.add_argument("-n", help="Minimum edge weight [Number of input assemblies]", default=0, type=int)
    synteny_parser.add_argument("-p", help="Output prefix [out]",
                                default="out", type=str, required=False)
    synteny_parser.add_argument("-k", help="Kmer size used for minimizer step", required=True, type=int)
    synteny_parser.add_argument("-w", help="Window size used for minimizers", required=True, type=int)
    synteny_parser.add_argument("--btllib_t", help="Number of threads for btllib wrapper functions "\
                                "(computing minimizers, reading fasta file) [4]", type=int, default=4)
    synteny_parser.add_argument("--w-rounds", help="decreasing list of 'w' values to use for refining ends",
                                default=[100, 10, 5], nargs="+", type=int)
    synteny_parser.add_argument("--dev", action="store_true", help="Developer mode - retain intermediate files")
    synteny_parser.add_argument("-v", "--version", action='version', version='ntJoin v1.1.1')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()

def main():
    "Run ntJoin"
    args = parse_arguments()
    if args.mode == "scaffold":
        NtjoinScaffolder(args).main_scaffolder()
    elif args.mode == "synteny":
        NtjoinSynteny(args).main_synteny()
    else:
        raise ValueError(f"Unexpected mode: {args.mode}")

if __name__ == "__main__":
    main()
