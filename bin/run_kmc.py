#!/usr/bin/env python3
"""
Runs KMC kmer-counting on an assembly, basing the maximum memory on the supplied genome size (bp)
Written by Lauren Coombe (@lcoombe)
"""

import argparse
import re
import shlex
import subprocess
import sys
from read_fasta import read_fasta

def calculate_max_mem(genome_size, max_mem_arg):
    "Calculate the max memory for KMC based on genome size or specified max memory"
    assert genome_size is not None or max_mem_arg is not None
    if genome_size is not None:
        max_mem = genome_size/1e9 * 2.0
        max_mem = max(max_mem, 1) # KMC require at least 1GB of RAM
        return max_mem
    mem_gb = 0
    gb_match = re.search(r'(\d+)G', max_mem_arg, flags=re.IGNORECASE)
    mb_match = re.search(r'(\d+)M', max_mem_arg, flags=re.IGNORECASE)
    kb_match = re.search(r'(\d+)k', max_mem_arg, flags=re.IGNORECASE)
    bytes_match = re.search(r'^(\d+)$', max_mem_arg)
    if gb_match:
        mem_gb = int(gb_match.group(1))
    elif mb_match:
        mem_gb = int(mb_match.group(1))/1024
    elif kb_match:
        mem_gb = int(kb_match.group(1))/(1024*1024)
    elif bytes_match:
        mem_gb = int(bytes_match.group(1))/(1024*1024*1024)
    else:
        print("ERROR: could not parse max mem format:", max_mem_arg)
        sys.exit(1)
    max_mem = max(mem_gb, 1)
    max_mem = min(max_mem, 1024)
    return max_mem


def main():
    "Run KMC kmer counting"
    parser = argparse.ArgumentParser(description="Run KMC kmer counting on an assembly",
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog="Either -g or -m must be set")
    parser.add_argument("FASTA", help="Fasta file")
    parser.add_argument("-l", help="Minimum length threshold (bp) [500]", default=500, type=int)
    parser.add_argument("-g", help="Approximate genome size (bp)", type=float)
    parser.add_argument("-m", help="Max memory (SI units)", type=str)
    parser.add_argument("-k", help="Kmer size", required=True, type=int)
    parser.add_argument("-t", help="Number of threads [4]", type=int, default=4)
    parser.add_argument("-c", help="Lower threshold for kmer counting [2]", default=2, type=int)
    parser.add_argument("-p", help="Output prefix [fastafile.k<k>.w<w>.kmc]", type=str, required=False)
    args = parser.parse_args()

    if args.g is None and args.m is None:
        print("ERROR: Either -g or -m must be set")
        sys.exit(1)

    # Make TMP directory
    tmpdir = args.FASTA + "k" + str(args.k) + ".w" + str(args.l) + ".kmc.tmp"
    cmd = "mkdir -p " + tmpdir
    print(cmd)
    cmd_shlex = shlex.split(cmd)
    ret = subprocess.call(cmd_shlex)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd_shlex)

    # Filter the fasta file based on length threshold
    filtered_fasta_name = args.FASTA + "k" + str(args.k) + ".w" + str(args.l)+ "." + str(args.l) + "plus.fa"
    filtered_fasta = open(filtered_fasta_name, 'w')
    with open(args.FASTA, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            if len(seq) >= args.l:
                filtered_fasta.write(">" + header + "\n" + seq + "\n")
    filtered_fasta.close()

    # Run KMC (1st step)
    max_mem = calculate_max_mem(args.g, args.m)
    kmc_out_prefix = args.p if args.p is not None else "%s.k%d.w%d.kmc" % (args.FASTA, args.k, args.l)
    cmd = "kmc -ci%d -k%d -m%d -t%d -fm %s %s %s" % \
          (args.c, args.k, max_mem, args.t, filtered_fasta_name, kmc_out_prefix, tmpdir)
    print(cmd)
    cmd_shlex = shlex.split(cmd)
    ret = subprocess.call(cmd_shlex)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd_shlex)

    # Run KMC (2nd step)
    cmd = "kmc_dump %s %s" % (kmc_out_prefix, kmc_out_prefix + ".tsv")
    print(cmd)
    cmd_shlex = shlex.split(cmd)
    ret = subprocess.call(cmd_shlex)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd_shlex)

    # Clean-up tmp files
    cmd = "rm -r %s %s" % (tmpdir, filtered_fasta_name)
    cmd_shlex = shlex.split(cmd)
    ret = subprocess.call(cmd_shlex)
    if ret != 0:
        raise subprocess.CalledProcessError(ret, cmd_shlex)


if __name__ == "__main__":
    main()
