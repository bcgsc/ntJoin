#!/usr/bin/env python3
"""
ntJoin: Identifying synteny between genome assemblies using minimizer graphs
Written by Lauren Coombe @lcoombe
"""

from collections import namedtuple
import re
import sys
import ntjoin_utils
import pybedtools

# Named tuples
Minimizer = namedtuple("Minimizer", ["mx", "position"])

# Regexes
fai_re = re.compile(r'^(\S+).k\d+.w\d+.tsv')

class SyntenyBlock:
    "A Synteny Block between the input assemblies"
    def __init__(self, k, *assemblies):
        "Instantiate a dictionary to keep track of assembly blocks for this synteny block"
        self.assembly_blocks = {assembly: AssemblyBlock() for assembly in assemblies}
        self.k = k # k-mer size used for minimizers, needed to adjust end coordinates


    def continue_block(self, mx, list_mx_info):
        "Given minimizer and preliminary blocks, return if synteny block should extend, else False"
        return all(mx_dict[mx][0] == self.assembly_blocks[assembly].contig_id \
            for assembly, mx_dict in list_mx_info.items())

    def extend_block(self, mx, list_mx_jnfo):
        "Extend the synteny block by extending each assembly block"
        for assembly, mx_dict in list_mx_jnfo.items():
            ctg, pos = mx_dict[mx]
            assert self.assembly_blocks[assembly].contig_id == ctg
            self.assembly_blocks[assembly].minimizers.append(Minimizer(mx, int(pos)))

    def start_block(self, mx, list_mx_info):
        "Start the new synteny block"
        for assembly, mx_dict in list_mx_info.items():
            ctg, pos = mx_dict[mx]
            self.assembly_blocks[assembly].contig_id = ctg
            self.assembly_blocks[assembly].minimizers.append(Minimizer(mx, int(pos)))

    def determine_orientations(self):
        "Determine the orientations of each assembly block"
        for _, assembly_block in self.assembly_blocks.items():
            positions = [mx.position for mx in assembly_block.minimizers]
            if all(x < y for x, y in zip(positions, positions[1:])):
                assembly_block.ori = "+"
            elif all(x > y for x, y in zip(positions, positions[1:])):
                assembly_block.ori = "-"
            else:
                assembly_block.ori = "?"

    def all_oriented(self):
        "Return true if all of the assembly blocks in the synteny block are oriented"
        return all(assembly_block.ori in ["+", "-"] for _, assembly_block in self.assembly_blocks.items())

    def get_block_string(self, num):
        "Given the specified synteny block ID, print the synteny blocks"
        return_str = ""
        for assembly, assembly_block in self.assembly_blocks.items():
            start_pos = assembly_block.get_block_start()
            end_pos = assembly_block.get_block_end() + self.k
            block_string = f"{num}\t{assembly}\t{assembly_block.contig_id}\t{start_pos}" \
                f"\t{end_pos}\t{assembly_block.ori}\n"
            return_str += block_string
        return return_str


class AssemblyBlock:
    "An assembly block for a given assembly. The AssemblyBlock objects per assembly make up a SyntenyBlock"
    def __init__(self):
        "Instantiate the AssemblyBlock"
        self.contig_id = None
        self.minimizers = []
        self.ori = None

    def get_block_start(self):
        "Get the starting coordinate of the assembly block"
        return min(self.minimizers[0].position, self.minimizers[-1].position)

    def get_block_end(self):
        "Get the end coordinate of the assembly block"
        return max(self.minimizers[0].position, self.minimizers[-1].position)


def find_synteny_blocks(path, list_mx_info, k):
    "Given a path (sequence of mx), print the order/orientation/regions of contigs for an assembly"
    out_blocks = []  # List of SyntenyBlock
    prelim_blocks = SyntenyBlock(k, *list(list_mx_info.keys()))
    past_start_flag = False
    for mx in path:
        if prelim_blocks.continue_block(mx, list_mx_info):
            prelim_blocks.extend_block(mx, list_mx_info)
        else:
            # This is either the first mx, or we are past a stretch of repeating contigs
            if past_start_flag:
                prelim_blocks.determine_orientations()
                if prelim_blocks.all_oriented():
                    out_blocks.append(prelim_blocks)
            prelim_blocks = SyntenyBlock(k, *list(list_mx_info.keys()))
            prelim_blocks.start_block(mx, list_mx_info)

    prelim_blocks.determine_orientations()
    if prelim_blocks.all_oriented():
        out_blocks.append(prelim_blocks)

    return out_blocks

def find_fa_name(assembly_mx_name):
    "Given the mx file name, return the corresponding fai file name"
    if fai_match := re.search(fai_re, assembly_mx_name):
        return f"{fai_match.group(1)}"
    print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
    print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
    sys.exit(1)

def get_synteny_bed_lists(paths, w):
    "Given a set of synteny blocks, return a dictionary with a Bed interval lists per contig, per assembly"
    synteny_beds = {}
    for subcomponent in paths:
        for block in subcomponent:
            for assembly, assembly_block in block.assembly_blocks.items():
                if assembly not in synteny_beds:
                    synteny_beds[assembly] = {}
                if assembly_block.contig_id not in synteny_beds[assembly]:
                    synteny_beds[assembly][assembly_block.contig_id] = []
                synteny_beds[assembly][assembly_block.contig_id].append(
                    ntjoin_utils.Bed(assembly_block.contig_id,
                                    assembly_block.get_block_start() + w,
                                    assembly_block.get_block_end() + block.k - w))

    return synteny_beds

def mask_assemblies_with_synteny_extents(synteny_beds):
    "Mask each reference assembly with determined synteny blocks"
    mx_to_fa_dict = {}
    for assembly, contig_dict in synteny_beds.items():
        bed_str = [f"{ctg}\t{bed.start}\t{bed.end}\tSYNTENY" for ctg in contig_dict \
                    for bed in contig_dict[ctg]]
        bed_str = "\n".join(bed_str)
        synteny_bed = pybedtools.BedTool(bed_str, from_string=True).sort()
        fa_filename = find_fa_name(assembly)
        synteny_bed.mask_fasta(fi=fa_filename, fo=f"{fa_filename}_masked.fa")
        mx_to_fa_dict[assembly] = f"{fa_filename}_masked.fa"
    return mx_to_fa_dict

def generate_new_minimizers(tsv_to_fa_dict, k, w, t):
    "Given the masked fasta files, generate minimizers at new w for each"
    list_mx_info = {}
    list_mxs = {}
    for assembly_tsv, assembly_masked in tsv_to_fa_dict.items():
        indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, k, int(w/10), t)
        mx_info, mxs_filt = ntjoin_utils.read_minimizers(indexlr_filename)
        list_mx_info[assembly_tsv] = mx_info
        list_mxs[assembly_tsv] = mxs_filt
    return list_mx_info, list_mxs


def generate_additional_minimizers(paths, w, t):
    "Given the existing synteny blocks, generate minimizers for increased block resolution"
    k = paths[0][0].k
    synteny_beds = get_synteny_bed_lists(paths, w)
    mx_to_fa_dict = mask_assemblies_with_synteny_extents(synteny_beds)
    list_mx_info, list_mxs = generate_new_minimizers(mx_to_fa_dict, k, w, t)
