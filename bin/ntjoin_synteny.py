#!/usr/bin/env python3
"""
ntJoin: Identifying synteny between genome assemblies using minimizer graphs
Written by Lauren Coombe @lcoombe
"""

from collections import namedtuple, defaultdict
import re
import shlex
import subprocess
import sys
import intervaltree
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

    def get_block_terminal_mx(self):
        "Return the terminal minimizer hashes for the assembly block"
        return self.contig_id, self.minimizers[0], self.minimizers[-1]

    def get_block_internal_mx_hashes(self):
        "Return the internal minimizer hashes for the assembly block"
        return [mx_pos.mx for mx_pos in self.minimizers[1:-1]]


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

def get_synteny_bed_lists(paths):
    "Given a set of synteny blocks, return a dictionary with a BED interval list per contig, per assembly"
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
                                    assembly_block.get_block_start(),
                                    assembly_block.get_block_end() + block.k))

    return synteny_beds

def mask_assemblies_with_synteny_extents(synteny_beds, w):
    "Mask each reference assembly with determined synteny blocks"
    mx_to_fa_dict = {}
    for assembly, contig_dict in synteny_beds.items():
        bed_str = [f"{ctg}\t{bed.start}\t{bed.end}\tSYNTENY" for ctg in contig_dict \
                    for bed in contig_dict[ctg] if bed.end - bed.start > 2*w]
        bed_str = "\n".join(bed_str)
        fa_filename = find_fa_name(assembly)
        synteny_bed = pybedtools.BedTool(bed_str, from_string=True).slop(g=f"{fa_filename}.fai", l=-1*w, r=-1*w).sort()
        synteny_bed.mask_fasta(fi=fa_filename, fo=f"{fa_filename}_masked.fa")
        mx_to_fa_dict[assembly] = f"{fa_filename}_masked.fa"
    return mx_to_fa_dict

def delete_w_iteration_files(*filenames):
    "Delete the given files for the specific lower w iteration"
    for filename in filenames:
        cmd = shlex.split(f"rm {filename}")
        ret_code = subprocess.call(cmd)
        assert ret_code == 0

def generate_new_minimizers(tsv_to_fa_dict, k, w, t, retain_files=False):
    "Given the masked fasta files, generate minimizers at new w for each"
    list_mxs = {}
    new_list_mxs_info = {}
    for assembly_tsv, assembly_masked in tsv_to_fa_dict.items():
        indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, k, w, t)
        mx_info, mxs_filt = ntjoin_utils.read_minimizers(indexlr_filename)
        new_list_mxs_info[assembly_tsv] = mx_info
        list_mxs[assembly_tsv] = mxs_filt
        if not retain_files:
            delete_w_iteration_files(indexlr_filename, assembly_masked)
    return list_mxs, new_list_mxs_info

def update_interval_tree(trees, assembly_name, ctg, mx1, mx2):
    "Update the given dictionary of trees with the new extent"
    start_pos = min(mx1.position, mx2.position)
    end_pos = max(mx1.position, mx2.position)
    if assembly_name not in trees or ctg not in trees[assembly_name]:
        trees[assembly_name][ctg] = intervaltree.IntervalTree()
    if trees[assembly_name][ctg][start_pos+1:end_pos]: # Checking that this doesn't overlap with anything
        print("WARNING: detected overlapping segments:", assembly_name, ctg, start_pos+1, end_pos,
                file=sys.stderr)
    trees[assembly_name][ctg][start_pos+1:end_pos] = (mx1, mx2)

def find_mx_in_blocks(paths):
    "Given the synteny blocks, find the minimizers at the terminal ends of each block, and internal"
    terminal_mxs = set()
    internal_mxs = set()
    intervaltrees = defaultdict(dict) # assembly -> contig -> IntervalTree of synteny block extents

    for subcomponent in paths:
        for block in subcomponent:
            curr_mx_len = len(terminal_mxs)
            for assembly, assembly_block in block.assembly_blocks.items():
                contig, mx1, mx2 = assembly_block.get_block_terminal_mx()
                terminal_mxs.add(mx1.mx)
                terminal_mxs.add(mx2.mx)
                update_interval_tree(intervaltrees, assembly, contig, mx1, mx2)
                internal = assembly_block.get_block_internal_mx_hashes()
                internal_mxs = internal_mxs.union(internal)
            assert len(terminal_mxs) == (curr_mx_len + 2)
    return terminal_mxs, internal_mxs, intervaltrees

def check_non_overlapping(paths):
    "Given the paths, do final check to ensure intervals are not overlapping, will print warnings if that's the case"
    intervaltrees = defaultdict(dict) # assembly -> contig -> IntervalTree of synteny block extents
    for subcomponent in paths:
        for block in subcomponent:
            for assembly, assembly_block in block.assembly_blocks.items():
                contig, mx1, mx2 = assembly_block.get_block_terminal_mx()
                update_interval_tree(intervaltrees, assembly, contig, mx1, mx2)


def filter_minimizers_synteny_blocks(list_mxs, black_list, intervaltrees, list_mx_info):
    "Filter minimizers found in the mx black list"
    return_mxs = {}
    for assembly in list_mxs:
        assembly_mxs_filtered = []
        for mx_list in list_mxs[assembly]:
            new_list = []
            for mx in mx_list:
                ctg, pos = list_mx_info[assembly][mx]
                i_tree = intervaltrees[assembly][ctg]
                if new_list:
                    prev_pos = list_mx_info[assembly][new_list[-1]][1]
                    start = min(prev_pos, pos)
                    end = max(prev_pos, pos)
                    if i_tree[start:end]: # Split the mx adjacency if it spans over a known synteny block
                        assembly_mxs_filtered.append(new_list)
                        new_list = []
                if mx not in black_list and not i_tree[pos]:
                    new_list.append(mx)
            assembly_mxs_filtered.append(new_list)

        return_mxs[assembly] = assembly_mxs_filtered
    return return_mxs

def update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info):
    "Update the directory containing mx -> contig, position associations"
    valid_mxs = set({mx for _, list_mx_val in list_mxs.items() \
                    for list_mx in list_mx_val for mx in list_mx})
    for assembly, mx_dict in new_list_mx_info.items():
        for mx in mx_dict:
            if mx in valid_mxs and mx not in list_mx_info[assembly]:
                list_mx_info[assembly][mx] = mx_dict[mx]


def generate_additional_minimizers(paths, new_w, prev_w, t, list_mx_info, dev=False):
    "Given the existing synteny blocks, generate minimizers for increased block resolution"
    k = paths[0][0].k
    synteny_beds = get_synteny_bed_lists(paths)
    mx_to_fa_dict = mask_assemblies_with_synteny_extents(synteny_beds, prev_w)
    list_mxs, new_list_mx_info = generate_new_minimizers(mx_to_fa_dict, k, new_w, t, retain_files=dev)
    terminal_mx, internal_mx, interval_trees = find_mx_in_blocks(paths)
    list_mxs = filter_minimizers_synteny_blocks(list_mxs, internal_mx, interval_trees, new_list_mx_info)
    list_mxs = ntjoin_utils.filter_minimizers(list_mxs) # Filter for mx in all assemblies
    update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info)
    return list_mxs, terminal_mx
