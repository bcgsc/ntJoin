#!/usr/bin/env python3
"""
ntJoin: Identifying synteny between genome assemblies using minimizer graphs
Written by Lauren Coombe @lcoombe
"""

from collections import defaultdict
import datetime
import re
import shlex
import subprocess
import sys
import intervaltree
import ntjoin_utils
import ntjoin
import pybedtools
from synteny_block import SyntenyBlock


# Regexes
fai_re = re.compile(r'^(\S+).k\d+.w\d+.tsv')

class NtjoinSynteny(ntjoin.Ntjoin):
    "Instance for ntJoin synteny mode"

    def __init__(self, args):
        super().__init__(args)
        self.weights_list = [1] * len(self.args.FILES)
        self.args.t = 1
        self.print_parameters_synteny()

    def print_parameters_synteny(self):
        "Pring the set parameters for the ntJoin synteny run"
        if self.args.n == 0:
            self.args.n = len(self.args.FILES)
        print("Running ntJoin synteny detection...")
        print("Parameters:")
        print("\tMinimizer TSV files: ", self.args.FILES)
        print("\t-n", self.args.n)
        print("\t-p", self.args.p)
        print("\t-k", self.args.k)
        print("\t-w", self.args.w)
        print("\t--btllib_t", self.args.btllib_t)
        print("\t--w-rounds", self.args.w_rounds)


    def find_synteny_blocks(self, path):
        "Given a path (sequence of mx), print the order/orientation/regions of contigs for an assembly"
        out_blocks = []  # List of SyntenyBlock
        prelim_blocks = SyntenyBlock(self.args.k, *list(self.list_mx_info.keys()))
        past_start_flag = False
        for mx in path:
            if prelim_blocks.continue_block(mx, self.list_mx_info):
                prelim_blocks.extend_block(mx, self.list_mx_info)
            else:
                # This is either the first mx, or we are past a stretch of repeating contigs
                if past_start_flag:
                    prelim_blocks.determine_orientations()
                    if prelim_blocks.all_oriented():
                        out_blocks.append(prelim_blocks)
                prelim_blocks = SyntenyBlock(self.args.k, *list(self.list_mx_info.keys()))
                prelim_blocks.start_block(mx, self.list_mx_info)

        prelim_blocks.determine_orientations()
        if prelim_blocks.all_oriented():
            out_blocks.append(prelim_blocks)

        return out_blocks

    @staticmethod
    def find_fa_name(assembly_mx_name):
        "Given the mx file name, return the corresponding fai file name"
        if fai_match := re.search(fai_re, assembly_mx_name):
            return f"{fai_match.group(1)}"
        print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
        print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
        sys.exit(1)

    @staticmethod
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

    def mask_assemblies_with_synteny_extents(self, synteny_beds, w):
        "Mask each reference assembly with determined synteny blocks"
        mx_to_fa_dict = {}
        for assembly, contig_dict in synteny_beds.items():
            bed_str = [f"{ctg}\t{bed.start}\t{bed.end}\tSYNTENY" for ctg in contig_dict \
                        for bed in contig_dict[ctg] if bed.end - bed.start > 2*w]
            bed_str = "\n".join(bed_str)
            fa_filename = self.find_fa_name(assembly)
            synteny_bed = pybedtools.BedTool(bed_str, from_string=True).slop(g=f"{fa_filename}.fai",
                                                                             l=-1*w, r=-1*w).sort()
            synteny_bed.mask_fasta(fi=fa_filename, fo=f"{fa_filename}_masked.fa")
            mx_to_fa_dict[assembly] = f"{fa_filename}_masked.fa"
        return mx_to_fa_dict

    @staticmethod
    def delete_w_iteration_files(*filenames):
        "Delete the given files for the specific lower w iteration"
        for filename in filenames:
            cmd = shlex.split(f"rm {filename}")
            ret_code = subprocess.call(cmd)
            assert ret_code == 0

    def generate_new_minimizers(self, tsv_to_fa_dict, w, retain_files=False):
        "Given the masked fasta files, generate minimizers at new w for each"
        list_mxs = {}
        new_list_mxs_info = {}
        for assembly_tsv, assembly_masked in tsv_to_fa_dict.items():
            indexlr_filename = ntjoin_utils.run_indexlr(assembly_masked, self.args.k, w, self.args.btllib_t)
            mx_info, mxs_filt = ntjoin_utils.read_minimizers(indexlr_filename)
            new_list_mxs_info[assembly_tsv] = mx_info
            list_mxs[assembly_tsv] = mxs_filt
            if not retain_files:
                self.delete_w_iteration_files(indexlr_filename, assembly_masked)
        return list_mxs, new_list_mxs_info

    @staticmethod
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

    def find_mx_in_blocks(self, paths):
        "Given the synteny blocks, find the minimizers at the terminal ends of each block, and internal mxs"
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
                    self.update_interval_tree(intervaltrees, assembly, contig, mx1, mx2)
                    internal = assembly_block.get_block_internal_mx_hashes()
                    internal_mxs = internal_mxs.union(internal)
                assert len(terminal_mxs) == (curr_mx_len + 2)
        return terminal_mxs, internal_mxs, intervaltrees

    def check_non_overlapping(self, paths):
        "Given the paths, do final check to ensure intervals are not overlapping, will print warnings if so"
        intervaltrees = defaultdict(dict) # assembly -> contig -> IntervalTree of synteny block extents
        for subcomponent in paths:
            for block in subcomponent:
                for assembly, assembly_block in block.assembly_blocks.items():
                    contig, mx1, mx2 = assembly_block.get_block_terminal_mx()
                    self.update_interval_tree(intervaltrees, assembly, contig, mx1, mx2)

    @staticmethod
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

    @staticmethod
    def update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info):
        "Update the directory containing mx -> contig, position associations"
        valid_mxs = set({mx for _, list_mx_val in list_mxs.items() \
                        for list_mx in list_mx_val for mx in list_mx})
        for assembly, mx_dict in new_list_mx_info.items():
            for mx in mx_dict:
                if mx in valid_mxs and mx not in list_mx_info[assembly]:
                    list_mx_info[assembly][mx] = mx_dict[mx]

    def refine_block_coordinates(self, paths):
        "Ready to start refining the synteny block coordinates"
        prev_w = self.args.w
        for new_w in self.args.w_rounds:
            print(datetime.datetime.today(), ": Extending synteny blocks with w =", new_w, file=sys.stdout)
            new_list_mxs, terminal_mxs = self.generate_additional_minimizers(
                    paths, new_w, prev_w, self.list_mx_info, self.args.dev)
            graph = ntjoin_utils.build_graph(new_list_mxs, self.weights, graph=self.graph, black_list=terminal_mxs)
            graph = self.filter_graph_global(graph)
            paths = self.find_paths_synteny_blocks(self.find_paths())
            with open(f"{self.args.p}.synteny_blocks.extended.tsv", 'w', encoding="utf-8") as outfile:
                block_num = 0
                for subcomponent in paths:
                    for block in subcomponent:
                        outfile.write(block.get_block_string(block_num))
                        block_num += 1
            prev_w = new_w

        print(datetime.datetime.today(), ": Done extended synteny blocks", file=sys.stdout)
        with open(f"{self.args.p}.synteny_blocks.extended.tsv", 'w', encoding="utf-8") as outfile:
            block_num = 0
            for subcomponent in paths:
                for block in subcomponent:
                    outfile.write(block.get_block_string(block_num))
                    block_num += 1

    def generate_additional_minimizers(self, paths, new_w, prev_w, list_mx_info, dev=False):
        "Given the existing synteny blocks, generate minimizers for increased block resolution"
        synteny_beds = self.get_synteny_bed_lists(paths)
        mx_to_fa_dict = self.mask_assemblies_with_synteny_extents(synteny_beds, prev_w)
        list_mxs, new_list_mx_info = self.generate_new_minimizers(mx_to_fa_dict, new_w, retain_files=dev)
        terminal_mx, internal_mx, interval_trees = self.find_mx_in_blocks(paths)
        list_mxs = self.filter_minimizers_synteny_blocks(list_mxs, internal_mx, interval_trees, new_list_mx_info)
        list_mxs = ntjoin_utils.filter_minimizers(list_mxs) # Filter for mx in all assemblies
        self.update_list_mx_info(list_mxs, list_mx_info, new_list_mx_info)
        return list_mxs, terminal_mx

    def find_paths_synteny_blocks(self, paths):
        "Given a list of paths, return a list of representative synteny blocks"
        return [self.find_synteny_blocks(blocks) for path in paths for blocks, _ in path]

    def load_minimizers(self):
        "Read the minimizers for the synteny mode"
        weights = {}  # Dictionary: assembly -> weight
        for assembly in self.args.FILES:
            mxs_info, mxs = ntjoin_utils.read_minimizers(assembly)
            self.list_mx_info[assembly] = mxs_info
            self.list_mxs[assembly] = mxs
            weights[assembly] = self.weights_list.pop(0)
        self.weights = weights

    def main_synteny(self):
        "Run the steps for ntJoin synteny mode"
        print("Running ntJoin synteny detection ...\n", file=sys.stdout)

        # Run the common ntJoin steps
        self.load_minimizers()

        self.make_minimizer_graph_and_paths()

        paths = self.ntjoin_find_paths()

        paths = self.find_paths_synteny_blocks(paths)

        with open(f"{self.args.p}.synteny_blocks.tsv", 'w', encoding="utf-8") as outfile:
            block_num = 0
            for subcomponent in paths:
                for block in subcomponent:
                    outfile.write(block.get_block_string(block_num))
                    block_num += 1
        print(datetime.datetime.today(), ": Done initial synteny blocks", file=sys.stdout)
        self.refine_block_coordinates(paths)

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)
