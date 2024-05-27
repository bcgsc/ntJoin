#!/usr/bin/env python3
"""
ntJoin: Scaffolding assemblies using reference assemblies and minimizer graphs
Written by Lauren Coombe (@lcoombe)
"""

import datetime
import re
from collections import Counter, defaultdict
import shlex
import subprocess
import sys
import warnings
import pybedtools
import pymannkendall as mk
import btllib
from packaging import version
from read_fasta import read_fasta
import ntjoin_utils
import ntjoin_overlap
import ntjoin
from path_node import PathNode
from overlap_region import OverlapRegion
warnings.simplefilter(action='ignore', category=RuntimeWarning)


class NtjoinScaffolder(ntjoin.Ntjoin):
    "ntJoin: Scaffolding assemblies using reference assemblies and minimizer graphs"

    def determine_orientation(self, positions):
        "Given a list of minimizer positions, determine the orientation of the contig"
        if len(positions) > 1:
            if all(x < y for x, y in zip(positions, positions[1:])):
                return "+"
            if all(x > y for x, y in zip(positions, positions[1:])):
                return "-"
            if self.args.mkt:
                mkt_result = mk.original_test(positions)
                if mkt_result.h and mkt_result.p <= 0.05:
                    return "+" if mkt_result.trend == "increasing" else "-"
            else:
                tally = Counter([x < y for x, y in zip(positions, positions[1:])])
                positive_perc = tally[True]/float(len(positions)-1)*100
                negative_perc = 100 - positive_perc
                if positive_perc >= self.args.m:
                    return "+"
                if negative_perc >= self.args.m:
                    return "-"

        return "?"

    @staticmethod
    def calc_start_coord(positions, ctg_min_mx):
        "Calculates the minimum coordinate for a contig region in a path"
        if min(positions) == ctg_min_mx:
            return 0
        return min(positions)


    def calc_end_coord(self, positions, ctg_max_mx, ctg_len):
        "Calculates the maximum coordinate for a contig region in a path"
        if max(positions) == ctg_max_mx:
            return ctg_len
        return max(positions) + self.args.k


    def calculate_gap_size(self, u, v, graph, cur_assembly):
        "Calculates the estimated distance between two contigs"
        u_mx = u.terminal_mx
        v_mx = v.first_mx

        # Don't attempt gap estimation when don't know orientation
        if u.ori == "?" or v.ori == "?":
            return 0, 0

        # Find the assemblies that have a path between these mx
        # Are situations where there is not a direct edge if an unoriented contig was in-between
        path = graph.get_shortest_paths(u_mx, v_mx, output="vpath")[0]
        supporting_assemblies = set.intersection(
            *map(set, [graph.es()[ntjoin_utils.edge_index(graph, s, t)]['support']
                       for s, t in zip(path, path[1:])]))
        if not supporting_assemblies:
            return self.args.g, self.args.g

        distances = [abs(self.list_mx_info[assembly][v_mx][1] - self.list_mx_info[assembly][u_mx][1])
                     for assembly in supporting_assemblies]
        mean_dist = int(sum(distances)/len(distances)) - self.args.k
        # Correct for the overhanging sequence before/after terminal minimizers
        if u.ori == "+":
            a = u.end - self.list_mx_info[cur_assembly][u_mx][1] - self.args.k
        else:
            a = self.list_mx_info[cur_assembly][u_mx][1] - u.start
        if v.ori == "+":
            b = self.list_mx_info[cur_assembly][v_mx][1] - v.start
        else:
            b = v.end - self.list_mx_info[cur_assembly][v_mx][1] - self.args.k

        try:
            assert a >= 0
            assert b >= 0
        except AssertionError as assert_error:
            print("ERROR: Gap distance estimation less than 0", "Vertex 1:", u, "Vertex 2:", v,
                  sep="\n")
            print("Minimizer positions:", self.list_mx_info[cur_assembly][u_mx][1],
                  self.list_mx_info[cur_assembly][v_mx][1])
            print("Estimated distance: ", mean_dist)
            raise ValueError from assert_error

        gap_size = max(mean_dist - a - b, self.args.g)
        if self.args.G > 0:
            gap_size = min(gap_size, self.args.G)

        return gap_size, mean_dist - a - b

    @staticmethod
    def is_new_region_overlapping(start, end, node_i, node_j, incorporated_segments_ctg):
        "Checks if the specified region overlaps any existing regions in incorporated segments"
        for segment in incorporated_segments_ctg:
            if start <= segment.end and segment.start <= end and \
                    (segment.start != node_i.start and segment.end != node_i.end) and \
                    (segment.start != node_j.start and segment.end != node_j.end):
                return True
        return False


    def merge_relocations(self, path, incorporated_segments):
        "If a path has adjacent collinear intervals of the same contig, merge them"
        if len(path) < 2:
            return path
        return_path = [path[0]]
        for node_i, node_j in zip(path, path[1:]):
            if node_i.contig == node_j.contig:
                if node_i.ori == "+" and node_j.ori == "+" and node_i.end <= node_j.start:
                    if self.is_new_region_overlapping(node_i.start, node_j.end, node_i, node_j,
                                                      incorporated_segments[node_i.contig]):
                        return_path.append(node_j)
                        continue
                    incorporated_segments[node_i.contig].add(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                        start=return_path[-1].start,
                                                                        end=node_j.end))
                    incorporated_segments[node_i.contig].remove(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                           start=return_path[-1].start,
                                                                           end=return_path[-1].end))
                    incorporated_segments[node_j.contig].remove(ntjoin_utils.Bed(contig=node_j.contig,
                                                                           start=node_j.start,
                                                                           end=node_j.end))
                    return_path[-1].end = node_j.end
                    return_path[-1].terminal_mx = node_j.terminal_mx
                    return_path[-1].gap_size = node_j.gap_size
                elif node_i.ori == "-" and node_j.ori == "-" and node_i.start >= node_j.end:
                    if self.is_new_region_overlapping(node_j.start, node_i.end, node_i, node_j,
                                                      incorporated_segments[node_i.contig]):
                        return_path.append(node_j)
                        continue
                    incorporated_segments[node_i.contig].add(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                        start=node_j.start,
                                                                        end=return_path[-1].end))
                    incorporated_segments[node_i.contig].remove(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                           start=return_path[-1].start,
                                                                           end=return_path[-1].end))
                    incorporated_segments[node_j.contig].remove(ntjoin_utils.Bed(contig=node_j.contig,
                                                                           start=node_j.start,
                                                                           end=node_j.end))
                    return_path[-1].start = node_j.start
                    return_path[-1].first_mx = node_j.first_mx
                    return_path[-1].gap_size = node_j.gap_size
                else:
                    return_path.append(node_j)
            else:
                return_path.append(node_j)

        return return_path


    def format_path(self, path, assembly, component_graph):
        "Given a path (sequence of mx), print the order/orientation/regions of contigs for an assembly"
        out_path = []  # List of PathNode
        curr_ctg, prev_mx, first_mx = None, None, None
        positions = []
        for mx in path:
            ctg, pos = self.list_mx_info[assembly][mx]
            if ctg is curr_ctg:
                positions.append(pos)
            else:
                # This is either the first mx, or we are past a stretch of repeating contigs
                if curr_ctg is not None:
                    ori = self.determine_orientation(positions)
                    if ori != "?":  # Don't add to path if orientation couldn't be determined
                        out_path.append(PathNode(contig=curr_ctg, ori=ori,
                                                 start=self.calc_start_coord(positions,
                                                                             self.mx_extremes[curr_ctg][0]),
                                                 end=self.calc_end_coord(positions,
                                                                         self.mx_extremes[curr_ctg][1],
                                                                         self.scaffolds[curr_ctg].length),
                                                 contig_size=self.scaffolds[curr_ctg].length,
                                                 first_mx=first_mx,
                                                 terminal_mx=prev_mx))
                curr_ctg = ctg
                positions = [pos]
                first_mx = mx
            prev_mx = mx
        ori = self.determine_orientation(positions)
        if ori != "?":
            out_path.append(PathNode(contig=curr_ctg, ori=ori,
                                     start=self.calc_start_coord(positions,
                                                                 self.mx_extremes[curr_ctg][0]),
                                     end=self.calc_end_coord(positions,
                                                             self.mx_extremes[curr_ctg][1],
                                                             self.scaffolds[curr_ctg].length),
                                     contig_size=self.scaffolds[curr_ctg].length,
                                     first_mx=first_mx,
                                     terminal_mx=prev_mx))
        for u, v in zip(out_path, out_path[1:]):
            gap_size, raw_gap_size = self.calculate_gap_size(u, v, component_graph, assembly)
            u.set_gap_size(gap_size)
            u.set_raw_gap_size(raw_gap_size)

        return out_path

    @staticmethod
    def tally_incorporated_segments(incorporated_list, path):
        "Keep track of contig segments incorporated into path"
        if len(path) < 2:
            return
        for path_node in path:
            if path_node.contig not in incorporated_list:
                incorporated_list[path_node.contig] = set()
            incorporated_list[path_node.contig].add(ntjoin_utils.Bed(contig=path_node.contig,
                                                        start=path_node.start,
                                                        end=path_node.end))


    @staticmethod
    def is_best_region(path_nodes, query_node):
        "Given query node and list of nodes with the same contig ID, return True if query node is the largest region"
        max_len = 0
        max_node = None
        for node in path_nodes:
            if node.get_aligned_length() > max_len:
                max_len = node.get_aligned_length()
                max_node = node
        if query_node.get_aligned_length() == max_len and max_node.terminal_mx == query_node.terminal_mx:
            return True
        return False


    @staticmethod
    def is_node_full_sequence(node, scaffold):
        "Given a Path node and a scaffold, return True if that path contains the entire sequence of the scaffold"
        return node.get_aligned_length() >= scaffold.length


    @staticmethod
    def is_subsumed(i, path, contig_regions):
        "Returns True if a contig is subsumed"
        if i == 0 or i >= (len(path)-1):
            return False
        (prev_node, next_node) = path[i-1], path[i+1]
        if prev_node.contig == next_node.contig and prev_node.ori == next_node.ori and\
            min(prev_node.start, next_node.start) == 0 and\
                max(prev_node.end, next_node.end) == prev_node.contig_size and\
                len(contig_regions[prev_node.contig]) == 2:
            return True
        return False

    def adjust_paths(self, paths, scaffolds, incorporated_segments):
        "Given the found paths, removes duplicate regions to avoid cutting sequences (no_cut=True option)"
        contig_regions = {}  # contig_id -> [list of PathNode]
        for path in paths:
            for node in path:
                if node.contig not in contig_regions:
                    contig_regions[node.contig] = []
                contig_regions[node.contig].append(node)

        intermediate_paths = []
        for path in paths:
            new_path = []
            for i, node in enumerate(path):
                if not self.is_subsumed(i, path, contig_regions):
                    new_path.append(node)
            new_path = self.merge_relocations(new_path, incorporated_segments)
            intermediate_paths.append(new_path)

        new_paths = []
        for path in intermediate_paths:
            new_path = []
            for i, node in enumerate(path):
                if (len(contig_regions[node.contig]) > 1 \
                    and self.is_best_region(contig_regions[node.contig], node)) \
                        or len(contig_regions[node.contig]) == 1 \
                        and not self.is_node_full_sequence(node, scaffolds[node.contig]):
                    node.start = 0
                    node.end = scaffolds[node.contig].length
                    new_path.append(node)
                elif len(contig_regions[node.contig]) > 1 \
                        and not self.is_best_region(contig_regions[node.contig], node):
                    if 0 < i < len(path)-1 and new_path:
                        new_path[-1].gap_size += (node.get_aligned_length())
                        if self.args.G > 0:
                            new_path[-1].gap_size = min(self.args.G, new_path[-1].gap_size)
                else:
                    new_path.append(node)

            new_paths.append(new_path)
        return new_paths


    def read_fasta_file(self, filename):
        "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
        print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
        scaffolds = {}
        try:
            with btllib.SeqReader(filename, btllib.SeqReaderFlag.LONG_MODE,
                                  self.args.btllib_t) as fin:
                for rec in fin:
                    scaffolds[rec.id] = ntjoin_utils.Scaffold(id=rec.id, length=len(rec.seq), sequence=rec.seq)
        except FileNotFoundError:
            print("ERROR: File", filename, "not found.")
            print("Minimizer TSV file must follow the naming convention:")
            print("\tassembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering,\n"
                  "and assembly.fa is the scaffolds fasta file")
            sys.exit(1)
        return scaffolds


    @staticmethod
    def get_fasta_segment(path_node, sequence):
        "Given a PathNode and the contig sequence, return the corresponding sequence"
        if path_node.ori == "-":
            return ntjoin_utils.reverse_complement(sequence[path_node.start:path_node.end]) + \
                   "N"*path_node.gap_size
        return sequence[path_node.start:path_node.end] + "N"*path_node.gap_size

    @staticmethod
    def format_bedtools_genome(scaffolds):
        "Format a BED file and genome dictionary for bedtools"
        bed_str = "\n".join([f"{scaffold}\t{0}\t{scaffolds[scaffold].length}"
                             for scaffold in scaffolds])
        bed = pybedtools.BedTool(bed_str, from_string=True)

        genome_dict = {scaffold: (0, scaffolds[scaffold].length) for scaffold in scaffolds}

        return bed, genome_dict

    @staticmethod
    def write_agp(agp_file, scaffold_id, path_str):
        "Write the given path string in AGP format to file"
        contig_re = re.compile(r'(\S+)([\+\-])\:(\d+)-(\d+)')
        gap_re = re.compile(r'(\d+)N')
        format_layout = ("{}\t" * 9).strip()

        ctg_coord_start = 1
        ctg_part_count = 1
        for component in path_str.split():
            contig_match = re.search(contig_re, component)
            gap_match = re.search(gap_re, component)
            if contig_match:
                (contig_id, contig_ori, contig_start, contig_end) = (contig_match.group(1), contig_match.group(2),
                                                                     int(contig_match.group(3))+1,
                                                                     int(contig_match.group(4)))

                length_segment = contig_end - contig_start + 1
                out_line = format_layout.format(scaffold_id, ctg_coord_start, ctg_coord_start + length_segment - 1,
                                                ctg_part_count, "W",
                                                contig_id, contig_start, contig_end, contig_ori)
            elif gap_match:
                length_segment = int(gap_match.group(1))
                out_line = format_layout.format(scaffold_id, ctg_coord_start, ctg_coord_start + length_segment - 1,
                                                ctg_part_count, "N",
                                                length_segment, "scaffold", "yes", "align_genus")
            else:
                print("ERROR: Path string is not formatted correctly: " + path_str)
                sys.exit(1)
            agp_file.write(out_line + "\n")
            ctg_coord_start = ctg_coord_start + length_segment
            ctg_part_count += 1

    @staticmethod
    def write_agp_unassigned(agpfile, header, seq):
        "Write unassigned contig to AGP file"
        header_re = re.compile(r'((\S+)\:(\d+)-(\d+))')
        format_layout = ("{}\t" * 9).strip()

        len_diff_start, len_diff_end = 0, 0
        sequence_start_strip = seq.strip().lstrip("Nn") # Strip from 5'
        if len(sequence_start_strip) != len(seq):
            len_diff_start = len(seq) - len(sequence_start_strip)
        sequence_end_strip = sequence_start_strip.rstrip("Nn") # Strip from 3'
        if len(sequence_end_strip) != len(sequence_start_strip):
            len_diff_end = len(sequence_start_strip) - len(sequence_end_strip)

        if not sequence_end_strip: # Just return if the sequence was all N (so empty)
            return

        header_match = re.search(header_re, header)
        agp = None
        if header_match:
            agp = ntjoin_utils.Agp(new_id=header_match.group(1), contig=header_match.group(2),
                      start=int(header_match.group(3)) + 1 + len_diff_start,
                      end=int(header_match.group(4)) - len_diff_end)
        assert len(seq.strip().strip("Nn")) == agp.end - agp.start + 1
        out_str = format_layout.format(agp.new_id, 1, agp.end - agp.start + 1,
                                       1, "W", agp.contig, agp.start, agp.end, "+")
        agpfile.write(out_str + "\n")

    @staticmethod
    def join_sequences(sequences_list, path, path_segments):
        "Join the sequences for a contig, adjusting the path coordinates if terminal Ns are stripped"
        sequence_start_strip = sequences_list[0].lstrip("Nn") # Strip from 5'
        if len(sequence_start_strip) != len(sequences_list[0]):
            len_diff = len(sequences_list[0]) - len(sequence_start_strip)
            sequences_list[0] = sequence_start_strip
            for i, node in enumerate(path):
                if node.contig == path_segments[0].contig and \
                                node.start == path_segments[0].start and \
                                node.end == path_segments[0].end:
                    if node.ori == "+":
                        node.start += len_diff
                    else:
                        node.end -= len_diff
                    assert len(sequence_start_strip) - node.gap_size == node.end - node.start
                    break

        sequence_end_strip = sequences_list[-1].rstrip("Nn") # Strip from 3'
        if len(sequence_end_strip) != len(sequences_list[-1]):
            len_diff = len(sequences_list[-1]) - len(sequence_end_strip)
            sequences_list[-1] = sequence_end_strip
            for i in reversed(range(len(path))):
                if path[i].contig == path_segments[-1].contig and \
                                path[i].start == path_segments[-1].start and \
                                path[i].end == path_segments[-1].end:
                    if path[i].ori == "+":
                        path[i].end -= len_diff
                    else:
                        path[i].start += len_diff
                    assert len(sequence_end_strip) == path[i].end - path[i].start
                    break

        return "".join(sequences_list)

    @staticmethod
    def check_terminal_node_gap_zero(path):
        "Ensure that the terminal PathNode has gap size of 0"
        for i in reversed(range(len(path))):
            if path[i].ori != "?":
                if path[i].gap_size != 0:
                    path[i].set_gap_size(0)
                break

    @staticmethod
    def remove_overlapping_regions(path, intersecting_regions):
        "Remove any regions that are overlapping, adjusting if needed"
        new_path = []
        for path_node in path:
            if path_node.contig in intersecting_regions:
                node_bed = ntjoin_utils.Bed(contig=path_node.contig, start=path_node.start, end=path_node.end)
                if node_bed in intersecting_regions[path_node.contig]:
                    new_bed = intersecting_regions[path_node.contig][node_bed]
                    if new_bed is None:
                        continue
                    if new_bed != node_bed:
                        path_node.start = new_bed.start
                        path_node.end = new_bed.end
            new_path.append(path_node)

        return new_path

    def adjust_for_trimming(self, fasta_filename, paths):
        "Go through path, trim the segments if overlapping"
        ct = 0
        mx_info = defaultdict(dict)  # path_index -> mx -> pos
        mxs = {}  # path_index -> [mx]
        cur_path_index = 0
        if not paths:
            return # If list of paths is empty
        cur_valid_segments = {f"{node.contig}_{node.start}_{node.end}"
                                  for node in paths[cur_path_index]}
        with btllib.Indexlr(fasta_filename, self.args.overlap_k, self.args.overlap_w,
                            btllib.IndexlrFlag.LONG_MODE, self.args.btllib_t) as minimizers:
            for mx_entry in minimizers:
                if mx_entry.id in cur_valid_segments:
                    self.tally_minimizers_overlap(ct, cur_path_index, mx_entry, mx_info, mxs, paths)
                    ct += 1
                else:
                    assert len(mx_info) == len(paths[cur_path_index])
                    ntjoin_overlap.merge_overlapping_path(paths[cur_path_index], mxs, mx_info)

                    ct = 0
                    mx_info = defaultdict(dict)  # path_index -> mx -> pos
                    mxs = {}  # path_index -> [mx]
                    cur_path_index += 1
                    cur_valid_segments = {f"{node.contig}_{node.start}_{node.end}"
                                              for node in paths[cur_path_index]}

                    if mx_entry.id in cur_valid_segments:
                        self.tally_minimizers_overlap(ct, cur_path_index, mx_entry, mx_info, mxs, paths)
                        ct += 1
        # Don't miss last path
        ntjoin_overlap.merge_overlapping_path(paths[cur_path_index], mxs, mx_info)

    @staticmethod
    def tally_minimizers_overlap(ct, cur_path_index, mx_entry, mx_info, mxs, paths):
        "Tally minimizer info for the given path"
        mxs[ct] = []
        dup_mxs = set()  # Set of minimizers identified as duplicates
        for mx_pos_strand in mx_entry.minimizers:
            mx, pos = str(mx_pos_strand.out_hash), mx_pos_strand.pos
            if not ntjoin_overlap.is_in_valid_region(pos, ct, paths[cur_path_index]):
                continue
            if ct in mx_info and mx in mx_info[ct]:  # This is a duplicate
                dup_mxs.add(mx)
            else:
                mx_info[ct][mx] = int(pos)
                mxs[ct].append(mx)
        mx_info[ct] = {mx: mx_info[ct][mx] for mx in mx_info[ct] if mx not in dup_mxs}
        mxs[ct] = [[mx for mx in mxs[ct] if mx not in dup_mxs and mx in mx_info[ct]]]


    def get_adjusted_sequence(self, sequence, node):
        "Return sequence adjusted for overlap trimming"
        return_sequence = sequence[node.start_adjust:node.get_end_adjusted_coordinate()]
        if node.gap_size > 0:
            if node.get_end_adjusted_coordinate() == node.get_aligned_length():
                # no trimming was done on the end, keep calculated gap size
                return return_sequence + "N"*node.gap_size
            return return_sequence + self.args.overlap_gap*"N"
        return return_sequence


    def print_scaffolds(self, paths, intersecting_regions, prev_incorporated_segments):
        "Given the paths, print out the scaffolds fasta"
        print(datetime.datetime.today(), ": Printing output scaffolds", file=sys.stdout)
        assembly = self.args.s

        min_match = re.search(r'^(\S+)(.k\d+.w\d+)\.tsv', assembly)
        assembly_fa, params = min_match.group(1), min_match.group(2)

        outfile = open(assembly_fa + params + ".n" + # pylint: disable=consider-using-with
                       str(self.args.n) + ".assigned.scaffolds.fa", 'w', encoding="utf-8")
        pathfile = open(self.args.p + ".path", 'w', encoding="utf-8") # pylint: disable=consider-using-with
        if self.args.agp:
            agpfile = open(self.args.p + ".agp", "w", encoding="utf-8") # pylint: disable=consider-using-with
        incorporated_segments = []  # List of Bed entries

        ct = 0
        pathfile.write(assembly_fa + "\n")

        # Deal with merging relocations
        for i, path in enumerate(paths):
            new_path = self.merge_relocations(path, prev_incorporated_segments)
            new_path = self.remove_overlapping_regions(new_path, intersecting_regions)
            self.check_terminal_node_gap_zero(new_path)
            paths[i] = new_path

        # Trim overlaps if option turned on
        if self.args.overlap:
            path_segments_file = open(self.args.p + ".segments.fa", 'w', encoding="utf-8") # pylint: disable=consider-using-with
            filtered_paths = []
            for path in paths:
                sequences = []
                nodes = []
                for node in path:
                    if node.ori == "?":
                        continue
                    sequences.append(self.get_fasta_segment(node, self.scaffolds[node.contig].sequence))
                    nodes.append(node)
                if len(sequences) < 2:
                    continue
                out_coords = ntjoin_overlap.get_valid_regions(nodes, self.args.overlap_k, self.args.overlap_w)
                for seq, node, out_coords in zip(sequences, nodes, out_coords):
                    my_seq = seq.strip("Nn")
                    print(f"my_seq: {my_seq} ({len(my_seq)}bp)")
                    print(f"out_coords: {out_coords}")
                    my_seq = my_seq[:out_coords[0]] + "N"*(out_coords[1] - out_coords[0]) + my_seq[out_coords[1]:]
                    print(f"my_seq: {my_seq} ({len(my_seq)}bp) ; node: {node.get_aligned_length()}")
                    print("-------")
                    assert len(my_seq) == node.get_aligned_length()
                    path_segments_file.write(f">{node.contig}_{node.start}_{node.end} { node.raw_gap_size}\n{my_seq}\n")
                filtered_paths.append(nodes)
            path_segments_file.close()

            self.adjust_for_trimming(self.args.p + ".segments.fa", filtered_paths)

        for path in paths:
            sequences = []
            path_segments = []
            nodes = []

            for node in path:
                if node.ori == "?":
                    continue
                sequences.append(self.get_fasta_segment(node, self.scaffolds[node.contig].sequence))
                path_segments.append(ntjoin_utils.Bed(contig=node.contig, start=node.start,
                                                      end=node.end))
                nodes.append(node)

            if len(sequences) < 2:
                continue

            if self.args.overlap:
                sequences = [self.get_adjusted_sequence(sequence, nodes[i])
                             for i, sequence in enumerate(sequences)]

            ctg_id = "ntJoin" + str(ct)
            ctg_sequence = self.join_sequences(sequences, path, path_segments)

            outfile.write(f">{ctg_id}\n{ctg_sequence}\n")
            incorporated_segments.extend(path_segments)
            path_str = " ".join([f"{node.contig}{node.ori}:{node.get_adjusted_start()}-" \
                            f"{node.get_adjusted_end()} {node.gap_size}N" for node in path])
            path_str = re.sub(r'\s+\d+N$', r'', path_str)
            pathfile.write(f"{ctg_id}\t{path_str}\n")
            if self.args.agp:
                self.write_agp(agpfile, ctg_id, path_str)

            ct += 1
        outfile.close()

        if self.args.agp:
            self.print_unassigned(assembly, assembly_fa, incorporated_segments, params,
                                            agpfile=agpfile)
        else:
            self.print_unassigned(assembly, assembly_fa, incorporated_segments, params)

        pathfile.close()
        if self.args.agp:
            agpfile.close()
        if self.args.overlap:
            cmd_shlex = shlex.split(f"rm {self.args.p}.segments.fa")
            subprocess.call(cmd_shlex)

    def print_unassigned(self, assembly, assembly_fa, incorporated_segments, params, agpfile=None):
        "Also print out the sequences that were NOT scaffolded"
        incorporated_segments_str = "\n".join([f"{chrom}\t{s}\t{e}"
                                               for chrom, s, e in incorporated_segments])
        genome_bed, genome_dict = self.format_bedtools_genome(self.scaffolds)
        # Needed to deal with failure in complement step seen with pybedtools 0.9.1+
        if version.parse(pybedtools.__version__) < version.parse("0.9.1"):
            incorporated_segments_bed = pybedtools.BedTool(incorporated_segments_str,
                                                           from_string=True).sort()
        else:
            incorporated_segments_bed = pybedtools.BedTool(incorporated_segments_str,
                                                           from_string=True).sort(genome=genome_dict)
        missing_bed = genome_bed.complement(i=incorporated_segments_bed, g=genome_dict)
        missing_bed.saveas(self.args.p + "." + assembly + ".unassigned.bed")

        with open(f"{assembly_fa}{params}.n{str(self.args.n)}.unassigned.scaffolds.fa",
                 'w', encoding="utf-8") as outfile:
            cmd = f"bedtools getfasta -fi {assembly_fa} -bed {self.args.p}.{assembly}.unassigned.bed -fo -"
            cmd_shlex = shlex.split(cmd)
            out_fasta = subprocess.Popen(cmd_shlex, stdout=subprocess.PIPE, universal_newlines=True)
            for header, seq, _, _ in read_fasta(iter(out_fasta.stdout.readline, '')):
                if self.args.agp:
                    self.write_agp_unassigned(agpfile, header, seq)
                seq = seq.strip().strip("Nn")
                if seq:
                    outfile.write(f">{header}\n{seq}\n")
            out_fasta.wait()
            if out_fasta.returncode != 0:
                print("bedtools getfasta failed -- is bedtools on your PATH?")
                print(out_fasta.stderr)
                raise subprocess.CalledProcessError(out_fasta.returncode, cmd_shlex)

    @staticmethod
    def tally_intersecting_segments(incorporated_segments):
        "Tally ctgs with intersecting segments, and keep track of 'best'"
        incorporated_bed_list = []
        for _, bed_entry_list in incorporated_segments.items():
            for bed_entry in bed_entry_list:
                incorporated_bed_list.append(bed_entry)
        incorporated_bed_str = "\n".join([f"{chrom}\t{s}\t{e}"
                                          for chrom, s, e in incorporated_bed_list])
        incorporated_segments_bed = pybedtools.BedTool(incorporated_bed_str,
                                                       from_string=True).sort()
        bed_intersect = incorporated_segments_bed.intersect(b=incorporated_segments_bed,
                                                            c=True, wa=True)

        overlap_regions = {}

        for bed in bed_intersect:
            if bed.count > 1:
                if bed.chrom not in overlap_regions:
                    overlap_regions[bed.chrom] = OverlapRegion()
                overlap_regions[bed.chrom].add_region(ntjoin_utils.Bed(contig=bed.chrom, start=bed.start, end=bed.end))

        overlap_regions_fix = {}
        for overlap_contig, bed_region in overlap_regions.items():
            overlap_regions_fix[overlap_contig] = bed_region.find_non_overlapping()

        return overlap_regions_fix

    def find_mx_min_max(self, target):
        "Given a dictionary in the form mx->(ctg, pos), find the min/max mx position per ctg"
        mx_extremes = {} # ctg -> (min_pos, max_pos)
        for mx in self.list_mx_info[target]:
            try:
                self.graph.vs().find(mx)
            except ValueError:
                continue
            ctg, pos = self.list_mx_info[target][mx]
            if ctg in mx_extremes:
                mx_extremes[ctg] = (min(mx_extremes[ctg][0], pos),
                                    max(mx_extremes[ctg][1], pos))
            else:
                mx_extremes[ctg] = (pos, pos)
        return mx_extremes

    def format_adjust_paths(self, paths):
        "Format and adjust the paths for relocations, incorporated sections"
        return_paths = []
        incorporated_segments = {}
        for path_list in paths:
            for path, sub_graph in path_list:
                ctg_path = self.format_path(path, self.args.s, sub_graph)
                return_paths.append(ctg_path)
                self.tally_incorporated_segments(incorporated_segments, ctg_path)

        paths_return_merged = []
        for path in return_paths:
            path = self.merge_relocations(path, incorporated_segments)
            paths_return_merged.append(path)

        return paths_return_merged, incorporated_segments


    def print_parameters_scaffold(self):
        "Print the set parameters for the ntJoin scaffolding run"
        print("Running ntJoin scaffolding..")
        print("Parameters:")
        print("\tReference TSV files: ", self.args.FILES)
        print("\t-s ", self.args.s)
        print("\t-l ", self.args.l)
        print("\t-r ", self.args.r)
        print("\t-p ", self.args.p)
        print("\t-n ", self.args.n)
        print("\t-k ", self.args.k)
        print("\t-g ", self.args.g)
        print("\t-G ", self.args.G)
        print("\t-t ", self.args.t)
        if self.args.agp:
            print("\t--agp")
        if self.args.no_cut:
            print("\t--no_cut")
        if self.args.mkt:
            print("Orienting contigs with Mann-Kendall Test (more computationally intensive)\n")
        else:
            print("Orienting contigs using increasing/decreasing minimizer positions\n")
        if self.args.overlap:
            print("\t--overlap")
            print("\t--overlap_gap", self.args.overlap_gap)
            print("\t--overlap_k", self.args.overlap_k)
            print("\t--overlap_w", self.args.overlap_w)
            print("\t--btllib_t", self.args.btllib_t)

    def main_scaffolder(self):
        "Run ntJoin scaffolding stage"

        self.load_minimizers_scaffold()

        # Generate minimizer graph, and get paths through the graph
        self.make_minimizer_graph()

        self.graph = self.filter_graph_global(self.graph)

        self.mx_extremes = self.find_mx_min_max(self.args.s)

        # Load target scaffolds into memory
        min_match = re.search(r'^(\S+).k\d+.w\d+\.tsv', self.args.s)
        if not min_match:
            print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
            print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
            sys.exit(1)
        assembly_fa = min_match.group(1)
        self.scaffolds = self.read_fasta_file(assembly_fa)  # scaffold_id -> Scaffold

        # Find the paths through the graph
        paths = self.find_paths()

        # Format the paths to PathNodes, tally incorporated segments
        paths, incorporated_segments = self.format_adjust_paths(paths)

        if self.args.no_cut:
            paths = self.adjust_paths(paths, self.scaffolds, incorporated_segments)

        # Tally any regions that overlap
        intersecting_regions = self.tally_intersecting_segments(incorporated_segments)

        # Print the final scaffolds
        self.print_scaffolds(paths, intersecting_regions, incorporated_segments)

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)

    def set_weights(self):
        "Parse the supplied weights"
        weights = [float(w) for w in re.split(r'\s+', self.args.r)]
        if len(weights) != len(self.args.FILES):
            print("ERROR: The length of supplied reference weights (-r) and "
                  "number of assembly minimizer TSV inputs must be equal.")
            print("Supplied lengths of arguments:")
            print("Weights (-r):", len(weights), "Minimizer TSV files:", len(self.args.FILES), sep=" ")
            sys.exit(1)
        return weights

    def load_minimizers_scaffold(self):
        "Load in minimizers for ntJoin scaffolding mode"
        # Load in minimizers for references
        self.load_minimizers()
        # Now, add minimizers for reference
        mxs_info, mxs = ntjoin_utils.read_minimizers(self.args.s)
        self.list_mx_info[self.args.s] = mxs_info
        self.list_mxs[self.args.s] = mxs
        self.weights[self.args.s] = self.args.l

    def __init__(self, args):
        "Create an ntJoin instance for scaffolding"
        super().__init__(args)
        self.weights_list = self.set_weights()
        self.print_parameters_scaffold()
        self.mx_extremes = {} # ctg -> (min_pos, max_pos)
        self.scaffolds = {} # scaffold_id -> Scaffold
