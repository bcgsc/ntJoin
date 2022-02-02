#!/usr/bin/env python3
"""
ntJoin: Scaffolding assemblies using reference assemblies and minimizer graphs
Written by Lauren Coombe (@lcoombe)
"""

import argparse
import datetime
import multiprocessing
import re
from collections import Counter
from collections import defaultdict
import shlex
import subprocess
import sys
import warnings
import igraph as ig
import pybedtools
import pymannkendall as mk
import btllib
from read_fasta import read_fasta
import ntjoin_utils
import ntjoin_overlap
warnings.simplefilter(action='ignore', category=RuntimeWarning)


class Ntjoin:
    "ntJoin: Scaffolding assemblies using reference assemblies and minimizer graphs"

    # Helper functions for interfacing with python-igraph
    @staticmethod
    def vertex_index(graph, name):
        "Returns vertex index based on vertex name"
        return graph.vs.find(name).index

    @staticmethod
    def vertex_name(graph, index):
        "Returns vertex name based on vertex id"
        return graph.vs[index]['name']

    @staticmethod
    def edge_index(graph, source_name, target_name):
        "Returns graph edge index based on source/target names"
        return graph.get_eid(source_name, target_name)

    @staticmethod
    def set_edge_attributes(graph, edge_attributes):
        "Sets the edge attributes for a python-igraph graph"
        graph.es()["support"] = [edge_attributes[e]['support'] for e in sorted(edge_attributes.keys())]
        graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]

    @staticmethod
    def convert_path_index_to_name(graph, path):
        "Convert path of vertex indices to path of vertex names"
        return [Ntjoin.vertex_name(graph, vs) for vs in path]

    @staticmethod
    def read_minimizers(tsv_filename):
        "Read the minimizers from a file, removing duplicate minimizers"
        print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
        mx_info = {}  # mx -> (contig, position)
        mxs = []  # List of lists of minimizers
        dup_mxs = set()  # Set of minimizers identified as duplicates
        with open(tsv_filename, 'r') as tsv:
            for line in tsv:
                line = line.strip().split("\t")
                if len(line) > 1:
                    mx_pos_split = line[1].split(" ")
                    mxs.append([mx_pos.split(":")[0] for mx_pos in mx_pos_split])
                    for mx_pos in mx_pos_split:
                        mx, pos = mx_pos.split(":")
                        if mx in mx_info:  # This is a duplicate, add to dup set, don't add to dict
                            dup_mxs.add(mx)
                        else:
                            mx_info[mx] = (line[0], int(pos))

        mx_info = {mx: mx_info[mx] for mx in mx_info if mx not in dup_mxs}
        mxs_filt = []
        for mx_list in mxs:
            mx_list_filt = [mx for mx in mx_list if mx not in dup_mxs]
            mxs_filt.append(mx_list_filt)
        return mx_info, mxs_filt

    @staticmethod
    def calc_total_weight(list_files, weights):
        "Calculate the total weight of an edge given the assembly support"
        return sum([weights[f] for f in list_files])


    def build_graph(self, list_mxs, weights):
        "Builds an undirected graph: nodes=minimizers; edges=between adjacent minimizers"
        print(datetime.datetime.today(), ": Building graph", file=sys.stdout)
        graph = ig.Graph()

        vertices = set()
        edges = defaultdict(dict)  # source -> target -> [list assembly support]

        for assembly in list_mxs:
            for assembly_mx_list in list_mxs[assembly]:
                for i, j in zip(range(0, len(assembly_mx_list)),
                                range(1, len(assembly_mx_list))):
                    if assembly_mx_list[i] in edges and \
                            assembly_mx_list[j] in edges[assembly_mx_list[i]]:
                        edges[assembly_mx_list[i]][assembly_mx_list[j]].append(assembly)
                    elif assembly_mx_list[j] in edges and \
                            assembly_mx_list[i] in edges[assembly_mx_list[j]]:
                        edges[assembly_mx_list[j]][assembly_mx_list[i]].append(assembly)
                    else:
                        edges[assembly_mx_list[i]][assembly_mx_list[j]] = [assembly]
                    vertices.add(assembly_mx_list[i])
                if assembly_mx_list:
                    vertices.add(assembly_mx_list[-1])

        formatted_edges = [(s, t) for s in edges for t in edges[s]]

        print(datetime.datetime.today(), ": Adding vertices", file=sys.stdout)
        graph.add_vertices(list(vertices))

        print(datetime.datetime.today(), ": Adding edges", file=sys.stdout)
        graph.add_edges(formatted_edges)

        print(datetime.datetime.today(), ": Adding attributes", file=sys.stdout)
        edge_attributes = {self.edge_index(graph, s, t): {"support": edges[s][t],
                                                          "weight": self.calc_total_weight(edges[s][t],
                                                                                           weights)}
                           for s in edges for t in edges[s]}
        self.set_edge_attributes(graph, edge_attributes)

        return graph


    def print_graph(self, graph):
        "Prints the minimizer graph in dot format"
        out_graph = self.args.p + ".mx.dot"
        outfile = open(out_graph, 'w')
        print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

        outfile.write("graph G {\n")

        colours = ["red", "green", "blue", "purple", "orange",
                   "turquoise", "pink", "yellow", "orchid", "salmon"]
        list_files = list(Ntjoin.list_mx_info.keys())
        if len(list_files) > len(colours):
            colours = ["red"]*len(list_files)

        for node in graph.vs():
            mx_ctg_pos_labels = "\n".join([str(Ntjoin.list_mx_info[assembly][node['name']])
                                           for assembly in Ntjoin.list_mx_info])
            node_label = "\"%s\" [label=\"%s\n%s\"]" % (node['name'], node['name'], mx_ctg_pos_labels)
            outfile.write("%s\n" % node_label)

        for edge in graph.es():
            outfile.write("\"%s\" -- \"%s\"" %
                          (self.vertex_name(graph, edge.source),
                           self.vertex_name(graph, edge.target)))
            weight = edge['weight']
            support = edge['support']
            if len(support) == 1:
                colour = colours[list_files.index(support[0])]
            elif len(support) == 2:
                colour = "lightgrey"
            else:
                colour = "black"
            outfile.write(" [weight=%s color=%s]\n" % (weight, colour))

        outfile.write("}\n")

        print("\nfile_name\tnumber\tcolour")
        for i, filename in enumerate(list_files):
            print(filename, i, colours[i], sep="\t")
        print("")

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
            return 0

        # Find the assemblies that have a path between these mx
        # Are situations where there is not a direct edge if an unoriented contig was in-between
        path = graph.get_shortest_paths(u_mx, v_mx, output="vpath")[0]
        supporting_assemblies = set.intersection(
            *map(set, [graph.es()[self.edge_index(graph, s, t)]['support']
                       for s, t in zip(path, path[1:])]))
        if not supporting_assemblies:
            return self.args.g

        distances = [abs(Ntjoin.list_mx_info[assembly][v_mx][1] - Ntjoin.list_mx_info[assembly][u_mx][1])
                     for assembly in supporting_assemblies]
        mean_dist = int(sum(distances)/len(distances)) - self.args.k
        # Correct for the overhanging sequence before/after terminal minimizers
        if u.ori == "+":
            a = u.end - Ntjoin.list_mx_info[cur_assembly][u_mx][1] - self.args.k
        else:
            a = Ntjoin.list_mx_info[cur_assembly][u_mx][1] - u.start
        if v.ori == "+":
            b = Ntjoin.list_mx_info[cur_assembly][v_mx][1] - v.start
        else:
            b = v.end - Ntjoin.list_mx_info[cur_assembly][v_mx][1] - self.args.k

        try:
            assert a >= 0
            assert b >= 0
        except AssertionError as assert_error:
            print("ERROR: Gap distance estimation less than 0", "Vertex 1:", u, "Vertex 2:", v,
                  sep="\n")
            print("Minimizer positions:", Ntjoin.list_mx_info[cur_assembly][u_mx][1],
                  Ntjoin.list_mx_info[cur_assembly][v_mx][1])
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


    def merge_relocations(self, path):
        "If a path has adjacent collinear intervals of the same contig, merge them"
        if len(path) < 2:
            return path
        return_path = [path[0]]
        for node_i, node_j in zip(path, path[1:]):
            if node_i.contig == node_j.contig:
                if node_i.ori == "+" and node_j.ori == "+" and node_i.end <= node_j.start:
                    if self.is_new_region_overlapping(node_i.start, node_j.end, node_i, node_j,
                                                      Ntjoin.incorporated_segments[node_i.contig]):
                        return_path.append(node_j)
                        continue
                    Ntjoin.incorporated_segments[node_i.contig].add(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                        start=return_path[-1].start,
                                                                        end=node_j.end))
                    Ntjoin.incorporated_segments[node_i.contig].remove(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                           start=return_path[-1].start,
                                                                           end=return_path[-1].end))
                    Ntjoin.incorporated_segments[node_j.contig].remove(ntjoin_utils.Bed(contig=node_j.contig,
                                                                           start=node_j.start,
                                                                           end=node_j.end))
                    return_path[-1].end = node_j.end
                    return_path[-1].terminal_mx = node_j.terminal_mx
                    return_path[-1].gap_size = node_j.gap_size
                elif node_i.ori == "-" and node_j.ori == "-" and node_i.start >= node_j.end:
                    if self.is_new_region_overlapping(node_j.start, node_i.end, node_i, node_j,
                                                      Ntjoin.incorporated_segments[node_i.contig]):
                        return_path.append(node_j)
                        continue
                    Ntjoin.incorporated_segments[node_i.contig].add(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                        start=node_j.start,
                                                                        end=return_path[-1].end))
                    Ntjoin.incorporated_segments[node_i.contig].remove(ntjoin_utils.Bed(contig=return_path[-1].contig,
                                                                           start=return_path[-1].start,
                                                                           end=return_path[-1].end))
                    Ntjoin.incorporated_segments[node_j.contig].remove(ntjoin_utils.Bed(contig=node_j.contig,
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
            ctg, pos = Ntjoin.list_mx_info[assembly][mx]
            if ctg is curr_ctg:
                positions.append(pos)
            else:
                # This is either the first mx, or we are past a stretch of repeating contigs
                if curr_ctg is not None:
                    ori = self.determine_orientation(positions)
                    if ori != "?":  # Don't add to path if orientation couldn't be determined
                        out_path.append(ntjoin_utils.PathNode(contig=curr_ctg, ori=ori,
                                                 start=self.calc_start_coord(positions,
                                                                             Ntjoin.mx_extremes[curr_ctg][0]),
                                                 end=self.calc_end_coord(positions,
                                                                         Ntjoin.mx_extremes[curr_ctg][1],
                                                                         Ntjoin.scaffolds[curr_ctg].length),
                                                 contig_size=Ntjoin.scaffolds[curr_ctg].length,
                                                 first_mx=first_mx,
                                                 terminal_mx=prev_mx))
                curr_ctg = ctg
                positions = [pos]
                first_mx = mx
            prev_mx = mx
        ori = self.determine_orientation(positions)
        if ori != "?":
            out_path.append(ntjoin_utils.PathNode(contig=curr_ctg, ori=ori,
                                     start=self.calc_start_coord(positions,
                                                                 Ntjoin.mx_extremes[curr_ctg][0]),
                                     end=self.calc_end_coord(positions,
                                                             Ntjoin.mx_extremes[curr_ctg][1],
                                                             Ntjoin.scaffolds[curr_ctg].length),
                                     contig_size=Ntjoin.scaffolds[curr_ctg].length,
                                     first_mx=first_mx,
                                     terminal_mx=prev_mx))
        for u, v in zip(out_path, out_path[1:]):
            gap_size, raw_gap_size = self.calculate_gap_size(u, v, component_graph, assembly)
            u.set_gap_size(gap_size)
            u.set_raw_gap_size(raw_gap_size)

        return out_path

    @staticmethod
    def filter_graph(graph, min_weight):
        "Filter the graph by edge weights on edges incident to branch nodes"
        branch_nodes = [node.index for node in graph.vs() if node.degree() > 2]
        to_remove_edges = [edge for node in branch_nodes for edge in graph.incident(node)
                           if graph.es()[edge]['weight'] < min_weight]
        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    def filter_graph_global(self, graph):
        "Filter the graph globally based on minimum edge weight"
        print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout)
        if self.args.n <= min(Ntjoin.weights.values()):
            return graph
        to_remove_edges = [edge.index for edge in graph.es()
                           if edge['weight'] < self.args.n]
        new_graph = graph.copy()
        new_graph.delete_edges(to_remove_edges)
        return new_graph

    def determine_source_vertex(self, sources, graph):
        '''Given the possible sources of the graph, determine which is the source and the target
            Based on the assembly with the largest weight - orient others based on this assembly
        '''
        max_wt_asm = [assembly for assembly in Ntjoin.weights
                      if Ntjoin.weights[assembly] == max(Ntjoin.weights.values())].pop()
        list_mx_info_maxwt = Ntjoin.list_mx_info[max_wt_asm]
        min_pos = min([list_mx_info_maxwt[self.vertex_name(graph, s)][1] for s in sources])
        max_pos = max([list_mx_info_maxwt[self.vertex_name(graph, s)][1] for s in sources])
        source = [s for s in sources
                  if list_mx_info_maxwt[self.vertex_name(graph, s)][1] == min_pos].pop()
        target = [s for s in sources
                  if list_mx_info_maxwt[self.vertex_name(graph, s)][1] == max_pos].pop()
        return source, target

    @staticmethod
    def is_graph_linear(graph):
        "Given a graph, return True if all the components are linear"
        for component in graph.components():
            component_graph = graph.subgraph(component)
            if not all(u.degree() < 3 for u in component_graph.vs()):
                return False
        return True


    def find_paths_process(self, component):
        "Find paths given a component of the graph"
        return_paths = []
        min_edge_weight = self.args.n
        max_edge_weight = sum(Ntjoin.weights.values())
        component_graph = Ntjoin.gin.subgraph(component)
        while not self.is_graph_linear(component_graph) and \
                min_edge_weight <= max_edge_weight:
            component_graph = self.filter_graph(component_graph, min_edge_weight)
            min_edge_weight += 1

        for subcomponent in component_graph.components():
            subcomponent_graph = component_graph.subgraph(subcomponent)
            source_nodes = [node.index for node in subcomponent_graph.vs() if node.degree() == 1]
            if len(source_nodes) == 2:
                source, target = self.determine_source_vertex(source_nodes, subcomponent_graph)
                path = subcomponent_graph.get_shortest_paths(source, target)[0]
                num_edges = len(path) - 1
                if len(path) == len(subcomponent_graph.vs()) and \
                        num_edges == len(subcomponent_graph.es()) and len(path) == len(set(path)):
                    # All the nodes/edges from the graph are in the simple path, no repeated nodes
                    path = self.convert_path_index_to_name(subcomponent_graph, path)
                    ctg_path = self.format_path(path, self.args.s,
                                                subcomponent_graph)
                    return_paths.append(ctg_path)
        return return_paths

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

    def find_paths(self, graph):
        "Finds paths per input assembly file"
        print(datetime.datetime.today(), ": Finding paths", file=sys.stdout)
        Ntjoin.gin = graph
        components = graph.components()
        print("\nTotal number of components in graph:", len(components), "\n", sep=" ", file=sys.stdout)

        if self.args.t == 1:
            paths = [self.find_paths_process(component) for component in components]
        else:
            with multiprocessing.Pool(self.args.t) as pool:
                paths = pool.map(self.find_paths_process, components)

        paths_return = []
        incorporated_segments = {}
        for path_list in paths:
            for path in path_list:
                paths_return.append(path)
                self.tally_incorporated_segments(incorporated_segments, path)

        Ntjoin.incorporated_segments = incorporated_segments

        paths_return_merged = []
        for path in paths_return:
            path = self.merge_relocations(path)
            paths_return_merged.append(path)

        return paths_return_merged, incorporated_segments


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

    def adjust_paths(self, paths, scaffolds):
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
            new_path = self.merge_relocations(new_path)
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


    @staticmethod
    def read_fasta_file(filename):
        "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
        print(datetime.datetime.today(), ": Reading fasta file", filename, file=sys.stdout)
        scaffolds = {}
        try:
            with open(filename, 'r') as fasta:
                for header, seq, _, _ in read_fasta(fasta):
                    scaffolds[header] = ntjoin_utils.Scaffold(id=header, length=len(seq), sequence=seq)
        except FileNotFoundError:
            print("ERROR: File", filename, "not found.")
            print("Minimizer TSV file must follow the naming convention:")
            print("\tassembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering,\n"
                  "and assembly.fa is the scaffolds fasta file")
            sys.exit(1)
        return scaffolds

    @staticmethod
    def reverse_complement(sequence):
        "Reverse complements a given sequence"
        translation_table = str.maketrans(
            "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
            "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
        return sequence[::-1].translate(translation_table)

    @staticmethod
    def get_fasta_segment(path_node, sequence):
        "Given a PathNode and the contig sequence, return the corresponding sequence"
        if path_node.ori == "-":
            return Ntjoin.reverse_complement(sequence[path_node.start:path_node.end]) + \
                   "N"*path_node.gap_size
        return sequence[path_node.start:path_node.end] + "N"*path_node.gap_size

    @staticmethod
    def format_bedtools_genome(scaffolds):
        "Format a BED file and genome dictionary for bedtools"
        bed_str = "\n".join(["%s\t%d\t%d" % (scaffold, 0, scaffolds[scaffold].length)
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
                        path[i].start += len_diff
                    else:
                        path[i].end -= len_diff
                    assert len(sequence_start_strip) - path[i].gap_size == path[i].end - path[i].start
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

    def adjust_for_trimming(self, fasta_filename, path):
        "Go through path, trim the segments if overlapping"
        ct = 0
        mx_info = defaultdict(dict)  # path_index -> mx -> pos
        mxs = {}  # path_index -> [mx]
        with btllib.Indexlr(fasta_filename, self.args.overlap_k, self.args.overlap_w,
                            btllib.IndexlrFlag.LONG_MODE, self.args.overlap_t) as minimizers:
            for mx_entry in minimizers:
                mxs[ct] = []
                dup_mxs = set()  # Set of minimizers identified as duplicates
                for mx_pos_strand in mx_entry.minimizers:
                    mx, pos = str(mx_pos_strand.out_hash), mx_pos_strand.pos
                    if not ntjoin_overlap.is_in_valid_region(pos, ct, path):
                        continue
                    if ct in mx_info and mx in mx_info[ct]:  # This is a duplicate
                        dup_mxs.add(mx)
                    else:
                        mx_info[ct][mx] = int(pos)
                        mxs[ct].append(mx)
                mx_info[ct] = {mx: mx_info[ct][mx] for mx in mx_info[ct] if mx not in dup_mxs}
                mxs[ct] = [[mx for mx in mxs[ct] if mx not in dup_mxs and mx in mx_info[ct]]]
                ct += 1

        for i in range(0, len(path) - 1):
            source = path[i]
            if source.raw_gap_size < 0:
                ntjoin_overlap.merge_overlapping(mxs, mx_info, i, i+1, path)



    @staticmethod
    def update_graph_tally(path, vertices, edges):
        "Update graph vertices/edges with given path"
        for s, t in zip(path, path[1:]):
            formatted_s, formatted_t = "{}{}".format(s.contig, s.ori), "{}{}".format(t.contig, t.ori)
            vertices.add(formatted_s)
            vertices.add(formatted_t)
            edges.add(ntjoin_utils.EdgeGraph(source=formatted_s,
                                             target=formatted_t,
                                             raw_gap_est=s.raw_gap_size))

    def print_scaffold_graph(self, vertices, edges):
        "Print scaffold graph generated from path file"
        outfile = open(self.args.p + ".scaffold.dot", 'w')
        outfile.write("digraph G {\n")

        for node in vertices:
            node_label = "\"{scaffold}\" [l={length}]\n". \
                format(scaffold=node, length=Ntjoin.scaffolds[node.strip("+-")].length)
            outfile.write(node_label)

        for edge in edges:
            edge_str = "\"{source}\" -> \"{target}\" [d={d} e={e} n={n}]\n". \
                format(source=edge.source, target=edge.target,
                       d=edge.raw_gap_est, e=100, n=1)
            outfile.write(edge_str)

        outfile.write("}\n")
        outfile.close()

    def get_adjusted_sequence(self, sequence, node):
        "Return sequence adjusted for overlap trimming"
        return_sequence = sequence[node.start_adjust:node.get_end_adjusted_coordinate()]
        if node.gap_size > 0:
            return return_sequence + self.args.overlap_gap*"N"
        return return_sequence


    def print_scaffolds(self, paths, intersecting_regions):
        "Given the paths, print out the scaffolds fasta"
        print(datetime.datetime.today(), ": Printing output scaffolds", file=sys.stdout)
        assembly = self.args.s

        min_match = re.search(r'^(\S+)(.k\d+.w\d+)\.tsv', assembly)
        assembly_fa, params = min_match.group(1), min_match.group(2)

        outfile = open(assembly_fa + params + ".n" +
                       str(self.args.n) + ".assigned.scaffolds.fa", 'w')
        pathfile = open(self.args.p + ".path", 'w')
        if self.args.agp:
            agpfile = open(self.args.p + ".agp", "w")
        if self.args.overlap:
            vertices = set()
            edges = set() # Set of EdgeGraph
        incorporated_segments = []  # List of Bed entries

        ct = 0
        pathfile.write(assembly_fa + "\n")
        for path in paths:
            sequences = []
            nodes = []
            path_segments = []

            path = self.merge_relocations(path)

            path = self.remove_overlapping_regions(path, intersecting_regions)

            self.check_terminal_node_gap_zero(path)

            for node in path:
                if node.ori == "?":
                    continue
                sequences.append(self.get_fasta_segment(node, Ntjoin.scaffolds[node.contig].sequence))
                path_segments.append(ntjoin_utils.Bed(contig=node.contig, start=node.start,
                                                      end=node.end))
                nodes.append(node)
            if len(sequences) < 2:
                continue

            if self.args.overlap:
                path_segments_file = open(self.args.p + ".segments.fa", 'w')
                for seq, node in zip(sequences, nodes):  # !! TODO: limit to overlapping section?
                    path_segments_file.write(">{}_{}-{} {}\n{}\n".format(node.contig, node.start,
                                                                         node.end, node.raw_gap_size, seq.strip()))
                path_segments_file.close()

            if self.args.overlap:
                self.adjust_for_trimming(self.args.p + ".segments.fa", nodes)
                sequences = [self.get_adjusted_sequence(sequence, nodes[i])
                             for i, sequence in enumerate(sequences)]
                cmd = shlex.split("rm {}".format(self.args.p + ".segments.fa"))
                return_code = subprocess.call(cmd)
                assert return_code == 0

            ## For debugging
            # !! TODO remove
            for seq, path in zip(sequences, path_segments):
                print(">{}_{}-{}\n{}".format(path.contig, path.start, path.end, seq), file=sys.stderr)

            ctg_id = "ntJoin" + str(ct)
            ctg_sequence = self.join_sequences(sequences, path, path_segments)

            outfile.write(">%s\n%s\n" %
                          (ctg_id, ctg_sequence))
            incorporated_segments.extend(path_segments)
            path_str = " ".join(["%s%s:%d-%d %dN" %
                                 (node.contig, node.ori, node.get_adjusted_start(), node.get_adjusted_end(),
                                  node.gap_size) for node in path])
            path_str = re.sub(r'\s+\d+N$', r'', path_str)
            pathfile.write("%s\t%s\n" % (ctg_id, path_str))
            if self.args.agp:
                self.write_agp(agpfile, ctg_id, path_str)
            if self.args.overlap:
                self.update_graph_tally(path, vertices, edges)
            ct += 1
        outfile.close()

        if self.args.agp:
            outfile = self.print_unassigned(assembly, assembly_fa, incorporated_segments, outfile, params, agpfile=agpfile)
        else:
            outfile = self.print_unassigned(assembly, assembly_fa, incorporated_segments, outfile, params)

        if self.args.overlap:
            self.print_scaffold_graph(vertices, edges)

        outfile.close()
        pathfile.close()
        if self.args.agp:
            agpfile.close()

    def print_unassigned(self, assembly, assembly_fa, incorporated_segments, outfile, params, agpfile=None):
        "Also print out the sequences that were NOT scaffolded"
        incorporated_segments_str = "\n".join(["%s\t%s\t%s" % (chrom, s, e)
                                               for chrom, s, e in incorporated_segments])
        incorporated_segments_bed = pybedtools.BedTool(incorporated_segments_str,
                                                       from_string=True).sort()
        genome_bed, genome_dict = self.format_bedtools_genome(Ntjoin.scaffolds)
        missing_bed = genome_bed.complement(i=incorporated_segments_bed, g=genome_dict)
        missing_bed.saveas(self.args.p + "." + assembly + ".unassigned.bed")
        outfile = open(assembly_fa + params + ".n" +
                       str(self.args.n) + ".unassigned.scaffolds.fa", 'w')
        cmd = "bedtools getfasta -fi %s -bed %s -fo -" % \
              (assembly_fa, self.args.p + "." + assembly + ".unassigned.bed")
        cmd_shlex = shlex.split(cmd)
        out_fasta = subprocess.Popen(cmd_shlex, stdout=subprocess.PIPE, universal_newlines=True)
        for header, seq, _, _ in read_fasta(iter(out_fasta.stdout.readline, '')):
            if self.args.agp:
                self.write_agp_unassigned(agpfile, header, seq)
            seq = seq.strip().strip("Nn")
            if seq:
                outfile.write(">{header}\n{seq}\n".format(header=header, seq=seq))
        out_fasta.wait()
        if out_fasta.returncode != 0:
            print("bedtools getfasta failed -- is bedtools on your PATH?")
            print(out_fasta.stderr)
            raise subprocess.CalledProcessError(out_fasta.returncode, cmd_shlex)
        return outfile

    @staticmethod
    def tally_intersecting_segments():
        "Tally ctgs with intersecting segments, and keep track of 'best'"
        incorporated_bed_list = []
        for ctg in Ntjoin.incorporated_segments:
            for bed_entry in Ntjoin.incorporated_segments[ctg]:
                incorporated_bed_list.append(bed_entry)
        incorporated_bed_str = "\n".join(["%s\t%s\t%s" % (chrom, s, e)
                                          for chrom, s, e in incorporated_bed_list])
        incorporated_segments_bed = pybedtools.BedTool(incorporated_bed_str,
                                                       from_string=True).sort()
        bed_intersect = incorporated_segments_bed.intersect(b=incorporated_segments_bed,
                                                            c=True, wa=True)

        overlap_regions = {}

        for bed in bed_intersect:
            if bed.count > 1:
                if bed.chrom not in overlap_regions:
                    overlap_regions[bed.chrom] = ntjoin_utils.OverlapRegion()
                overlap_regions[bed.chrom].add_region(ntjoin_utils.Bed(contig=bed.chrom, start=bed.start, end=bed.end))

        overlap_regions_fix = {}
        for overlap_contig in overlap_regions:
            overlap_regions_fix[overlap_contig] = overlap_regions[overlap_contig].find_non_overlapping()

        return overlap_regions_fix

    @staticmethod
    def find_mx_min_max(graph, target):
        "Given a dictionary in the form mx->(ctg, pos), find the min/max mx position per ctg"
        mx_extremes = {} # ctg -> (min_pos, max_pos)
        for mx in Ntjoin.list_mx_info[target]:
            try:
                graph.vs().find(mx)
            except ValueError:
                continue
            ctg, pos = Ntjoin.list_mx_info[target][mx]
            if ctg in mx_extremes:
                mx_extremes[ctg] = (min(mx_extremes[ctg][0], pos),
                                    max(mx_extremes[ctg][1], pos))
            else:
                mx_extremes[ctg] = (pos, pos)
        return mx_extremes


    @staticmethod
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
        parser.add_argument("-v", "--version", action='version', version='ntJoin v1.0.8')
        parser.add_argument("--agp", help="Output AGP file describing scaffolds", action="store_true")
        parser.add_argument("--no_cut", help="Do not cut input contigs, place in most representative path",
                            action="store_true")
        parser.add_argument("--overlap", help="Print scaffold graph form of paths, including putative overlaps",
                            action="store_true")
        parser.add_argument("--overlap_gap", help="Length of gap introduced between overlapping, trimmed segments",
                            type=int, default=20)
        parser.add_argument("--overlap_k", help="Kmer size used for overlap minimizer step",
                            type=int, default=15)
        parser.add_argument("--overlap_w", help="Window size used for overlap minimizer step",
                            type=int, default=10)
        parser.add_argument("--overlap_t", help="Number of threads for computing overlap minimizers",
                            type=int, default=4)
        return parser.parse_args()

    def print_parameters(self):
        "Print the set parameters for the ntJoin run"
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
            print("\t--overlap_t", self.args.overlap_t)

    def main(self):
        "Run ntJoin graph stage"
        print("Running ntJoin v1.0.8 ...\n")
        self.print_parameters()

        # Parse the weights of each input reference assembly
        input_weights = [float(w) for w in re.split(r'\s+', self.args.r)]
        if len(input_weights) != len(self.args.FILES):
            print("ERROR: The length of supplied reference weights (-r) and "
                  "number of assembly minimizer TSV inputs must be equal.")
            print("Supplied lengths of arguments:")
            print("Weights (-r):", len(input_weights), "Minimizer TSV files:", len(self.args.FILES), sep=" ")
            sys.exit(1)

        # Read in the minimizers for each assembly
        list_mx_info = {}  # Dictionary of dictionaries: assembly -> mx -> (contig, position)
        list_mxs = {}  # Dictionary: assembly -> [lists of mx]
        weights = {}  # Dictionary: assembly -> weight
        for assembly in self.args.FILES:
            mxs_info, mxs = self.read_minimizers(assembly)
            list_mx_info[assembly] = mxs_info
            list_mxs[assembly] = mxs
            weights[assembly] = input_weights.pop(0)
        mxs_info, mxs = self.read_minimizers(self.args.s)
        list_mx_info[self.args.s] = mxs_info
        list_mxs[self.args.s] = mxs
        weights[self.args.s] = self.args.l
        weight_str = "\n".join(["%s: %s" % (assembly, weights[assembly]) for assembly in weights])
        print("\nWeights of assemblies:\n", weight_str, "\n", sep="")

        Ntjoin.list_mx_info = list_mx_info
        Ntjoin.weights = weights

        # Filter minimizers - Keep only if found in all assemblies
        list_mxs = ntjoin_utils.filter_minimizers(list_mxs)

        # Build a graph: Nodes = mx; Edges between adjacent mx in the assemblies
        graph = self.build_graph(list_mxs, Ntjoin.weights)

        # Print the DOT graph
        self.print_graph(graph)

        # Filter the graph edges by minimum weight
        graph = self.filter_graph_global(graph)

        # Find the min and max pos of minimizers for target assembly, per ctg
        Ntjoin.mx_extremes = self.find_mx_min_max(graph, self.args.s)

        # Load target scaffolds into memory
        min_match = re.search(r'^(\S+).k\d+.w\d+\.tsv', self.args.s)
        if not min_match:
            print("ERROR: Target assembly minimizer TSV file must follow the naming convention:")
            print("\ttarget_assembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
            sys.exit(1)
        assembly_fa = min_match.group(1)
        scaffolds = self.read_fasta_file(assembly_fa)  # scaffold_id -> Scaffold

        Ntjoin.scaffolds = scaffolds

        # Find the paths through the graph
        paths, incorporated_segments = self.find_paths(graph)

        Ntjoin.incorporated_segments = incorporated_segments

        if self.args.no_cut:
            paths = self.adjust_paths(paths, scaffolds)

        # Tally any regions that overlap
        intersecting_regions = self.tally_intersecting_segments()

        # Print the final scaffolds
        self.print_scaffolds(paths, intersecting_regions)

        print(datetime.datetime.today(), ": DONE!", file=sys.stdout)

    def __init__(self):
        "Create an ntJoin instance"
        self.args = self.parse_arguments()

def main():
    "Run ntJoin"
    Ntjoin().main()

if __name__ == "__main__":
    main()
