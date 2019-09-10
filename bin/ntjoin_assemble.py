#!/usr/bin/env python3
"""
ntJoin: Scaffold multiple genome assemblies using minimizers
Written by Lauren Coombe (@lcoombe)
"""

import argparse
import datetime
import multiprocessing
import re
from collections import Counter
from collections import defaultdict
from collections import namedtuple
import shlex
import subprocess
import sys
import igraph as ig
import pybedtools
import pymannkendall as mk
from read_fasta import read_fasta


# Defining namedtuples
Bed = namedtuple("Bed", ["contig", "start", "end"])
Scaffold = namedtuple("Scaffold", ["id", "length", "sequence"])

# Defining helper classes
class PathNode:
    "Defines a node in a path"
    def __init__(self, contig, ori, start, end, contig_size,
                 first_mx, terminal_mx, gap_size=50):
        self.contig = contig
        self.ori = ori
        self.start = start
        self.end = end
        self.contig_size = contig_size
        self.first_mx = first_mx
        self.terminal_mx = terminal_mx
        self.gap_size = gap_size

    def set_gap_size(self, gap_size):
        "Set the gap size of the path node"
        self.gap_size = gap_size

    def get_aligned_length(self):
        "Get the aligned length based on start/end coordinates"
        return self.end - self.start

    def __str__(self):
        return "Contig:%s\tOrientation:%s\tStart-End:%d-%d\tLength:%s\tFirstmx:%s\tLastmx:%s" \
               % (self.contig, self.ori, self.start, self.end, self.contig_size,
                  self.first_mx, self.terminal_mx)
class Ntjoin:
    "ntJoin: Scaffolding assemblies with minimizers"

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
        print("Reading minimizers:", tsv_filename, datetime.datetime.today(), file=sys.stdout)
        mx_info = {}  # mx -> (contig, position)
        mxs = []  # List of lists of minimizers (ordered)
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
    def filter_minimizers(list_mxs):
        "Filters out minimizers that are not found in all assemblies"
        print("Filtering minimizers", datetime.datetime.today(), file=sys.stdout)
        list_mx_sets = [{mx for mx_list in list_mxs[assembly] for mx in mx_list}
                        for assembly in list_mxs]

        mx_intersection = set.intersection(*list_mx_sets)

        return_mxs = {}
        for assembly in list_mxs:
            assembly_mxs_filtered = [[mx for mx in mx_list if mx in mx_intersection]
                                     for mx_list in list_mxs[assembly]]
            return_mxs[assembly] = assembly_mxs_filtered

        return return_mxs

    @staticmethod
    def calc_total_weight(list_files, weights):
        "Calculate the total weight of an edge given the assembly support"
        return sum([weights[f] for f in list_files])


    def build_graph(self, list_mxs):
        "Builds an undirected graph: nodes=minimizers; edges=between adjacent minimizers"
        print("Building graph", datetime.datetime.today(), file=sys.stdout)
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

        print("Adding vertices", datetime.datetime.today(), file=sys.stdout)
        graph.add_vertices(list(vertices))

        print("Adding edges", datetime.datetime.today(), file=sys.stdout)
        graph.add_edges(formatted_edges)

        print("Adding attributes", datetime.datetime.today(), file=sys.stdout)
        edge_attributes = {self.edge_index(graph, s, t): {"support": edges[s][t],
                                                          "weight": self.calc_total_weight(edges[s][t],
                                                                                           Ntjoin.weights)}
                           for s in edges for t in edges[s]}
        self.set_edge_attributes(graph, edge_attributes)

        return graph


    def print_graph(self, graph):
        "Prints a graph in dot format"
        out_graph = self.args.p + ".mx.dot"
        outfile = open(out_graph, 'w')
        print("Printing graph", out_graph, datetime.datetime.today(), sep=" ", file=sys.stdout)

        outfile.write("graph G {\n")

        # TODO: Make this more general
        colours = ["red", "green", "blue", "purple", "orange",
                   "turquoise", "pink", "yellow", "orchid", "salmon"]
        list_files = list(Ntjoin.list_mx_info.keys())
        if len(list_files) > len(colours):
            colours = ["red"]*len(list_files)

        for node in graph.vs():
            files_labels = "\n".join([str(Ntjoin.list_mx_info[assembly][node['name']])
                                      for assembly in Ntjoin.list_mx_info])
            node_label = "\"%s\" [label=\"%s\n%s\"]" % (node['name'], node['name'], files_labels)
            outfile.write("%s\n" % node_label)

        for edge in graph.es():
            outfile.write("\"%s\" -- \"%s\"" %
                          (self.vertex_name(graph, edge.source),
                           self.vertex_name(graph, edge.target)))
            # For debugging only
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

        print("file_name\tnumber\tcolour")
        for i, filename in enumerate(list_files):
            print(filename, i, colours[i], sep="\t")

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
        except:
            print("ERROR: Gap distance estimation less than 0", "Vertex 1:", u, "Vertex 2:", v,
                  sep="\n")
            print("Minimizer positions:", Ntjoin.list_mx_info[cur_assembly][u_mx][1],
                  Ntjoin.list_mx_info[cur_assembly][v_mx][1])
            print("Estimated distance: ", mean_dist)
            raise AssertionError

        gap_size = max(mean_dist - a - b, self.args.g)
        return gap_size


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
                        out_path.append(PathNode(contig=curr_ctg, ori=ori,
                                                 start=self.calc_start_coord(positions,
                                                                             Ntjoin.mx_extremes[assembly][curr_ctg][0]),
                                                 end=self.calc_end_coord(positions,
                                                                         Ntjoin.mx_extremes[assembly][curr_ctg][1],
                                                                         Ntjoin.scaffolds[assembly][curr_ctg].length),
                                                 contig_size=Ntjoin.scaffolds[assembly][curr_ctg].length,
                                                 first_mx=first_mx,
                                                 terminal_mx=prev_mx))
                curr_ctg = ctg
                positions = [pos]
                first_mx = mx
            prev_mx = mx
        ori = self.determine_orientation(positions)
        out_path.append(PathNode(contig=curr_ctg, ori=ori,
                                 start=self.calc_start_coord(positions,
                                                             Ntjoin.mx_extremes[assembly][curr_ctg][0]),
                                 end=self.calc_end_coord(positions,
                                                         Ntjoin.mx_extremes[assembly][curr_ctg][1],
                                                         Ntjoin.scaffolds[assembly][curr_ctg].length),
                                 contig_size=Ntjoin.scaffolds[assembly][curr_ctg].length,
                                 first_mx=first_mx,
                                 terminal_mx=prev_mx))
        for u, v in zip(out_path, out_path[1:]):
            gap_size = self.calculate_gap_size(u, v, component_graph, assembly)
            u.gap_size = gap_size

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
        "Filter the graph globally based on min edge weight"
        print("Filtering the graph", datetime.datetime.today(), file=sys.stdout)
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
            if not all(u.degree() < 2 for u in component_graph.vs()):
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
                    for assembly in Ntjoin.list_mx_info:
                        ctg_path = self.format_path(path, assembly,
                                                    subcomponent_graph)
                        return_paths.append((assembly, ctg_path))
        return return_paths


    def find_paths(self, graph):
        "Finds paths per input assembly file"
        print("Finding paths", datetime.datetime.today(), file=sys.stdout)
        paths = {}
        for assembly in Ntjoin.list_mx_info:
            paths[assembly] = []
        Ntjoin.gin = graph
        components = graph.components()
        print("Total number of components in graph:", len(components), sep=" ", file=sys.stdout)

        with multiprocessing.Pool(self.args.t) as pool:
            path_tuples = pool.map(self.find_paths_process, components)

        for my_list in path_tuples:
            for assembly, path in my_list:
                paths[assembly].append(path)

        return paths

    @staticmethod
    def read_fasta_file(filename):
        "Read a fasta file into memory. Returns dictionary of scafID -> Scaffold"
        print("Reading fasta file", filename, datetime.datetime.today(), file=sys.stdout)
        scaffolds = {}
        try:
            with open(filename, 'r') as fasta:
                for header, seq, _, _ in read_fasta(fasta):
                    scaffolds[header] = Scaffold(id=header, length=len(seq), sequence=seq)
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
            return Ntjoin.reverse_complement(sequence[path_node.start:path_node.end+1]) + \
                   "N"*path_node.gap_size
        return sequence[path_node.start:path_node.end+1] + "N"*path_node.gap_size

    @staticmethod
    def format_bedtools_genome(scaffolds):
        "Format a BED file and genome dictionary for bedtools"
        bed_str = "\n".join(["%s\t%d\t%d" % (scaffold, 0, scaffolds[scaffold].length)
                             for scaffold in scaffolds])
        bed = pybedtools.BedTool(bed_str, from_string=True)

        genome_dict = {scaffold: (0, scaffolds[scaffold].length) for scaffold in scaffolds}

        return bed, genome_dict


    def print_scaffolds(self, paths):
        "Given the paths, print out the scaffolds fasta"
        print("Printing output scaffolds", datetime.datetime.today(), file=sys.stdout)
        pathfile = open(self.args.p + ".path", 'w')

        for assembly in paths:
            min_match = re.search(r'^(\S+)(.k\d+.w\d+)\.tsv', assembly)
            assembly_fa = min_match.group(1)
            outfile = open(assembly_fa + min_match.group(2) + ".n" +
                           str(self.args.n) + ".assigned.scaffolds.fa", 'w')
            all_scaffolds = Ntjoin.scaffolds[assembly]
            incorporated_segments = []  # List of Bed entries

            ct = 0
            pathfile.write(assembly + "\n")
            for path in paths[assembly]:
                sequences = []
                path_segments = []
                for node in path:
                    if node.ori == "?":
                        continue
                    sequences.append(self.get_fasta_segment(node, all_scaffolds[node.contig].sequence))
                    path_segments.append(Bed(contig=node.contig, start=node.start,
                                             end=node.end))
                if len(sequences) < 2:
                    continue
                outfile.write(">%s\n%s\n" %
                              ("mx" + str(ct),
                               "".join(sequences).strip("Nn")))
                incorporated_segments.extend(path_segments)
                path_str = " ".join(["%s%s:%d-%d %dN" %
                                     (tup.contig, tup.ori, tup.start, tup.end, tup.gap_size)
                                     for tup in path])
                path_str = re.sub(r'\s+\d+N$', r'', path_str)
                pathfile.write("%s\t%s\n" % ("mx" + str(ct), path_str))
                ct += 1
            outfile.close()

            # Also print out the sequences that were NOT scaffolded
            incorporated_segments_str = "\n".join(["%s\t%s\t%s" % (chrom, s, e)
                                                   for chrom, s, e in incorporated_segments])
            incorporated_segments_bed = pybedtools.BedTool(incorporated_segments_str,
                                                           from_string=True).sort()
            genome_bed, genome_dict = self.format_bedtools_genome(all_scaffolds)

            missing_bed = genome_bed.complement(i=incorporated_segments_bed, g=genome_dict)
            missing_bed.saveas(self.args.p + "." + assembly + ".unassigned.bed")

            outfile = open(assembly_fa + min_match.group(2) + ".n" +
                           str(self.args.n) + ".unassigned.scaffolds.fa", 'w')

            cmd = "bedtools getfasta -fi %s -bed %s -fo -" % \
                  (assembly_fa, self.args.p + "." + assembly + ".unassigned.bed")
            cmd_shlex = shlex.split(cmd)

            out_fasta = subprocess.Popen(cmd_shlex, stdout=subprocess.PIPE, universal_newlines=True)
            for line in iter(out_fasta.stdout.readline, ''):
                outfile.write(line)
            out_fasta.wait()
            if out_fasta.returncode != 0:
                print("bedtools getfasta failed -- is bedtools on your PATH?")
                print(out_fasta.stderr)
                raise subprocess.CalledProcessError(out_fasta.returncode, cmd_shlex)

            outfile.close()
        pathfile.close()

    @staticmethod
    def find_mx_min_max(graph):
        "Given a dictionary in the form assembly->mx->(ctg, pos), find the min/max mx position per ctg"
        mx_extremes = {} # assembly -> ctg -> (min_pos, max_pos)
        for assembly in Ntjoin.list_mx_info:
            mx_extremes[assembly] = {}
            for mx in Ntjoin.list_mx_info[assembly]:
                try:
                    graph.vs().find(mx)
                except ValueError:
                    continue
                ctg, pos = Ntjoin.list_mx_info[assembly][mx]
                if ctg in mx_extremes[assembly]:
                    mx_extremes[assembly][ctg] = (min(mx_extremes[assembly][ctg][0], pos),
                                                  max(mx_extremes[assembly][ctg][1], pos))
                else:
                    mx_extremes[assembly][ctg] = (pos, pos)
        return mx_extremes


    @staticmethod
    def parse_arguments():
        "Parse ntJoin arguments"
        parser = argparse.ArgumentParser(
            description="ntJoin: Scaffold multiple genome assemblies using minimizers",
            epilog="Note: Script expects each minimizer TSV file has a matching fasta file.\n"
                   "Ex: myscaffolds.fa.k32.w1000.tsv - myscaffolds.fa is matching fasta",
            formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument("FILES", nargs="+", help="Minimizer TSV files")
        parser.add_argument("-p", help="Output prefix [out]", default="out",
                            type=str, required=False)
        parser.add_argument("-n", help="Minimum edge weight [2]", default=2, type=int)
        parser.add_argument("-l",
                            help="List of assembly weights (in quotes, separated by spaces, "
                                 "in same order as minimizer TSV files)",
                            required=True, type=str)
        parser.add_argument("-k", help="k value used for minimizer step", required=True, type=int)
        parser.add_argument("-g", help="Minimum gap size", required=False, default=1, type=int)
        parser.add_argument("--mkt", help="Use Mann-Kendall Test to orient contigs (slower)",
                            action='store_true')
        parser.add_argument('-m', help="Require at least m % of minimizer positions to be "
                                       "increasing/decreasing to assign contig orientation [90]\n "
                                       "Note: Only used with --mkt is NOT specified", default=90, type=int)
        parser.add_argument('-t', help="Number of threads [4]", default=4, type=int)
        parser.add_argument("-v", "--version", action='version', version='ntJoin v0.0.1')
        return parser.parse_args()

    def main(self):
        "Run ntJoin graph stage"
        # Sanity checking of user's specified arguments
        print("Running ntJoin...")
        if self.args.mkt:
            print("Orienting contigs with Mann-Kendall Test (more computationally intensive)")
        else:
            print("Orienting contigs using increasing/decreasing minimizer positions")

        # Parse the weights of each input assembly
        input_weights = re.split(r'\s+', self.args.l)
        if len(input_weights) != len(self.args.FILES):
            print("ERROR: The length of supplied weights and number of assembly minimizer must be equal.")
            print("Supplied lengths of arguments:")
            print("Weights (-l):", len(input_weights), "Minimizer TSV files:", len(self.args.FILES), sep=" ")
            sys.exit(1)

        # Read in the minimizers for each assembly
        list_mx_info = {}  # Dictionary of dictionaries: assembly -> mx -> (contig, position)
        list_mxs = {}  # Dictionary of lists of minimizers: assembly -> [lists of mx]
        weights = {}  # Dictionary: assembly -> weight
        for assembly in self.args.FILES:
            mxs_info, mxs = self.read_minimizers(assembly)
            list_mx_info[assembly] = mxs_info
            list_mxs[assembly] = mxs
            weights[assembly] = float(input_weights.pop(0))
        Ntjoin.list_mx_info = list_mx_info
        Ntjoin.weights = weights

        # Filter minimizers - Keep only if unique in an assembly and found in all assemblies
        list_mxs = self.filter_minimizers(list_mxs)

        # Build a graph: Nodes = mx; Edges between adjacent mx in the assemblies
        graph = self.build_graph(list_mxs)

        # Print the DOT graph
        self.print_graph(graph)

        # Filter the graph edges by minimum weight
        graph = self.filter_graph_global(graph)

        # Find the min and max pos of minimizers per assembly, per ctg
        Ntjoin.mx_extremes = self.find_mx_min_max(graph)

        # Load scaffolds into memory
        scaffolds = {} # assembly -> scaffold_id -> Scaffold
        for assembly in self.args.FILES:
            min_match = re.search(r'^(\S+).k\d+.w\d+\.tsv', assembly)
            if not min_match:
                print("ERROR: Minimizer TSV file must follow the naming convention:")
                print("\tassembly.fa.k<k>.w<w>.tsv, where <k> and <w> are parameters used for minimizering")
                sys.exit(1)
            assembly_fa = min_match.group(1)
            scaffolds[assembly] = self.read_fasta_file(assembly_fa)
        Ntjoin.scaffolds = scaffolds

        paths = self.find_paths(graph)

        self.print_scaffolds(paths)

    def __init__(self):
        "Create an ntJoin instance"
        self.args = self.parse_arguments()

def main():
    "Run ntJoin"
    Ntjoin().main()

if __name__ == "__main__":
    main()
