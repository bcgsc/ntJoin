#!/usr/bin/env python3
"""
Scaffold multiple genome assemblies using minimizers
Written by Lauren Coombe
"""

import argparse
import re
from collections import defaultdict
from collections import namedtuple
import networkx as nx
import pybedtools
import shlex
import subprocess
import sys
from read_fasta import read_fasta


# Defining namedtuples
PathNode = namedtuple("PathNode", ["contig", "ori", "start", "end"])
Bed = namedtuple("Bed", ["contig", "start", "end"])


def read_minimizers(tsv_filename):
    "Read all the minimizers from a file into a dictionary, removing duplicate minimizers"
    mx_info = {}  # mx -> (contig, position)
    mxs = []  # List of lists of minimizers, which have ordering information
    dup_mxs = set()  # Set of minimizers seen to be duplicates
    with open(tsv_filename, 'r') as tsv:
        for line in tsv:
            line = line.strip().split("\t")
            if len(line) > 1:
                mxs.append([mx_pos.split(":")[0] for mx_pos in line[1].split(" ")])
                for mx_pos in line[1].split(" "):
                    mx, pos = mx_pos.split(":")
                    if mx in mx_info:  # This is a duplicate, add to dup set, don't add to dict
                        dup_mxs.add(mx)
                    else:
                        mx_info[mx] = (line[0], int(pos))

    mx_info = {key: mx_info[key] for key in mx_info if key not in dup_mxs}
    mxs_filt = []
    for mx_list in mxs:
        mx_list_filt = [mx for mx in mx_list if mx not in dup_mxs]
        mxs_filt.append(mx_list_filt)
    return mx_info, mxs_filt


def filter_minimizers(list_mxs):
    "Filters out minimizers from each dictionary of lists that are not found in all other sets"
    list_sets = []
    for assembly in list_mxs:
        list_sets.append(set([mx for mx_list in list_mxs[assembly] for mx in mx_list]))

    mx_intersection = set.intersection(*list_sets)

    return_mxs = {}
    for assembly in list_mxs:
        assembly_mxs_filtered = [[mx for mx in mx_list if mx in mx_intersection] \
                                 for mx_list in list_mxs[assembly]]
        return_mxs[assembly] = assembly_mxs_filtered

    return return_mxs


def calc_total_weight(list_files, weights):
    "Given a list of supporting files for an edges, and dictionary specifying their weights, calculate the total weight"
    return sum([weights[f] for f in list_files])


def build_graph(list_mxs, weights):
    "Builds an undirected graph: minimizers=nodes; edges=between adjacent minimizers"
    graph = nx.Graph()

    edges = defaultdict(dict)  # source -> target -> [list assembly support]

    for assembly in list_mxs:
        for assembly_mx_list in list_mxs[assembly]:
            for i in range(0, len(assembly_mx_list)-1):
                if assembly_mx_list[i] in edges and \
                        assembly_mx_list[i+1] in edges[assembly_mx_list[i]]:
                    edges[assembly_mx_list[i]][assembly_mx_list[i+1]].append(assembly)
                elif assembly_mx_list[i+1] in edges and \
                        assembly_mx_list[i] in edges[assembly_mx_list[i+1]]:
                    edges[assembly_mx_list[i+1]][assembly_mx_list[i]].append(assembly)
                else:
                    edges[assembly_mx_list[i]][assembly_mx_list[i+1]] = [assembly]

    formatted_edges = [(s, t) for s in edges for t in edges[s]]
    edge_attributes = {(s, t): {"support": edges[s][t],
                                "weight": calc_total_weight(edges[s][t], weights)}
                       for s in edges for t in edges[s]}

    graph.add_edges_from(formatted_edges)
    nx.set_edge_attributes(graph, edge_attributes)
    return graph


def print_graph(graph, prefix, list_mxs_info):
    "Prints a graph in dot format"
    out_graph = prefix + ".mx.dot"
    outfile = open(out_graph, 'w')

    outfile.write("graph G {\n")

    # TODO: Make this more general
    colours = ["red", "green", "blue", "purple", "orange", "turquoise", "pink"]
    list_files = list(list_mxs_info.keys())

    for node in graph.nodes:
        files_labels = "\n".join([str(list_mxs_info[assembly][node]) for assembly in list_mxs_info])
        node_label = "\"%s\" [label=\"%s\n%s\"]" % (node, node, files_labels)
        outfile.write("%s\n" % node_label)

    for edge in graph.edges.data():
        outfile.write("\"%s\" -- \"%s\"" %
                      (edge[0], edge[1]))
        # For debugging only
        weight = edge[2]['support']
        if len(weight) == 1:
            colour = colours[list_files.index(weight[0])]
        elif len(weight) == 2:
            colour = "lightgrey"
        else:
            colour = "black"
        outfile.write(" [color=%s]\n" % colour)

    outfile.write("}\n")

    print("file_name\tnumber\tcolour")
    for i, filename in enumerate(list_files):
        print(filename, i, colours[i], sep="\t")


def determine_orientation(positions):
    "Given a list of minimizer positions, determine the orientation of the contig"
    if len(positions) == 1:
        return "?"
    if all(x < y for x, y in zip(positions, positions[1:])):
        return "+"
    if all(x > y for x, y in zip(positions, positions[1:])):
        return "-"
    return "?"


def calc_min_coord(positions, ctg_min_mx):
    "Calculates the minimum coordinate for a contig region in a path"
    if min(positions) == ctg_min_mx:
        return 0
    else:
        return min(positions)


def calc_max_coord(positions, ctg_max_mx, ctg_len):
    "Calculates the maximum coordinate for a contig region in a path"
    if max(positions) == ctg_max_mx:
        return ctg_len
    else:
        return max(positions)


def format_path(tuple_paths, mx_extremes, scaffolds):
    "Given a list of tuples (ctg, position), print out the order, orientation, and blocks of the contigs"
    out_path = []  # List of PathNode
    curr_ctg = None
    positions = []
    for tup in tuple_paths:
        if tup[0] is curr_ctg:
            positions.append(tup[1])
        else:
            # This is either the first tuple, or we are past a stretch of repeating contigs
            if curr_ctg is not None:
                ori = determine_orientation(positions)
                out_path.append(PathNode(contig=curr_ctg, ori=ori,
                                         start=calc_min_coord(positions, mx_extremes[curr_ctg][0]),
                                         end=calc_max_coord(positions, mx_extremes[curr_ctg][1],
                                                            len(scaffolds[curr_ctg]))))
            curr_ctg = tup[0]
            positions = [tup[1]]
    ori = determine_orientation(positions)
    out_path.append(PathNode(contig=curr_ctg, ori=ori,
                             start=calc_min_coord(positions, mx_extremes[curr_ctg][0]),
                             end=calc_max_coord(positions, mx_extremes[curr_ctg][1],
                                                len(scaffolds[curr_ctg]))))
    return out_path


def read_dot(dotfile_name):
    "Given a dot file, reads into a graph data structure"
    graph = nx.Graph(nx.drawing.nx_pydot.read_dot(dotfile_name))
    for _, _, eprop in graph.edges.data():
        if eprop['color'] == "black":
            eprop['weight'] = [1, 2, 3]
        elif eprop['color'] == 'lightgrey':
            eprop['weight'] = [1, 2]
        else:
            eprop['weight'] = [1]
    return graph


def filter_graph(graph, min_weight):
    "Filter the graph by edge weights"
    to_remove_edges = [(u, v) for u, v, e_prop in graph.edges.data() if e_prop['weight'] < min_weight]
    new_graph = graph.copy()
    new_graph.remove_edges_from(to_remove_edges)
    to_remove_nodes = [u for u in graph.nodes if graph.degree(u) > 2]
    new_graph.remove_nodes_from(to_remove_nodes)
    return new_graph


def find_paths(graph, list_mx_info, mx_extremes, scaffolds):
    "Finds paths per input assembly file"
    paths = {}
    skipped, total = 0, 0
    for assembly in list_mx_info:
        paths[assembly] = []
    for component in nx.connected_components(graph):
        component_graph = nx.subgraph(graph, component)
        source_nodes = [node for node in component_graph.nodes if component_graph.degree(node) == 1]
        if len(source_nodes) == 2:
            path = nx.shortest_path(component_graph, source_nodes[0], source_nodes[1])
            num_edges = len(path) - 1
            if len(path) == len(component_graph.nodes()) and \
                    num_edges == len(component_graph.edges()) and len(path) == len(set(path)):
                # All the nodes/edges from the graph are in the simple path, no repeated nodes
                for assembly in list_mx_info:
                    file_name, list_mx = assembly, list_mx_info[assembly]
                    tuple_paths = [list_mx[mx] for mx in path]
                    ctg_path = format_path(tuple_paths, mx_extremes[assembly], scaffolds[assembly])
                    paths[file_name].append(ctg_path)
                total += 1
            else:
                print("WARNING: Component with node", list(component.nodes)[0],
                      "was skipped.", sep=" ")
                skipped += 1
        else:
            print("WARNING: Component with node", list(component_graph.nodes)[0], "was skipped.", sep=" ")
            skipped += 1

    print("Warning: ", skipped, " paths of", total ,"were skipped", sep=" ")

    return paths


def read_fasta_file(filename):
    "Read a fasta file into memory"
    scaffolds = {}
    with open(filename, 'r') as fasta:
        for header, seq, _, _ in read_fasta(fasta):
            scaffolds[header] = seq
    return scaffolds


def reverse_complement(sequence):
    "Reverse complements a given sequence"
    mappings = {"A": "T", "C": "G", "G": "C", "T": "A",
                "M": "K", "R": "Y", "W": "W", "S": "S",
                "Y": "R", "K": "M", "V": "B", "H": "D",
                "D": "H", "B": "V", "N": "N"}
    new_sequence = ""
    for char in reversed(sequence):
        new_sequence += mappings[char.upper()]
    return new_sequence

def get_fasta_segment(path_node, sequence, k):
    "Given a PathNode, and the contig sequence, return the segment with the right bounds and orientation"
    if path_node.ori == "-":
        return reverse_complement(sequence[path_node.start:path_node.end+k+1])
    else:
        return sequence[path_node.start:path_node.end+k+1]


def format_bedtools_genome(scaffolds):
    "Format a BED file and genome dictionary for bedtools"
    bed_str = "\n".join(["%s\t%d\t%d" % (scaffold, 0, len(scaffolds[scaffold])) for scaffold in scaffolds])
    bed = pybedtools.BedTool(bed_str, from_string=True)

    genome_dict = {scaffold: (0, len(scaffolds[scaffold])) for scaffold in scaffolds}

    return bed, genome_dict


def print_scaffolds(paths, prefix, gap_size, k):
    "Given the paths, print out the scaffolds fasta"
    pathfile = open(prefix + ".path", 'w')

    for assembly in paths:
        min_match = re.search(r'^(\S+)\.tsv', assembly)
        assembly_fa = min_match.group(1)
        outfile = open(assembly_fa + ".scaffolds.fa", 'w')
        all_scaffolds = read_fasta_file(assembly_fa)
        incorporated_segments = []  # List of Bed entries

        ct = 0
        pathfile.write(assembly + "\n")
        for path in paths[assembly]:
            sequences = []
            path_segments = []
            for node in path:
                if node.ori == "?":
                    continue
                sequences.append(get_fasta_segment(node, all_scaffolds[node.contig], k))
                path_segments.append(Bed(contig=node.contig, start=node.start,
                                                 end=node.end))
            if len(sequences) < 2:
                continue
            outfile.write(">%s\n%s\n" %
                          ("mx" + str(ct),
                           ("N"*gap_size).join([seq for seq in sequences])))
            incorporated_segments.extend(path_segments)
            pathfile.write("%s\t%s\n" % ("mx" + str(ct),
                                         " ".join(["%s%s:%d-%d" % (tup.contig, tup.ori, tup.start, tup.end)
                                                   for tup in path])))
            ct += 1

        # Also print out the sequences that were NOT scaffolded
        incorporated_segments_str = "\n".join(["%s\t%s\t%s" % (chr, s, e) for chr, s, e in incorporated_segments])
        incorporated_segments_bed = pybedtools.BedTool(incorporated_segments_str, from_string=True).sort()
        genome_bed, genome_dict = format_bedtools_genome(all_scaffolds)

        missing_bed = genome_bed.complement(i=incorporated_segments_bed, g=genome_dict)
        missing_bed.saveas(prefix + ".unassigned.bed")

        cmd = "bedtools getfasta -fi %s -bed %s" % (assembly_fa, prefix + ".unassigned.bed")
        cmd_shlex = shlex.split(cmd)

        out_fasta = subprocess.Popen(cmd_shlex, stdout=subprocess.PIPE, universal_newlines=True)
        for line in iter(out_fasta.stdout.readline, ''):
            outfile.write(line)

        outfile.close()
    pathfile.close()


def find_mx_min_max(list_mx_info, graph):
    "Given a dictionary in the form assembly -> mx -> (ctg, pos), find the min and max mx position per ctg"
    mx_extremes = {} # assembly -> ctg -> (min_pos, max_pos)
    for assembly in list_mx_info:
        mx_extremes[assembly] = {}
        for mx in list_mx_info[assembly]:
            if mx not in graph:
                continue
            ctg, pos = list_mx_info[assembly][mx]
            if ctg in mx_extremes[assembly]:
                mx_extremes[assembly][ctg] = (min(mx_extremes[assembly][ctg][0], pos),
                                              max(mx_extremes[assembly][ctg][1], pos))
            else:
                mx_extremes[assembly][ctg] = (pos, pos)
    return mx_extremes


def main():
    "Run minimizer scaffolder"
    parser = argparse.ArgumentParser(
        description="Scaffold multiple genome assemblies using minimizers")
    parser.add_argument("FILES", nargs="+", help="Minimizer TSV files")
    parser.add_argument("-p", help="Output prefix [out]", default="out",
                        type=str, required=False)
    parser.add_argument("-g", help="Gap size [50]", default=50, type=int)
    parser.add_argument("-n", help="Minimum edge weight [2]", default=2, type=int)
    parser.add_argument("-l", help="List of assembly weights", required=True, type=str)
    parser.add_argument("-k", help="k value used for minimizering", required=True, type=int)
    args = parser.parse_args()

    # Parse the weights
    input_weights = re.split(r'\s+', args.l)
    if len(input_weights) != len(args.FILES):
        print("ERROR: The length of supplied weights and number of assembly files must be equal.")
        print("Supplied lengths:", len(input_weights), len(args.FILES), sep=" ")
        sys.exit(1)

    # Read in the minimizers
    list_mx_info = {}  # Dictionary of dictionaries of form mx -> (ctg, pos)
    list_mxs = {}  # Dictionary of Lists of minimizers (1 list per assembly file)
    weights = {}  # Dictionary of form file -> weight
    for assembly in args.FILES:
        mxs_info, mxs = read_minimizers(assembly)
        list_mx_info[assembly] = mxs_info
        list_mxs[assembly] = mxs
        weights[assembly] = float(input_weights.pop(0))

    # Filter minimizers - Keep only if unique in an assembly and found in all assemblies
    list_mxs = filter_minimizers(list_mxs)

    # Build a graph: Nodes = mx; Edges between adjacent mx in the assemblies
    graph = build_graph(list_mxs, weights)

    # Print the DOT graph
    print_graph(graph, args.p + "-before", list_mx_info)

    # Filter the graph edges + nodes
    graph = filter_graph(graph, args.n)

    # Print the DOT graph
    print_graph(graph, args.p, list_mx_info)

    # Find the min and max pos of minimizers per assembly, per ctg
    mx_extremes = find_mx_min_max(list_mx_info, graph)
    print(mx_extremes)

    # Load scaffolds into memory
    scaffolds = {}
    for assembly in args.FILES:
        min_match = re.search(r'^(\S+)\.tsv', assembly) # TODO: Make this more general
        assembly_fa = min_match.group(1)
        scaffolds[assembly] = read_fasta_file(assembly_fa)

    paths = find_paths(graph, list_mx_info, mx_extremes, scaffolds)

    print_scaffolds(paths, args.p, args.g, args.k)


if __name__ == "__main__":
    main()
