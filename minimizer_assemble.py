#!/usr/bin/env python3
"""
Scaffold multiple genome assemblies using minimizers
Written by Lauren Coombe
"""

import argparse
import re
from collections import defaultdict
import networkx as nx



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


def build_graph(list_mxs):
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

    formatted_edges = [(s, t, edges[s][t]) for s in edges for t in edges[s]]

    graph.add_weighted_edges_from(formatted_edges)
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
        weight = edge[2]['weight']
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


def format_path(tuple_paths):
    "Given a list of tuples (ctg, position), print out the order and orientation of the contigs"
    out_path = [] # List of (ctg, orientation)
    curr_ctg = None
    positions = []
    for tup in tuple_paths:
        if tup[0] is curr_ctg:
            positions.append(tup[1])
        else:
            # This is either the first tuple, or we are past a stretch of repeating contigs
            if curr_ctg is not None:
                ori = determine_orientation(positions)
                out_path.append((curr_ctg, ori))
            curr_ctg = tup[0]
            positions = [tup[1]]
    ori = determine_orientation(positions)
    out_path.append((curr_ctg, ori))
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
    to_remove_edges = [(u, v) for u, v, e_prop in graph.edges.data() if len(e_prop['weight']) < min_weight]
    new_graph = graph.copy()
    new_graph.remove_edges_from(to_remove_edges)
    to_remove_nodes = [u for u in graph.nodes if graph.degree(u) > 2]
    new_graph.remove_nodes_from(to_remove_nodes)
    return new_graph


def find_paths(graph, list_mx_info):
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
                    ctg_path = format_path(tuple_paths)
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


def read_fasta(filename):
    "Read a fasta file into memory"
    scaffolds = {}
    header = None
    header_match = re.compile(r'^>(\S+)')
    with open(filename, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            line_match = re.search(header_match, line)
            if line_match:
                header = line_match.group(1)
            else:
                scaffolds[header] = line
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


def print_scaffolds(paths, prefix, gap_size):
    "Given the paths, print out the scaffolds fasta"
    pathfile = open(prefix + ".path", 'w')

    for assembly in paths:
        min_match = re.search(r'^(\S+)\.tsv', assembly)
        assembly_fa = min_match.group(1)
        outfile = open(assembly_fa + ".scaffolds.fa", 'w')
        all_scaffolds = read_fasta(assembly_fa)

        scaffolded = set()  # Track the pieces incorporated into a scaffold
        ct = 0
        pathfile.write(assembly + "\n")
        for path in paths[assembly]:
            sequences = []
            for ctg, ori in path:
                if ori == "?":
                    continue
                elif ori == "+":
                    sequences.append((ctg, all_scaffolds[ctg]))
                else:
                    sequences.append((ctg, reverse_complement(all_scaffolds[ctg])))
                #scaffolded.add(ctg)
            if len(sequences) < 2:
                continue
            outfile.write(">%s\n%s\n" %
                          ("mx" + str(ct),
                           ("N"*gap_size).join([seq[1] for seq in sequences])))
            scaffolded = scaffolded.union(set([seq[0] for seq in sequences]))
            pathfile.write("%s\t%s\n" % ("mx" + str(ct), " ".join(["".join(tup) for tup in path])))
            ct += 1
        # Also print out the sequences that were NOT scaffolded
        for ctg in all_scaffolds:
            if ctg not in scaffolded:
                outfile.write(">%s\n%s\n" % (ctg, all_scaffolds[ctg]))

        outfile.close()
    pathfile.close()


def main():
    "Run minimizer scaffolder"
    parser = argparse.ArgumentParser(
        description="Scaffold multiple genome assemblies using minimizers")
    parser.add_argument("FILES", nargs="+", help="Minimizer TSV files")
    parser.add_argument("-p", help="Output prefix [out]", default="out",
                        type=str, required=False)
    parser.add_argument("-g", help="Gap size [50]", default=50, type=int)
    parser.add_argument("-w", help="Minimum number of assemblies supporting an edge [2]", default=2, type=int)
    args = parser.parse_args()

    # Read in the minimizers
    list_mx_info = {}  # Dictionary of dictionaries of form mx -> (ctg, pos)
    list_mxs = {}  # Dictionary of Lists of minimizers (1 list per assembly file)
    for assembly in args.FILES:
        #list_files.append(assembly)
        mxs_info, mxs = read_minimizers(assembly)
        list_mx_info[assembly] = mxs_info
        list_mxs[assembly] = mxs

    list_mxs = filter_minimizers(list_mxs)

    graph = build_graph(list_mxs)

    print_graph(graph, args.p + "-before", list_mx_info)

    graph = filter_graph(graph, args.w)

    print_graph(graph, args.p, list_mx_info)

    paths = find_paths(graph, list_mx_info)

    print_scaffolds(paths, args.p, args.g)


if __name__ == "__main__":
    main()
