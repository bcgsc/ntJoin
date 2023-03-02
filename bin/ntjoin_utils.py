#!/usr/bin/env python3
"""
ntJoin helper functions
Written by Lauren Coombe (@lcoombe)
"""

import datetime
from collections import namedtuple, defaultdict
import shlex
import subprocess
import sys
import os
import igraph as ig


# Defining namedtuples
Bed = namedtuple("Bed", ["contig", "start", "end"])
Agp = namedtuple("Unassigned_bed", ["new_id", "contig", "start", "end"])
Scaffold = namedtuple("Scaffold", ["id", "length", "sequence"])
EdgeGraph = namedtuple("EdgeGraph", ["source", "target", "raw_gap_est"])
Minimizer = namedtuple("Minimizer", ["mx", "position"])

class HiddenPrints:
    "Adapted from: https://stackoverflow.com/questions/8391411/how-to-block-calls-to-print"
    def __init__(self):
        self._original_stdout = sys.stdout

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w', encoding="utf-8")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# Helper functions for interfacing with python-igraph
def vertex_index(graph, name):
    "Returns vertex index based on vertex name"
    return graph.vs.find(name).index

def vertex_name(graph, index):
    "Returns vertex name based on vertex id"
    return graph.vs[index]['name']

def edge_index(graph, source_name, target_name):
    "Returns graph edge index based on source/target names"
    return graph.get_eid(source_name, target_name)

def set_edge_attributes(graph, edge_attributes):
    "Sets the edge attributes for a python-igraph graph"
    graph.es()["support"] = [edge_attributes[e]['support'] for e in sorted(edge_attributes.keys())]
    graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]

def calc_total_weight(list_files, weights):
    "Calculate the total weight of an edge given the assembly support"
    return sum([weights[f] for f in list_files])

def build_graph(list_mxs, weights, graph=None, black_list=None):
    "Builds an undirected graph: nodes=minimizers; edges=between adjacent minimizers"
    print(datetime.datetime.today(), ": Building graph", file=sys.stdout)

    if graph is None:
        graph = ig.Graph()
        prev_edge_attributes = {}
    else:
        prev_edge_attributes = {e.index: {"support": e['support'],
                                            "weight": e['weight']} for e in graph.es()}

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
                if black_list is None or assembly_mx_list[i] not in black_list:
                    vertices.add(assembly_mx_list[i])
            if assembly_mx_list:
                if black_list is None or assembly_mx_list[-1] not in black_list:
                    vertices.add(assembly_mx_list[-1])

    formatted_edges = [(s, t) for s in edges for t in edges[s]]

    print(datetime.datetime.today(), ": Adding vertices", file=sys.stdout)
    graph.add_vertices(list(vertices))

    print(datetime.datetime.today(), ": Adding edges", file=sys.stdout)
    graph.add_edges(formatted_edges)

    print(datetime.datetime.today(), ": Adding attributes", file=sys.stdout)
    edge_attributes = {edge_index(graph, s, t): {"support": edges[s][t],
                                                "weight": calc_total_weight(edges[s][t],
                                                                            weights)}
                        for s in edges for t in edges[s]}
    edge_attributes.update(prev_edge_attributes)
    set_edge_attributes(graph, edge_attributes)

    return graph

# Other helper functions

def reverse_complement(sequence):
    "Reverse complements a given sequence"
    translation_table = str.maketrans(
        "ACGTUNMRWSYKVHDBacgtunmrwsykvhdb",
        "TGCAANKYWSRMBDHVtgcaankywsrmbdhv")
    return sequence[::-1].translate(translation_table)

def filter_minimizers(list_mxs):
    "Filters out minimizers that are not found in all assemblies"
    print(datetime.datetime.today(), ": Filtering minimizers", file=sys.stdout)
    list_mx_sets = [{mx for mx_list in list_mxs[assembly] for mx in mx_list}
                    for assembly in list_mxs]
    mx_intersection = set.intersection(*list_mx_sets)

    return_mxs = {}
    for assembly in list_mxs:
        assembly_mxs_filtered = [[mx for mx in mx_list if mx in mx_intersection]
                                 for mx_list in list_mxs[assembly]]
        return_mxs[assembly] = assembly_mxs_filtered

    return return_mxs

def read_minimizers(tsv_filename):
    "Read the minimizers from a file, removing duplicate minimizers"
    print(datetime.datetime.today(), ": Reading minimizers", tsv_filename, file=sys.stdout)
    mx_info = {}  # mx -> (contig, position)
    mxs = []  # List of lists of minimizers
    dup_mxs = set()  # Set of minimizers identified as duplicates
    with open(tsv_filename, 'r', encoding="utf-8") as tsv:
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

    mx_info = {mx: mx_entry_info for mx, mx_entry_info in mx_info.items() if mx not in dup_mxs}
    mxs_filt = []
    for mx_list in mxs:
        mx_list_filt = [mx for mx in mx_list if mx not in dup_mxs]
        mxs_filt.append(mx_list_filt)
    return mx_info, mxs_filt

def run_indexlr(assembly, k, w, t):
    "Run indexlr on the given assembly with the specified k and w"
    cmd = f"indexlr {assembly} --long --pos -k{k} -w{w} -t{t} -o {assembly}.k{k}.w{w}.tsv"
    cmd = shlex.split(cmd)
    ret_code = subprocess.call(cmd)
    assert ret_code == 0
    return f"{assembly}.k{k}.w{w}.tsv"
