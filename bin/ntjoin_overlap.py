#!/usr/bin/env python3
"""
Detect and trim overlapping adjacent sequences
"""
from collections import namedtuple
import numpy as np
import ntjoin_utils
import ntjoin_assemble


MappedPathInfo = namedtuple("MappedPathInfo",
                            ["mapped_region_length", "mid_mx", "median_length_from_end"])

def merge_overlapping(list_mxs, list_mx_info, source, target, nodes):
    "Find the cut points for overlapping adjacent contigs"

    weights = {source: 1, target: 1}
    list_mxs_pair = {source: list_mxs[source], target: list_mxs[target]}
    list_mxs_pair = filter_minimizers_position(list_mxs_pair, source, target, list_mx_info, nodes)

    with ntjoin_utils.HiddenPrints():
        graph = ntjoin_assemble.Ntjoin.build_graph(ntjoin_assemble.Ntjoin(), list_mxs_pair, weights)
    graph = filter_graph_global(graph, 2)

    paths_components = []
    for component in graph.components():
        component_graph = graph.subgraph(component)
        source_nodes = [node.index for node in component_graph.vs() if node.degree() == 1]
        singleton_nodes = [node.index for node in component_graph.vs() if node.degree() == 0]
        if len(source_nodes) == 2:
            source_node, target_node = source_nodes
            if ntjoin_assemble.Ntjoin.vertex_name(component_graph, source_node) > \
                    ntjoin_assemble.Ntjoin.vertex_name(component_graph, target_node):
                source_node, target_node = target_node, source_node
            paths = component_graph.get_shortest_paths(source_node, target_node)
            assert len(paths) == 1
            path = [ntjoin_assemble.Ntjoin.vertex_name(component_graph, mx) for mx in paths[0]]
            start_mx, end_mx = path[0], path[-1]
            source_start, target_start = [list_mx_info[assembly][start_mx]
                                          for assembly in [source, target]]
            source_end, target_end = [list_mx_info[assembly][end_mx]
                                      for assembly in [source, target]]
            source_align_len = abs(source_start - source_end)
            target_align_len = abs(target_start - target_end)

            mid_mx = path[int(len(path)/2)]
            mid_mx_dist_end_source = get_dist_from_end(source,
                                                       list_mx_info[source][mid_mx],
                                                       nodes[source].get_aligned_length())
            mid_mx_dist_end_target = get_dist_from_end(target,
                                                       list_mx_info[target][mid_mx],
                                                       nodes[target].get_aligned_length(), target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=np.median([source_align_len,
                                                                                   target_align_len]),
                                                   mid_mx=mid_mx,
                                                   median_length_from_end=np.median(
                                                       [mid_mx_dist_end_source, mid_mx_dist_end_target])))
        elif singleton_nodes:
            assert len(singleton_nodes) == 1
            mid_mx = ntjoin_assemble.Ntjoin.vertex_name(component_graph, singleton_nodes[0])
            mid_mx_dist_end_source = get_dist_from_end(source, list_mx_info[source][mid_mx],
                                                       nodes[source].get_aligned_length())
            mid_mx_dist_end_target = get_dist_from_end(target, list_mx_info[target][mid_mx],
                                                       nodes[target].get_aligned_length(), target=True)
            paths_components.append(MappedPathInfo(mapped_region_length=1, mid_mx=mid_mx,
                                                   median_length_from_end=np.median([mid_mx_dist_end_source,
                                                                                     mid_mx_dist_end_target])))
        else:
            print("NOTE: non-singleton, {} source nodes".format(len(source_nodes)))
    if not paths_components:
        return False
    path = sorted(paths_components, key=lambda x: (x.mapped_region_length, x.median_length_from_end,
                                                   x.mid_mx), reverse=True)[0]
    source_cut, target_cut = list_mx_info[source][path.mid_mx], list_mx_info[target][path.mid_mx]

    if source_cut is None or target_cut is None:
        return False

    nodes[source].end_adjust = source_cut
    nodes[target].start_adjust = target_cut

    return True

def is_in_valid_region(pos, index, nodes):
    "Return true if the position is in the valid overlap region, else False"
    if index > 0 and pos < nodes[index-1].raw_gap_size*-1:
        return True
    if pos >= nodes[index].get_aligned_length() + nodes[index].raw_gap_size:
        return True
    return False

def is_in_valid_end(pos, index, nodes, source=True):
    "Return true if the mx is in a valid end"
    if source:
        return pos >= nodes[index].get_aligned_length() + nodes[index].raw_gap_size
    return pos < nodes[index].raw_gap_size*-1

def filter_minimizers_position(list_mxs_pair, source, target,
                               list_mx_info, nodes):
    "Filter to keep minimizers in particular positions"
    list_mxs_pair_return = {source: [[mx for mx in list_mxs_pair[source][0]
                                     if is_in_valid_end(list_mx_info[source][mx], source, nodes, source=True)]],
                            target: [[mx for mx in list_mxs_pair[target][0]
                                     if is_in_valid_end(list_mx_info[target][mx], source, nodes, source=False)]]}

    with ntjoin_utils.HiddenPrints():
        list_mxs_pair_return = ntjoin_utils.filter_minimizers(list_mxs_pair_return)

    return list_mxs_pair_return


def filter_graph_global(graph, n):
    "Filter the graph globally based on minimum edge weight"
    to_remove_edges = [edge.index for edge in graph.es()
                       if edge['weight'] < n]
    new_graph = graph.copy()
    new_graph.delete_edges(to_remove_edges)
    return new_graph

def get_dist_from_end(ori, pos, scaf_len, target=False):
    "Given the orientation, calculate the dist of the mx from the scaffold end (return -ve value)"
    if (ori == "+" and not target) or (ori == "-" and target):
        return (scaf_len - pos)*-1
    return pos*-1

def calc_total_weight(list_files, weights):
    "Calculate the total weight of an edge given the assembly support"
    return sum([weights[f] for f in list_files])

def set_edge_attributes(graph, edge_attributes):
    "Sets the edge attributes for a python-igraph graph"
    graph.es()["support"] = [edge_attributes[e]['support'] for e in sorted(edge_attributes.keys())]
    graph.es()["weight"] = [edge_attributes[e]['weight'] for e in sorted(edge_attributes.keys())]
