#!/usr/bin/env python3
"""
ntJoin: Scaffolding assemblies using reference assemblies and minimizer graphs
Written by Lauren Coombe (@lcoombe)
"""

import datetime
import multiprocessing
import sys
import warnings
import ntjoin_utils
warnings.simplefilter(action='ignore', category=RuntimeWarning)


class Ntjoin:
    "ntJoin: Scaffolding and analyzing synteny in assemblies using reference assemblies and minimizer graphs"

    @staticmethod
    def convert_path_index_to_name(graph, path):
        "Convert path of vertex indices to path of vertex names"
        return [ntjoin_utils.vertex_name(graph, vs) for vs in path]

    def print_graph(self, graph, out_prefix=None):
        "Prints the minimizer graph in dot format"
        if out_prefix is None:
            out_graph = self.args.p + ".mx.dot"
        else:
            out_graph = out_prefix + "mx.dot"

        with open(out_graph, 'w',  encoding="utf-8") as outfile:
            print(datetime.datetime.today(), ": Printing graph", out_graph, sep=" ", file=sys.stdout)

            outfile.write("graph G {\n")

            colours = ["red", "green", "blue", "purple", "orange",
                    "turquoise", "pink", "yellow", "orchid", "salmon"]
            list_files = list(self.list_mx_info.keys())
            if len(list_files) > len(colours):
                colours = ["red"]*len(list_files)

            for node in graph.vs():
                mx_ctg_pos_labels = "\n".join([f"{file_name}_{asm_mx_info[node['name']]}"
                                            for file_name, asm_mx_info in self.list_mx_info.items()])
                node_label = f"\"{node['name']}\" [label=\"{node['name']}\n{mx_ctg_pos_labels}\"]"
                outfile.write(f"{node_label}\n")

            for edge in graph.es():
                outfile.write(f"\"{ntjoin_utils.vertex_name(graph, edge.source)}\" --" \
                            f"\"{ntjoin_utils.vertex_name(graph, edge.target)}\"")
                weight = edge['weight']
                support = edge['support']
                if len(support) == 1:
                    colour = colours[list_files.index(support[0])]
                elif len(support) == 2:
                    colour = "lightgrey"
                else:
                    colour = "black"
                outfile.write(f" [weight={weight} color={colour}]\n")

            outfile.write("}\n")

        print("\nfile_name\tnumber\tcolour")
        for i, filename in enumerate(list_files):
            print(filename, i, colours[i], sep="\t")
        print("", flush=True)


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
        print(datetime.datetime.today(), ": Filtering the graph", file=sys.stdout, flush=True)
        if self.args.n <= min(self.weights.values()):
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
        max_wt_asm = [assembly for assembly, asm_weight in self.weights.items()
                      if asm_weight == max(self.weights.values())].pop()
        list_mx_info_maxwt = self.list_mx_info[max_wt_asm]
        min_pos = min((list_mx_info_maxwt[ntjoin_utils.vertex_name(graph, s)][1] for s in sources))
        max_pos = max((list_mx_info_maxwt[ntjoin_utils.vertex_name(graph, s)][1] for s in sources))
        source = [s for s in sources
                  if list_mx_info_maxwt[ntjoin_utils.vertex_name(graph, s)][1] == min_pos].pop()
        target = [s for s in sources
                  if list_mx_info_maxwt[ntjoin_utils.vertex_name(graph, s)][1] == max_pos].pop()
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
        max_edge_weight = sum(self.weights.values())
        component_graph = self.graph.subgraph(component)
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
                    return_paths.append((path, subcomponent_graph))
        return return_paths


    def find_paths(self):
        "Finds paths through the minimizer graph"
        print(datetime.datetime.today(), ": Finding paths", file=sys.stdout)
        components = self.graph.components()
        print("\nTotal number of components in graph:", len(components), "\n", sep=" ", file=sys.stdout, flush=True)

        if self.args.t == 1:
            paths = [self.find_paths_process(component) for component in components]
        else:
            with multiprocessing.Pool(self.args.t) as pool:
                paths = pool.map(self.find_paths_process, components)

        return paths

    def load_minimizers(self, repeat_bf=False):
        "Load in minimizers for ntJoin scaffolding mode"
        weights = {}  # Dictionary: assembly -> weight
        for assembly in self.args.FILES:
            mxs_info, mxs = ntjoin_utils.read_minimizers(assembly, repeat_bf)
            self.list_mx_info[assembly] = mxs_info
            self.list_mxs[assembly] = mxs
            weights[assembly] = self.weights_list.pop(0)
        self.weights = weights


    def make_minimizer_graph(self):
        "Run ntJoin graph stage"
        print(datetime.datetime.today(), ": Generating minimizer graph ...\n")

        # Print the weights of the input assemblies
        weight_str = "\n".join([f"{assembly}: {asm_weight}" for assembly, asm_weight in self.weights.items()])
        print("\nWeights of assemblies:\n", weight_str, "\n", sep="", flush=True)

        # Filter minimizers - Keep only if found in all assemblies
        list_mxs = ntjoin_utils.filter_minimizers(self.list_mxs)

        # Build a graph: Nodes = mx; Edges between adjacent mx in the assemblies
        self.graph = ntjoin_utils.build_graph(list_mxs, self.weights)

        # Print the DOT graph
        self.print_graph(self.graph)


    def ntjoin_find_paths(self):
        "Find the paths through the graph"
        paths = self.find_paths()
        return paths

    def __init__(self, args):
        "Create an ntJoin instance"
        self.list_mx_info = {}  # Dictionary of dictionaries: assembly -> mx -> (contig, position)
        self.list_mxs = {}  # Dictionary: assembly -> [lists of mx]
        self.graph = None
        self.args = args
        self.weights = {}
        self.weights_list = []
