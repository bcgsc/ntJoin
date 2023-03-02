#!/usr/bin/env python3
'''
Represents an ntJoin synteny block
'''

from assembly_block import AssemblyBlock
from ntjoin_utils import Minimizer

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
