#!/usr/bin/env python3

'''
Representing a ntJoin node in a scaffolding path
'''

class OrientationError(ValueError):
    "Orientation type error"
    def __init__(self):
        self.message = "Orientation must be + or -"
        super().__init__(self.message)

class PathNode:
    "Defines a node in a path of contig regions"
    def __init__(self, contig, ori, start, end, contig_size,
                 first_mx, terminal_mx, gap_size=0, raw_gap_size=0):
        self.contig = contig
        self.ori = ori
        self.start = start
        self.end = end
        self.contig_size = contig_size
        self.first_mx = first_mx
        self.terminal_mx = terminal_mx
        self.gap_size = gap_size
        self.raw_gap_size = raw_gap_size
        self.start_adjust = 0
        self.end_adjust = 0  # Adjust for trimming

    def set_gap_size(self, gap_size):
        "Set the gap size of the path node"
        self.gap_size = gap_size

    def set_raw_gap_size(self, raw_gap_size):
        "Set the 'raw' gap size. Equal to gap_size if > min_gap_size"
        self.raw_gap_size = raw_gap_size

    def get_aligned_length(self):
        "Get the aligned length based on start/end coordinates"
        return self.end - self.start

    def get_end_adjusted_coordinate(self):
        "Return the adjusted end coordinate"
        if self.end_adjust == 0:
            return self.get_aligned_length()
        return self.end_adjust

    def get_adjusted_start(self):
        "Return the start coordinate of segment, adjusted for any trimming"
        if self.ori == "+":
            return self.start + self.start_adjust
        if self.ori == "-":
            return self.start + (self.get_aligned_length() - self.get_end_adjusted_coordinate())
        raise OrientationError()

    def get_adjusted_end(self):
        "Return the end coordinate of the segment, adjusted for any trimming"
        if self.ori == "+":
            return self.end - (self.get_aligned_length() - self.get_end_adjusted_coordinate())
        if self.ori == "-":
            return self.end - self.start_adjust
        raise OrientationError()

    def __str__(self):
        return f"Contig:{self.contig}\tOrientation:{self.ori}\tStart-End:{self.start}-{self.end}\t"\
                f"Length:{self.contig_size}\tFirstmx:{self.first_mx}\tLastmx:{self.terminal_mx}\t" \
                f"Adjusted_start-end:{self.start_adjust}-{self.end_adjust}"
