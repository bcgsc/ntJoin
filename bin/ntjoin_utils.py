#!/usr/bin/env python3
"""
ntJoin helper functions
Written by Lauren Coombe (@lcoombe)
"""

import datetime
from collections import namedtuple
import sys
import os


# Defining namedtuples
Bed = namedtuple("Bed", ["contig", "start", "end"])
Agp = namedtuple("Unassigned_bed", ["new_id", "contig", "start", "end"])
Scaffold = namedtuple("Scaffold", ["id", "length", "sequence"])
EdgeGraph = namedtuple("EdgeGraph", ["source", "target", "raw_gap_est"])

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

# Helper functions
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


# Defining helper classes
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

class OrientationError(ValueError):
    "Orientation type error"
    def __init__(self):
        self.message = "Orientation must be + or -"
        super().__init__(self.message)

class OverlapRegion:
    "Overlapping regions in a contig to fix"
    def __init__(self):
        self.regions = []
        self.best_region = None

    def add_region(self, bed_region):
        "Add a new region to the overlapping set"
        if self.best_region is None or \
                        (bed_region.end - bed_region.start) > \
                        (self.best_region.end - self.best_region.start):
            self.best_region = bed_region
        self.regions.append(bed_region)
        assert bed_region.contig == self.best_region.contig

    @staticmethod
    def are_overlapping(region1, region2):
        "Returns True if the given regions are overlapping"
        return region1.start <= region2.end and region2.start <= region1.end

    @staticmethod
    def is_subsumed(region1, region2):
        "Returns True is region 1 is subsumed in region2"
        return region1.start >= region2.start and region1.end <= region2.end

    def find_non_overlapping(self):
        "Given overlapping regions, resolve to remove overlaps"
        return_regions = {} # Bed -> (replacement Bed) or None
        if not self.regions or self.best_region is None:
            return None
        for bed_region in self.regions:
            if bed_region == self.best_region:
                return_regions[bed_region] = bed_region
            elif self.is_subsumed(bed_region, self.best_region):
                # Subsumed region in the best region
                return_regions[bed_region] = None
            elif self.are_overlapping(bed_region, self.best_region):
                # Overlaps with best region, but isn't subsumed
                if bed_region.start <= self.best_region.start:
                    new_region = Bed(contig=bed_region.contig, start=bed_region.start,
                                     end=self.best_region.start - 1)
                elif bed_region.end >= self.best_region.end:
                    new_region = Bed(contig=bed_region.contig, start=self.best_region.end + 1,
                                     end=bed_region.end)
                return_regions[bed_region] = new_region
            else:
                return_regions[bed_region] = bed_region

        # Double check if any still overlaps. If so, adjust smaller of the overlapping regions.
        existing_overlaps = True
        while existing_overlaps:
            sorted_regions = sorted([(b, a) for b, a in return_regions.items() if a is not None], key=lambda x: x[1])
            i, j = 0, 1
            existing_overlaps = False
            while j < len(sorted_regions):
                region1_before, region2_before = sorted_regions[i][0], sorted_regions[j][0]
                region1_after, region2_after = sorted_regions[i][1], sorted_regions[j][1]
                if region1_after is None or region2_after is None:
                    i += 1
                    j += 1
                    continue
                if self.are_overlapping(region1_after, region2_after):
                    existing_overlaps = True
                    if self.is_subsumed(region1_after, region2_after):
                        # Region 1 is subsumed in region 2 - Remove region 1
                        return_regions[region1_before] = None
                    elif self.is_subsumed(region2_after, region1_after):
                        # Region 2 is subsumed in region 1 - Remove region 2
                        return_regions[region2_before] = None
                    elif (region1_after.end - region1_after.start) > (region2_after.end - region2_after.start):
                        # Adjust region 2 start
                        return_regions[region2_before] = Bed(contig=region2_after.contig, start=region1_after.end + 1,
                                                             end=region2_after.end)
                    elif (region1_after.end - region1_after.start) <= (region2_after.end - region2_after.start):
                        # Adjust region 1 end
                        return_regions[region1_before] = Bed(contig=region1_after.contig, start=region1_after.start,
                                                             end=region2_after.start - 1)
                    else:
                        print("Unexpected case!")
                        print(region1_before, region2_before, region1_after, region1_after)

                i += 1
                j += 1

        return return_regions
