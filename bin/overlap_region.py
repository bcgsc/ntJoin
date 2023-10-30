#!/usr/bin/env python3
'''
Representing an overlapping region in a contig to resolve
'''
from ntjoin_utils import Bed

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
