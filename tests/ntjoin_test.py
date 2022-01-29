"""Tests for ntJoin"""

import shlex
import subprocess
import re

def launch_ntjoin(cmd, prefix):
    cmd_shlex = shlex.split(cmd)
    return_code = subprocess.call(cmd_shlex)
    assert return_code == 0
    return_paths = []
    with open(prefix + ".path", 'r') as paths:
        for line in paths:
            path_match = re.search(r'^ntJoin', line)
            if path_match:
                return_paths.append(line.strip())
    return return_paths

def run_ntjoin(ref1, target, prefix, window=1000, n=2):
    "Run ntJoin with a pair of files"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 references=\'{ref}\' reference_weights=\'2\' " \
          "prefix={prefix} k=32 w={w} n={n}".format(target=target, ref=ref1, prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

def run_ntjoin_nocut(ref1, target, prefix, window=1000, n=2):
    "Run ntJoin with a pair of files"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 references=\'{ref}\' reference_weights=\'2\' " \
          "prefix={prefix} k=32 w={w} n={n} no_cut=True".format(target=target, ref=ref1, prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

def run_ntjoin_multiple(ref1, ref2, target, prefix, window=1000, n=2):
    "Run ntJoin with a target and 2 references"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 references=\'{ref1} {ref2}\' reference_weights=\'2 2\' " \
          "prefix={prefix} k=32 w={w} n={n}".format(target=target, ref1=ref1, ref2=ref2, prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

def run_ntjoin_agp(ref1, target, prefix, window=1000, n=2):
    "Run ntJoin with a pair of files"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 references=\'{ref}\' reference_weights=\'2\' " \
          "prefix={prefix} k=32 w={w} n={n} agp=True".format(target=target, ref=ref1, prefix=prefix, w=window, n=n)
    cmd_shlex = shlex.split(cmd)
    return_code = subprocess.call(cmd_shlex)
    assert return_code == 0
    return_agp = []
    with open(prefix + ".agp", 'r') as agp:
        for line in agp:
            return_agp.append(line.strip())
    return return_agp

def run_ntjoin_config(config_file, target, prefix, window=1000, n=2):
    "Run ntJoin with a target and reference config file"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 " \
          "reference_config={config} prefix={prefix} k=32 w={w} n={n}".format(target=target, config=config_file,
                                                                              prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

def run_ntjoin_config_extra(config_file, target, prefix, window=1000, n=2):
    "Run ntJoin with a target and reference config file, overriding reference and reference_weights variables"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 reference=na reference_weights=na " \
          "reference_config={config} prefix={prefix} k=32 w={w} n={n}".format(target=target, config=config_file,
                                                                              prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

def run_ntjoin_overlap(ref1, target, prefix, window=1000, n=2):
    "Run ntJoin with a pair of files"
    cmd = "../ntJoin assemble -B target={target} target_weight=1 references=\'{ref}\' reference_weights=\'2\' " \
          "prefix={prefix} k=32 w={w} n={n} overlap=True".format(target=target, ref=ref1, prefix=prefix, w=window, n=n)
    return_paths = launch_ntjoin(cmd, prefix)
    return return_paths

# Following 4 tests to check for the expected PATHs given 2 pieces that should be merged
#     together based on the reference in different orientations
#     - Pieces are the reference piece split, with ~20bp in between

def test_mx_f_f():
    "Testing ntJoin with assembly + reference, fwd-fwd orientation"
    paths = run_ntjoin("ref.fa", "scaf.f-f.fa", "f-f_test")
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:0-1981 20N 2_f+:0-2329"

def test_mx_f_f_termN():
    "Testing stripping terminal Ns in output scaffolds"
    paths = run_ntjoin("ref.fa", "scaf.f-f.termN.fa", "f-f_test")
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:4-1985 20N 2_f+:0-2329"

def test_mx_f_r():
    "Testing ntJoin with assembly + reference, fwd-rev orientation"
    paths = run_ntjoin("ref.fa", "scaf.f-r.fa", "f-r_test")
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:0-1981 20N 2_r-:0-2329"


def test_mx_r_f():
    "Testing ntJoin with assembly + reference, rev-fwd orientation"
    paths = run_ntjoin("ref.fa", "scaf.r-f.fa", "r-f_test")
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_r-:0-1981 20N 2_f+:0-2329"


def test_mx_r_r():
    "Testing ntJoin with assembly + reference, rev-rev orientation"
    paths = run_ntjoin("ref.fa", "scaf.r-r.fa", "r-r_test")
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_r-:0-1981 20N 2_r-:0-2329"

# Test checks for the expected gap length and sequence orientation for a
# test with 2 expected output paths
def test_gap_dist_multiple():
    "Testing ntJoin with assembly + reference, estimated gap length"
    paths = run_ntjoin("ref.multiple.fa", "scaf.multiple.fa", "gap-dist_test", window=500, n=1)
    assert len(paths) == 2
    assert paths[0] != paths[1]
    expected_paths = ["2_1_p+:0-2492 100N 2_2_n-:0-2574", "1_1_p+:0-1744 124N 1_2_p+:0-1844"]
    assert paths.pop().split("\t")[1] in expected_paths
    assert paths.pop().split("\t")[1] in expected_paths


# Tests for gap distance estimation, misassembled scaffolds
# Expected that 2 input scaffolds will be broken and joined based on the reference.
# Testing orientations of joins: +/+ -/- +/- -/+
def test_regions_ff_rr():
    "Testing ntJoin correcting misassemblies, joins in fwd-fwd and rev-rev"
    paths = run_ntjoin("ref.multiple.fa", "scaf.misassembled.f-f.r-r.fa", "regions-ff-rr_test", window=500, n=1)
    assert len(paths) == 2
    assert paths[0] != paths[1]
    expected_paths = ["2_1n-1_2p-:0-2176 20N 1_1p-2_2n-:2010-4489", "1_1p-2_2n+:0-1541 468N 2_1n-1_2p+:2676-4379"]
    assert paths.pop().split("\t")[1] in expected_paths
    assert paths.pop().split("\t")[1] in expected_paths

def test_regions_ff_rr_nocut():
    "Testing ntJoin correcting misassemblies, joins in fwd-fwd and rev-rev"
    paths = run_ntjoin_nocut("ref.multiple.fa", "scaf.misassembled.f-f.r-r.fa", "regions-ff-rr-nocut_test", window=500, n=1)
    assert len(paths) == 1
    assert paths[0].split("\t")[1] == "2_1n-1_2p-:0-4379 20N 1_1p-2_2n-:0-4489"

def test_regions_fr_rf():
    "Testing ntJoin correcting misassemblies, joins in fwd-rev and rev-fwd"
    paths = run_ntjoin("ref.multiple.fa", "scaf.misassembled.f-r.r-f.fa", "regions-fr-rf_test", 500, n=2)
    assert len(paths) == 2
    assert paths[0] != paths[1]
    expected_paths = ["2_1n-1_2n-:0-2176 212N 1_1p-2_2p+:2017-4489", "1_1p-2_2p+:0-1617 198N 2_1n-1_2n-:2675-4379"]
    assert paths.pop().split("\t")[1] in expected_paths
    assert paths.pop().split("\t")[1] in expected_paths

def test_regions_fr_rf_config():
    "Testing ntJoin correcting misassemblies, joins in fwd-rev and rev-fwd, using config file"
    paths = run_ntjoin_config("test_config_single.csv", "scaf.misassembled.f-r.r-f.fa", "regions-fr-rf_test", 500, n=2)
    assert len(paths) == 2
    assert paths[0] != paths[1]
    expected_paths = ["2_1n-1_2n-:0-2176 212N 1_1p-2_2p+:2017-4489", "1_1p-2_2p+:0-1617 198N 2_1n-1_2n-:2675-4379"]
    assert paths.pop().split("\t")[1] in expected_paths
    assert paths.pop().split("\t")[1] in expected_paths

def test_regions_3():
    "Testing ntJoin with target + 2 references"
    paths = run_ntjoin_multiple("ref.fa", "scaf.f-f.copy.fa", "scaf.f-f.fa", "f-f-f_test", n=1)
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:0-1981 20N 2_f+:0-2329"

def test_regions_3_config():
    "Testing ntJoin with target + 2 references, using config file"
    paths = run_ntjoin_config("test_config_multiple.csv", "scaf.f-f.fa", "f-f-f_test", n=1)
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:0-1981 20N 2_f+:0-2329"

def test_regions_3_config_extra():
    "Testing ntJoin with target + 2 references, using config file, command having extra parameters"
    paths = run_ntjoin_config_extra("test_config_multiple.csv", "scaf.f-f.fa", "f-f-f_test", n=1)
    assert len(paths) == 1
    assert paths.pop() == "ntJoin0\t1_f+:0-1981 20N 2_f+:0-2329"

# Testing AGP output
def test_mx_r_f_agp():
    "Testing ntJoin with assembly + reference, rev-fwd orientation - AGP output"
    agp = run_ntjoin_agp("ref.fa", "scaf.r-f.fa", "r-f_test")
    assert len(agp) == 3
    assert agp[0].strip() == "ntJoin0\t1\t1981\t1\tW\t1_r\t1\t1981\t-"
    assert agp[1].strip() == "ntJoin0\t1982\t2001\t2\tN\t20\tscaffold\tyes\talign_genus"
    assert agp[2].strip() == "ntJoin0\t2002\t4330\t3\tW\t2_f\t1\t2329\t+"

# Testing AGP output
def test_mx_f_f_agp():
    "Testing ntJoin with assembly + reference, fwd-fwd orientation, with terminal gaps - AGP output + unassigned"
    agp = run_ntjoin_agp("ref.fa", "scaf.f-f.termN.unassigned.fa", "f-f_test")
    assert len(agp) == 4
    assert agp[0].strip() == "ntJoin0\t1\t1981\t1\tW\t1_f\t5\t1985\t+"
    assert agp[1].strip() == "ntJoin0\t1982\t2001\t2\tN\t20\tscaffold\tyes\talign_genus"
    assert agp[2].strip() == "ntJoin0\t2002\t4330\t3\tW\t2_f\t1\t2329\t+"
    assert agp[3].strip() == "unassigned:0-14\t1\t8\t1\tW\tunassigned\t3\t10\t+"

# Testing overlap code
def test_mx_f_f_overlap():
    "Testing ntJoin with assembly + reference, fwd-fwd orientation, overlap code on"
    paths = run_ntjoin_overlap("ref.fa", "scaf.f-f.overlapping.fa", "f-f_test_overlap")
    assert paths.pop() == "ntJoin0\t1+:0-2037 20N 2+:38-2331"

def test_mx_f_r_overlap():
    "Testing ntJoin with assembly + reference, fwd-rev orientation, overlap code on"
    paths = run_ntjoin_overlap("ref.fa", "scaf.f-r.overlapping.fa", "f-r_test_overlap")
    assert paths.pop() == "ntJoin0\t1+:0-2037 20N 2-:0-2293"

def test_mx_r_r_overlap():
    "Testing ntJoin with assembly + reference, rev-rev orientation, overlap code on"
    paths = run_ntjoin_overlap("ref.fa", "scaf.r-r.overlapping.fa", "f-r_test_overlap")
    assert paths.pop() == "ntJoin0\t1-:62-2099 20N 2-:0-2293"
