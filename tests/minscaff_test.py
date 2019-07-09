import shlex
import subprocess
import re
import sys
import os

def run_ntJoin(file1, file2, prefix):
    cmd = "../minimizer_assembler-make assemble -B list_files=\'" + file1 + " " + file2 + "\' " \
          "list_weights=\'2 1\' k=32 w=1000 n=2 prefix=" + prefix
    cmd_shlex = shlex.split(cmd)
    return_code = subprocess.call(cmd_shlex)
    assert return_code == 0
    with open(prefix + ".path", 'r') as paths:
        for line in paths:
            path_match = re.search(r'^mx', line)
            if path_match:
                return_line = line.strip()
    return return_line

def test_mx_f_f():
    path = run_ntJoin("ref.fa", "scaf.f-f.fa", "f-f_test")
    assert path == "mx0\t1_f+:0-1981 20N 2_f+:0-2329"


def test_mx_f_r():
    path = run_ntJoin("ref.fa", "scaf.f-r.fa", "f-r_test")
    assert path == "mx0\t1_f+:0-1981 20N 2_r-:0-2329"


def test_mx_r_f():
    path = run_ntJoin("ref.fa", "scaf.r-f.fa", "r-f_test")
    assert path == "mx0\t1_r-:0-1981 20N 2_f+:0-2329"


def test_mx_r_r():
    path = run_ntJoin("ref.fa", "scaf.r-r.fa", "r-f_test")
    assert path == "mx0\t1_r-:0-1981 20N 2_r-:0-2329"
