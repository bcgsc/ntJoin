import shlex
import subprocess
import re
import sys
import os

def test_mx():
    cmd = "../minimizer_assembler-make assemble -B list_files=\'ref.fa scaf.fa\' " \
          "list_weights=\'2 1\' k=32 w=1000 n=2 prefix=out_scaf.k32.w1000.n2"
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    with open("out_scaf.k32.w1000.n2.path", 'r') as paths:
        for line in paths:
            path_match = re.search(r'^mx', line)
            if path_match:
                assert line.strip() == "mx0\t1+:0-1981 2+:0-2329"
    os.remove("out_scaf.k32.w1000.n2.path")


def test_mx_rc():
    cmd = "../minimizer_assembler-make assemble -Bd list_files=\'ref.fa scaf.rc.fa\' " \
          "list_weights=\'2 1\' k=32 w=1000 n=2 prefix=out_scaf.rc.k32.w1000.n2"
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    with open("out_scaf.rc.k32.w1000.n2.path", 'r') as paths:
        for line in paths:
            print(line)
            path_match = re.search(r'^mx', line)
            if path_match:
                assert line.strip() == "mx0\t5+:0-1981 2-:0-2329"
    os.remove("out_scaf.k32.w1000.n2.path")
