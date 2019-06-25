import shlex
import subprocess
import re

def test_mx():
    cmd = "../minimizer_assemble.py -p test_min -n 1 -g 50 -l \"1 1\" -k 32 " \
          "ref.fa.k32.w1000.tsv scaf.fa.k32.w1000.tsv"
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    with open("test_min.path", 'r') as paths:
        for line in paths:
            path_match = re.search(r'^mx', line)
            if path_match:
                assert line.strip() == "mx0\t1+:0-1981 2+:0-2329"


def test_mx_rc():
    cmd = "../minimizer_assemble.py -p test_min_rc -n 1 -g 50 -l \"1 1\" -k 32 " \
          "ref.fa.k32.w1000.tsv scaf.rc.fa.k32.w1000.tsv"
    cmd_shlex = shlex.split(cmd)
    subprocess.call(cmd_shlex)
    with open("test_min_rc.path", 'r') as paths:
        for line in paths:
            path_match = re.search(r'^mx', line)
            if path_match:
                assert line.strip() == "mx0\t1+:0-1981 2-:0-2329"