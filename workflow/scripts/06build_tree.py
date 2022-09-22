import sys
import subprocess

aln = sys.argv[1]
treeout = sys.argv[2]

def run_cmd(cmd):
    sp = subprocess.Popen(cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        shell=True)
    (std_out, std_err) = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


runtree = ["fasttree -lg < " +
           aln + ">" + treeout]
run_cmd(runtree)
