import sys
import subprocess


queryfaamerged = sys.argv[1]
alnout = sys.argv[2]
trimout = sys.argv[3]


def run_cmd(cmd):
    sp = subprocess.Popen(cmd, \
        stderr=subprocess.PIPE, \
        stdout=subprocess.PIPE, \
        shell=True)

    (std_out, std_err) = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


mafftaln = ["mafft --quiet --thread 4 " + queryfaamerged + " > " + alnout]
trimaln = ["trimal -gt 0.1 -in " + alnout + " -out " + trimout]

for cmd in [mafftaln, trimaln]:
	run_cmd(cmd)
