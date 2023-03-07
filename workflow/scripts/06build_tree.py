import click
import subprocess


def run_cmd(cmd):
    sp = subprocess.Popen(cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        shell=True)
    (std_out, std_err) = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


@click.command()
@click.option("--aln", "-a", type=click.Path(exists=True), help='Alignment file in FASTA format')
@click.option("--treeout", "-t",  type=click.Path(), help='Output file for the tree')
def main(aln, treeout):
    runtree = f"fasttree -lg < {aln} > {treeout}"
    run_cmd(runtree)


if __name__ == "__main__":
    main()
