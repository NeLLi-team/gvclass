import subprocess
import click


def run_cmd(cmd):
    sp = subprocess.Popen(cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        shell=True)
    (std_out, std_err) = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


@click.command()
@click.option("--queryfaamerged", "-q", required=True, type=click.Path(exists=True), help="Input fasta file containing merged fasta sequences")
@click.option("--alnout", "-a", required=True, type=click.Path(), help="Output fasta file for aligned sequences")
@click.option("--trimout", "-t", required=True, type=click.Path(), help="Output fasta file for trimmed sequences")
def main(queryfaamerged, alnout, trimout):
    mafftaln = f"mafft --quiet --thread 4 {queryfaamerged} > {alnout}"
    trimaln = f"trimal -gt 0.1 -in {alnout} -out {trimout}"
    for cmd in [mafftaln, trimaln]:
        run_cmd(cmd)

if __name__ == "__main__":
    main()
