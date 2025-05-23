import subprocess
from typing import List

import click


def run_cmd(cmd: List[str]) -> None:
    """
    Run a command using subprocess and print the standard output and error.

    Args:
        cmd (List[str]): The command to run as a list of strings.
    """
    sp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    std_out, std_err = sp.communicate()
    print('std_err:', std_err.decode())
    print('std_out:', std_out.decode())


@click.command()
@click.option("--queryfaamerged", "-q", type=click.Path(exists=True), help="Input fasta file containing merged fasta sequences")
@click.option("--alnout", "-a", type=click.Path(), help="Output fasta file for aligned sequences")
@click.option("--trimout", "-t", type=click.Path(), help="Output fasta file for trimmed sequences")
@click.option("--mafftoption", "-o", default="linsi", help="MAFFT alignment option (e.g., 'auto', 'linsi', 'fftns')")
def main(queryfaamerged, alnout, trimout, mafftoption):
    # Count the number of sequences in the input file
    with open(queryfaamerged, "r") as f:
        sequence_count = sum(1 for line in f if line.startswith(">"))

    # Proceed only if there are at least 3 sequences
    if sequence_count >= 3:
        # Adjust the MAFFT command based on the specified option
        if mafftoption == "linsi":
            mafft_command = "mafft-linsi"
        elif mafftoption == "ginsi":
            mafft_command = "mafft-ginsi"
        else:
            # Use the generic mafft command if no specific option is provided
            mafft_command = "mafft --auto"

        # MAFFT alignment command
        mafftaln = f"{mafft_command} --quiet --thread 4 {queryfaamerged} > {alnout}"

        # Trimming of the alignment
        trimaln = f"trimal -gt 0.1 -in {alnout} -out {trimout}"

        for cmd in [mafftaln, trimaln]:
            run_cmd(cmd)
    else:
        print("Alignment not created as there is only one sequence in the input file.")

if __name__ == "__main__":
    main()