import os
import subprocess
from collections import defaultdict
from typing import List, Dict
import click
from Bio import SeqIO


def run_cmd(cmd: str) -> None:
    """
    Run a command using subprocess and print the standard output and error.

    Args:
        cmd (str): The command to run.
    """
    try:
        sp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = sp.communicate()
        print('std_out:', std_out.decode())
        print('std_err:', std_err.decode())
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}\n{e}")


def run_blastp(queryfaa: str, refdb: str, blastpout: str) -> None:
    """
    Run diamond blastp command with specified parameters.

    Args:
        queryfaa (str): Path to the query FASTA file.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
    """
    dblastp = [
        "diamond blastp --threads 8",
        "--quiet --outfmt 6",
        "--more-sensitive --query-cover 30",
        "--subject-cover 30",
        "-e 1e-5",
        f"-d {refdb} -q {queryfaa} -o {blastpout}"
    ]
    run_cmd(" ".join(dblastp))


def parse_blastp(blastpout: str) -> List[str]:
    """
    Parse the diamond blastp output and return the best hits.

    Args:
        blastpout (str): Path to the diamond blastp output file.

    Returns:
        List[str]: List of best hit subject IDs.
    """
    queryids_hits_dict: Dict[str, List[str]] = defaultdict(list)
    seen: List[str] = []

    try:
        num_queries = len(set(line.split()[0] for line in open(blastpout)))
        with open(blastpout, "r") as infile:
            for line in infile:
                queryname, subjectname, pident = line.split()[0], line.split()[1], float(line.split()[2])

                if subjectname not in seen and len(queryids_hits_dict[queryname]) <= 100 / num_queries and pident != 100:
                    queryids_hits_dict[queryname].append(subjectname)
                    seen.append(subjectname)
    except Exception as e:
        print(f"Error parsing BLAST output: {e}")
        return []

    besthits = [x for v in queryids_hits_dict.values() for x in v]
    return besthits


def reduce_ref(queryfaa: str, reffaa: str, reduced_reffaa_merged: str, refdb: str, blastpout: str) -> None:
    """
    Reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
    """
    run_blastp(queryfaa, refdb, blastpout)
    besthits = parse_blastp(blastpout)

    seenids: List[str] = []
    outseqs: List[SeqIO.SeqRecord] = []

    for seq_record in SeqIO.parse(queryfaa, "fasta"):
        if seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

    for seq_record in SeqIO.parse(reffaa, "fasta"):
        if seq_record.description in besthits and seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)

    SeqIO.write(outseqs, reduced_reffaa_merged, "fasta")


@click.command()
@click.option('--queryfaa', '-q', type=click.Path(exists=True), required=True, help='Query sequences from hmmsearch for a single marker')
@click.option('--reffaa', '-r', type=click.Path(exists=True), required=True, help='Reference sequences for a single marker')
@click.option('--refdb', '-d', type=click.Path(exists=True), required=True, help='Diamond formatted database for a single marker')
@click.option('--blastpout', '-b', type=click.Path(), required=True, help='Output file name for diamond blastp results')
@click.option('--reduced_reffaa_merged', '-o', type=click.Path(), required=True, help='Output file name for reduced reference sequences')
def main(queryfaa: str, reffaa: str, refdb: str, blastpout: str, reduced_reffaa_merged: str) -> None:
    """
    Main function to reduce the reference sequences based on the best hits from diamond blastp.

    Args:
        queryfaa (str): Path to the query FASTA file.
        reffaa (str): Path to the reference FASTA file.
        refdb (str): Path to the reference database.
        blastpout (str): Path to the output file for diamond blastp results.
        reduced_reffaa_merged (str): Path to the output file for reduced reference sequences.
    """
    if os.path.getsize(queryfaa) > 0:
        reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout)
    else:
        print(f"{queryfaa} is empty, no hits, omitting blast")


if __name__ == '__main__':
    main()