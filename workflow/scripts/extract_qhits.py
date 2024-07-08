import os
from collections import defaultdict
from typing import Dict, Set
import click
from Bio import SeqIO


def get_qhits(hmmout: str, outmodel: str) -> Dict[str, Set[str]]:
    """
    Parse HMM output file and extract query hits for the specified model.

    Args:
        hmmout (str): Path to the HMM output file.
        outmodel (str): Name of the output model.

    Returns:
        Dict[str, Set[str]]: Dictionary mapping model names to sets of query hits.
    """
    hits_model_dict: Dict[str, Set[str]] = defaultdict(set)

    try:
        with open(hmmout, "r") as infile:
            with click.progressbar(infile, label='Parsing HMM output') as bar:
                for line in bar:
                    if not line.startswith("#") and outmodel in line:
                        fields = line.split()
                        hits_model_dict[fields[3]].add(fields[0])
    except Exception as e:
        print(f"Error processing {hmmout}: {e}")

    return hits_model_dict


def export_hitsfaa(hits_model_dict: Dict[str, Set[str]], queryfaa: str, queryhitsfaa: str) -> None:
    """
    Export hits from the query protein FASTA file to the output FASTA file.

    Args:
        hits_model_dict (Dict[str, Set[str]]): Dictionary mapping model names to sets of query hits.
        queryfaa (str): Path to the query protein FASTA file.
        queryhitsfaa (str): Path to the output file containing hits.
    """
    query_ids = set()

    try:
        with click.progressbar(hits_model_dict.values(), label='Collecting query IDs') as bar1:
            for hits in bar1:
                query_ids |= set(hits)

        with open(queryhitsfaa, "w") as outfile:
            with click.progressbar(SeqIO.parse(queryfaa, "fasta"), label='Writing hits to file') as bar2:
                for seq_record in bar2:
                    if seq_record.description in query_ids:
                        SeqIO.write(seq_record, outfile, "fasta")
    except Exception as e:
        print(f"Error exporting hits to {queryhitsfaa}: {e}")


@click.command()
@click.option('--hmmout', '-h', type=click.Path(exists=True), required=True, help='Path to HMM output file.')
@click.option('--queryfaa', '-q', type=click.Path(exists=True), required=True, help='Path to query protein FASTA file.')
@click.option('--queryhitsfaa', '-o', type=click.Path(), required=True, help='Path to output file containing hits.')
def main(hmmout: str, queryfaa: str, queryhitsfaa: str) -> None:
    """
    Extract query hits from HMM output and export them to a FASTA file.

    Args:
        hmmout (str): Path to the HMM output file.
        queryfaa (str): Path to the query protein FASTA file.
        queryhitsfaa (str): Path to the output file containing hits.
    """
    outmodel = os.path.splitext(os.path.basename(queryhitsfaa))[0]
    hits_model_dict = get_qhits(hmmout, outmodel)
    export_hitsfaa(hits_model_dict, queryfaa, queryhitsfaa)


if __name__ == '__main__':
    main()