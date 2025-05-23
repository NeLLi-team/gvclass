import os
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
import click

def reformat_sequences(records: List[SeqRecord], input_file: str) -> List[SeqRecord]:
    """
    Reformat the sequence IDs and descriptions.

    Args:
        records (List[SeqRecord]): List of sequence records.
        input_file (str): Path to the input FASTA file.

    Returns:
        List[SeqRecord]: List of reformatted sequence records.
    """
    input_basename = os.path.splitext(os.path.basename(input_file))[0]
    for record in records:
        if "|" in record.id:
            record.id = f"{input_basename}|{record.id.split('|')[-1].replace(':', '__')}"
        else:
            record.id = f"{input_basename}|{record.id.replace(':', '__')}"
        record.description = ""
    return records

def calculate_stats(records: List[SeqRecord], input_file: str) -> pd.DataFrame:
    """
    Calculate statistics for the input FASTA file.

    Args:
        records (List[SeqRecord]): List of sequence records.
        input_file (str): Path to the input FASTA file.

    Returns:
        pd.DataFrame: DataFrame containing the calculated statistics.
    """
    total_length = sum(len(record.seq) for record in records)
    gc_count = sum(record.seq.count('G') + record.seq.count('C') for record in records)
    gc_percentage = (gc_count / total_length * 100) if total_length > 0 else 0
    stats_dict = {
        'query': os.path.splitext(os.path.basename(input_file))[0],
        'contigs': len(records),
        'LENbp': total_length,
        'GCperc': gc_percentage,
        'genecount': len(records),
        'CODINGperc': 'no_fna',
        'ttable': 'no_fna'
    }
    return pd.DataFrame([stats_dict])

@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Input FASTA file', required=True)
@click.option('--output', '-o', type=click.Path(), help='Output FASTA file', required=True)
@click.option('--stats', '-s', type=click.Path(), help='Genome stats outfile', required=True)
def main(input: str, output: str, stats: str) -> None:
    """
    Reformat sequences and calculate statistics for the input FASTA file.

    Args:
        input (str): Path to the input FASTA file.
        output (str): Path to the output FASTA file.
        stats (str): Path to the genome stats output file.
    """
    records = list(SeqIO.parse(input, "fasta"))
    reformatted_records = reformat_sequences(records, input)
    SeqIO.write(reformatted_records, output, "fasta")

    stats_df = calculate_stats(records, input)
    stats_df.to_csv(stats, sep="\t", index=False)

if __name__ == '__main__':
    main()