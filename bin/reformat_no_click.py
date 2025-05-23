#!/usr/bin/env python
import os
import argparse
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

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

def main():
    """
    Reformat sequences and calculate statistics for the input FASTA file.
    """
    parser = argparse.ArgumentParser(description='Reformat sequences and calculate statistics')
    parser.add_argument('-i', '--input', type=str, help='Input FASTA file', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output FASTA file', required=True)
    parser.add_argument('-s', '--stats', type=str, help='Genome stats outfile', required=True)
    
    args = parser.parse_args()
    
    records = list(SeqIO.parse(args.input, "fasta"))
    reformatted_records = reformat_sequences(records, args.input)
    SeqIO.write(reformatted_records, args.output, "fasta")

    stats_df = calculate_stats(records, args.input)
    stats_df.to_csv(args.stats, sep="\t", index=False)

if __name__ == '__main__':
    main()
