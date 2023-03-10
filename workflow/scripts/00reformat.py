import click
from Bio import SeqIO
import pandas as pd
import os


@click.command()
@click.option('--input', '-i', type=click.Path(exists=True), help='Input FASTA file')
@click.option('--output', '-o', type=click.Path(), help='Output FASTA file')
@click.option('--stats', '-s', type=click.Path(), help='Genome stats outfile')
def main(input, output, stats):
    # read in the input file
    records = list(SeqIO.parse(input, "fasta"))

    # modify the sequence IDs and descriptions
    protcount = 0
    for record in records:
        protcount += 1
        if "|" in record.id:
            record.id = record.id.split("|")[-1]
        record.id = f"{os.path.splitext(os.path.basename(input))[0]}|{record.id}"
        record.description = ""

    # write the reformatted sequences to the output file
    SeqIO.write(records, output, "fasta")

    # calculate statistics and write them to the stats file
    stats_dict = {
        'query': os.path.splitext(os.path.basename(input))[0],
        'contigs': 'no_fna',
        'LENbp': 'no_fna',
        'GCperc': 'no_fna',
        'genecount': protcount,
        'CODINGperc': 'no_fna',
        'ttable': 'no_fna'
    }
    pd.DataFrame.from_dict(stats_dict, orient='index').T.to_csv(stats, sep="\t", index=False)

if __name__ == '__main__':
    main()
