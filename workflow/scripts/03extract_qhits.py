import os
from Bio import SeqIO
import shutil
from collections import defaultdict
import click

@click.command()
@click.option('--hmmout', '-h', required=True, help='Path to HMM output file.')
@click.option('--queryfaa', '-q', required=True, help='Path to query protein FASTA file.')
@click.option('--queryhitsfaa', '-o', required=True, help='Path to output file containing hits.')
def main(hmmout, queryfaa, queryhitsfaa):
    outmodel = os.path.splitext(os.path.basename(queryhitsfaa))[0]
    hits_model_dict = get_qhits(hmmout, outmodel)
    export_hitsfaa(hits_model_dict, queryfaa, queryhitsfaa)

def export_hitsfaa(hits_model_dict, queryfaa, queryhitsfaa):
    query_ids = set()
    with click.progressbar(hits_model_dict.values(), label='Exporting hits to file') as bar1:
        for hits in bar1:
            query_ids |= set(hits)
    with open(queryhitsfaa, "w") as outfile:
        with click.progressbar(SeqIO.parse(queryfaa, "fasta"), label='Writing hits to file') as bar2:
            for seq_record in bar2:
                if seq_record.description in query_ids:
                    SeqIO.write(seq_record, outfile, "fasta")

def get_qhits(hmmout, outmodel):
    hits_model_dict = defaultdict(set)
    with open(hmmout, "r") as infile:
        with click.progressbar(infile, label='Parsing HMM output') as bar:
            for line in bar:
                if not line.startswith("#") and outmodel in line:
                    hits_model_dict[line.split()[3]].add(line.split()[0])
    return hits_model_dict

if __name__ == '__main__':
    main()
