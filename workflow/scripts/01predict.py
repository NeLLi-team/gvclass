import click
import pandas as pd
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.preprocessing import LabelEncoder


def get_label_mapping(mapping_file):
    """Returns a dictionary mapping encoded labels to their original values"""
    with open(mapping_file) as f:
        return {int(line.split('\t')[0]): line.strip().split('\t')[1] for line in f}


def get_best_kmers(kmers_file):
    """Returns a list of the best kmers"""
    with open(kmers_file) as f:
        return [line.strip() for line in f.readlines()]


def count_kmers(kmers, fna):
    """Counts the occurrences of kmers in a given FASTA file"""
    kmer_counts_genome_dict = {x:0 for x in kmers}
    total_seqs = sum(1 for _ in SeqIO.parse(fna, "fasta"))
    with click.progressbar(SeqIO.parse(fna, "fasta"), label='Counting kmers', length=total_seqs) as records:
        for seq_record in records:
            for kmer in kmers:
                kmer_counts_genome_dict[kmer] += Seq(seq_record.seq).count_overlap(kmer)
    return kmer_counts_genome_dict


@click.command()
@click.option('--queryseq', '-q', type=click.Path(exists=True), required=True, help='Path to query sequence')
@click.option('--model', '-m', type=click.Path(exists=True), required=True, help='Path to prebuilt model')
@click.option('--kmers', '-k', type=click.Path(exists=True), required=True, help='Path to kmer file')
@click.option('--labels', '-l', type=click.Path(exists=True), required=True, help='Path to labels')
@click.option('--out', '-o', type=click.Path(exists=False), required=True, help='Path to labels')
def main(queryseq, model, kmers, labels, out):
    # Load XGBoost model
    xgbmodel = pickle.load(open(model, 'rb'))

    # Load label mapping
    inv_mapping = get_label_mapping(labels)

    # Load and preprocess query sequence counts
    kmers_list = get_best_kmers(kmers)
    kmer_count_all_dict = {queryseq.split("/")[-1].split(".f")[0]: count_kmers(kmers_list, queryseq)}
    df_qcounts  =  pd.DataFrame.from_dict(kmer_count_all_dict).T.fillna(0)
    df_qcounts = df_qcounts.div(df_qcounts.sum(axis=1), axis=0)
    df_qcounts = df_qcounts.reset_index().drop(['index'], axis=1)

    # Ensure correct order of features
    cols_when_model_builds = xgbmodel.get_booster().feature_names
    df_qcounts = df_qcounts[cols_when_model_builds]

    # Make prediction and print result
    y_pred = xgbmodel.predict(df_qcounts.head(1))
    prediction = inv_mapping[y_pred[0]]
    with open(out, "w") as outfile:
        outfile.write(f"{queryseq.split('/')[-1].split('.f')[0]}\t{prediction}")


if __name__ == "__main__":
    main()
