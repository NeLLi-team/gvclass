import click
import os
import sys
import shutil
import glob
from subprocess import call
from Bio import SeqIO
import pandas as pd


def calc_stats(fnain, faain):
    """
    Calculate coding density, assembly size, GC%
    """
    cumulative_gc = 0
    cumulative_len = 0
    cumulative_len_aa = 0
    gene_count = 0
    contigs = 0

    for seq_record in SeqIO.parse(fnain, "fasta"):
        cumulative_gc += seq_record.seq.upper().count('G') + seq_record.seq.upper().count('C')
        cumulative_len += len(seq_record.seq)
        contigs += 1

    for seq_record in SeqIO.parse(faain, "fasta"):
        cumulative_len_aa += len(seq_record.seq)
        gene_count += 1

    return [contigs, cumulative_len, round(cumulative_gc / cumulative_len * 100, 2), gene_count, round(cumulative_len_aa * 3 / cumulative_len * 100, 2)]


def check_filename(fnafile):
    """
    Input has to be .fna with no additional . in filename
    """
    if len(fnafile.split(".")) > 2:
        print("fnafile name contains additional '.'")
        sys.exit()
    if fnafile.split(".")[1] != "fna":
        print("fnafile ending needs to be .fna")
        sys.exit()


def rename_header(bestcode, finalfaa, gffout):
    """
    Rename sequence header to ><filename>|<proteinid>
    """
    records = []
    infilebase = os.path.splitext(os.path.basename(bestcode))[0]

    for seq_record in SeqIO.parse(bestcode, "fasta"):
        headerbase = seq_record.description.split("|")[0]
        if headerbase != infilebase:
            seq_record.id = f"{infilebase}|{str(seq_record.id).split()[0]}"
            seq_record.description = ""
            records.append(seq_record)
        elif headerbase == infilebase:
            seq_record.id = str(seq_record.id).split()[0]
            seq_record.description = ""
            records.append(seq_record)

    SeqIO.write(records, finalfaa, "fasta")
    # keep only best gff
    shutil.copy(gffout + "_code" + bestcode.split("_code")[-1], gffout)


def run_genecalling_codes(fnafile, code, gffout):
    """
    Genecalling with prodigal and defined tt
    """
    faaout = gffout.replace(".gff", ".faa")
    if code > 0:
        run_command = f"prodigal -n -q -f gff -i {fnafile} -a {faaout}_code{code} -o {gffout}_code{code} -g {code}"
        call(run_command, shell=True)
    else:  # output for snakemake with -p meta
        run_command = f"prodigal -n -q -f gff -p meta -i {fnafile} -a {faaout}_codemeta -o {gffout}"
        call(run_command, shell=True)
        shutil.copy(gffout, gffout + "_codemeta")


@click.command()
@click.option('--fnafile', '-f', type=click.Path(exists=True), help='Input FNA file')
@click.option('--gffout', '-g', type=click.Path(), help='Output GFF file')
@click.option('--genecalling_statsout', '-gs', '--genecalling-statsout', type=click.Path(), help='Genecalling statistics file')
@click.option('--summary_statsout', '-ss', type=click.Path(), help='Genome stats outfile')
@click.option('--finalfaa', '-fa', type=click.Path(), help='Output FASTA file')
def main(fnafile, gffout, genecalling_statsout, summary_statsout, finalfaa):
    # check input
    check_filename(fnafile)
    # gene calling
    # create output for snakemake
    codes = [0, 1, 4, 6, 11] # pass these from snakemake config instead
    for code in codes:
        run_genecalling_codes(fnafile, code, gffout)
    stats_dict = {}
    faaoutdir = os.path.dirname(gffout)
    for faain in glob.glob(os.path.join(faaoutdir, "*.faa_code*")):
        genome_id = os.path.basename(faain)
        stats_dict[genome_id] = calc_stats(fnafile, faain)

    # gene calling stats to select ttable that yields highest coding density
    stats_df = pd.DataFrame.from_dict(stats_dict, columns=["contigs", "LENbp", "GCperc", "genecount", "CODINGperc"], orient="index")
    stats_df = stats_df.sort_values("CODINGperc", ascending=False)
    stats_df["ttable"] = stats_df.index.map(lambda x: x.split("_")[-1])
    stats_df.insert(0, "query", stats_df.index.map(lambda x: os.path.splitext(x)[0]))
    bestcode = os.path.join(faaoutdir, stats_df.index[0])  # genome id with code appended that had highest coding density
    stats_df.to_csv(genecalling_statsout, sep="\t", index=None)
    stats_df.head(1).to_csv(summary_statsout, sep="\t", index=None)

    # rename headers and cleanup
    rename_header(bestcode, finalfaa, gffout)
    for tempout in glob.glob(f"{gffout}_code*"):
        os.remove(tempout)
        faaout = tempout.replace(".gff", ".faa")
        if os.path.isfile(faaout):
            os.remove(faaout)


if __name__ == '__main__':
    main()
