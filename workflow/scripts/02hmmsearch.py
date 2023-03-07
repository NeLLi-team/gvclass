import click
import subprocess
import sys
import pandas as pd


def run_cmd(cmd):
    """Runs a command in the shell and prints standard output and error"""
    sp = subprocess.Popen(cmd, shell=True)
    std_out, std_err = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


def get_models(modelsin):
    """
    Generate list of all model names
    Some models might have no hits in hmmout, but they should be displayed in count matrix
    """
    models = []
    with open(modelsin, "r") as f:
        models = [line.split()[-1].rstrip() for line in f if line.startswith("NAME")]
    return models


def get_markercompleteness(models, hmmout, query):
    """Returns a dictionary of marker count for each model"""
    count_dict = {x: 0 for x in models}
    seen = []
    with open(hmmout, "r") as f:
        lines = [line.rstrip() for line in f if not line.startswith("#")]
        for line in lines:
            if line.split()[0] not in seen:
                count_dict[line.split()[3]] += 1
                seen.append(line.split()[0])
    count_dict = {query: count_dict}
    return count_dict


@click.command()
@click.option('--queryfaa', '-q', type=click.Path(exists=True), help='Input query FASTA file')
@click.option('--modelscombined', '-m', type=click.Path(exists=True), help='Combined HMM model file')
@click.option('--hmmout', '-h', type=click.Path(), help='Output file for HMM hits')
@click.option('--target', '-t', type=str, help='Target database (UNI56 or UNIREF90)')
@click.option('--countout', '-c', type=click.Path(), help='Output file for marker gene counts')
def main(queryfaa, modelscombined, hmmout, target, countout):
    ### hmmsearch ###
    cutoff = "-E 1e-10"
    if target == "UNI56":
        cutoff = "--cut_ga"
    hmmsearch = ["hmmsearch" +
                 " --noali" +
                 " " + cutoff + " " +
                 " --domtblout " + hmmout +
                 " --cpu 4 " + modelscombined +
                 " " + queryfaa]
    run_cmd(hmmsearch)

    ### get counts ###
    query = queryfaa.split("/")[-1].split(".")[0]
    models = get_models(modelscombined)
    count_dict = get_markercompleteness(models, hmmout, query)
    pd.DataFrame.from_dict(count_dict).T.to_csv(countout, sep="\t")


if __name__ == '__main__':
    main()
