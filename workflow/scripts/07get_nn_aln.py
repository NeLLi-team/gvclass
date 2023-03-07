import glob
from collections import Counter
from Bio import AlignIO
import click

def get_nn(aln, query):
    hd_dict_paralogs = {}
    for qseq in AlignIO.read(aln, "fasta"):
        if qseq.id.split("|")[0] == query:
            hd_dict = {}
            for rseq in AlignIO.read(aln, "fasta"):
                if rseq.id.split("|")[0] != query:
                    hd = hamming_distance(qseq.seq, rseq.seq)
                    if hd > 30:
                        hd_dict[rseq.id] = hd
            nn = max(hd_dict, key=hd_dict.get)
            hd_dict_paralogs[nn + ";;;" + qseq.id] = hd_dict[nn]
    return hd_dict_paralogs


def flattenkeys(dict1, labels_dict, taxlevel):
    if taxlevel == "order":
        level = 3
    elif taxlevel == "family":
        level = 2
    elif taxlevel == "genus":
        level = 1
    dictlist = [dict1[model] for model in dict1]
    flattened = list({k: v for d in dictlist for k, v in d.items()}.keys())
    return [labels_dict[tax.split("|")[0]][level] for tax in flattened]


def hamming_distance(qseq, rseq):
    """
    get hamming distance between two aligned sequences
    """
    # return id%
    hd = (sum(aa1 != aa2 for aa1, aa2 in zip(qseq, rseq)))
    return hd

@click.command()
@click.option("--query", "-q", required=True)
@click.option("--aln_outpath", "-a", required=True)
@click.option("--outname", "-o", required=True)
@click.option("--ncldv_labels", "-l", required=True)
def main(query, aln_outpath, outname, ncldv_labels):
    hd_dict_all = {}
    hd_dict_all_best = {}
    for aln in glob.glob(aln_outpath + "*.mafft01"):
        modelname = aln.split("/")[-1].split(".")[0]
        try:
            hd_dict_all[modelname] = get_nn(aln, query)
            hd_dict_all_best[modelname] = {max(hd_dict_all[modelname], key=hd_dict_all[modelname].get) : hd_dict_all[modelname][max(hd_dict_all[modelname], key=hd_dict_all[modelname].get)]}
        except:
            pass

    labels_dict = {}
    with open(ncldv_labels, "r") as infile:
        for line in infile:
            taxid, lineage = line.strip().split("\t")
            labels_dict[taxid] = lineage.split("|")

    final_aln_dict = Counter(flattenkeys(hd_dict_all_best, labels_dict, "order"))
    aln_result = "|".join([f"{x}:{y}" for x, y in final_aln_dict.items()])

    with open(outname, "w") as outfile:
        if len(aln_result) > 0:
            outfile.write(f"{query}\t{aln_result}\taln\n")
        else:
            outfile.write(f"{query}\tno_hits\taln\n")
        for model, value in hd_dict_all.items():
            value = dict(sorted(value.items(), key=lambda item: item[1], reverse=True))
            outfile.write("\n".join([f"{model}\t{k.split(';;;')[1]}\t{k.split(';;;')[0]}\t{'|'.join(labels_dict[k.split('|')[0]])}\t{v}" for k, v in value.items()]) + "\n")

if __name__ == "__main__":
    main()
