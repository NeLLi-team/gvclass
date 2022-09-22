import sys
from Bio import AlignIO
import glob
import pandas as pd
from collections import Counter


# To do: add query protein id to final output

query = sys.argv[1]
aln_outpath = sys.argv[2]
outname = sys.argv[3]
ncldv_labels = sys.argv[4]

def flattenkeys_to_list(dict1, labels_dict, taxlevel):
    if taxlevel=="order":
        level = 3
    elif taxlevel=="family":
        level = 2
    elif taxlevel=="genus":
        level = 1
    dictlist = [  dict1[model] for model in dict1.keys() ]
    flattened = list({k: v for d in dictlist for k, v in d.items()}.keys())
    return [labels_dict[tax.split("|")[0]][level] for tax in flattened]


def hamming_distance(qseq, rseq):
    """
    Get hamming distance between two aligned sequences
    """
    if len(qseq) != len(rseq):
        raise ValueError("Sequences are of unequal length")
    hd = (sum(aa1 != aa2 for aa1, aa2 in zip(qseq, rseq)))
    seqlen = len(qseq)
    # convert hamming distance to id perc
    idperc = (seqlen-hd)/seqlen*100
    return idperc


def get_nn(aln, query):
    """
    Nearest neighbor in alignment based on hamming distance
    """
    hd_dict_paralogs = {} # include nearest neighbor for all paralogs
    for qseq in AlignIO.read(aln, "fasta"):
        if qseq.id.split("|")[0] == query:
            hd_dict = {}
            for rseq in AlignIO.read(aln, "fasta"):
                if rseq.id.split("|")[0] != query:
                    if hamming_distance(qseq.seq, rseq.seq) > 30:
                        hd_dict[rseq.id] = hamming_distance(qseq.seq, rseq.seq)
            nn = max(hd_dict, key=hd_dict.get) # nearest neighbor
            hd_dict_paralogs[nn + ";;;" + qseq.id] = hd_dict[nn]
    #nnp = max(hd_dict_paralogs, key=hd_dict_paralogs.get) # get only ref that is most similar from all paralogs
    #return {nnp : hd_dict_paralogs[nnp]} # return only ref that is most similar
    return hd_dict_paralogs

# all similarities including paralogs
hd_dict_all = {}
# only single paralog that had highest similarity to reference
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
        taxid = line.split("\t")[0]
        #print (line)
        lineage = line.split("\t")[1].replace("\n","").split("|")
        labels_dict[taxid] = lineage

# this is the summary result, show number of best hits per tax string$
# cutoff: tree dist 2, aln 30 idperc --> needs further evaluation$
"""
for taxlevel in ["order", "family", "genus", "species"]
"""
final_aln_dict = Counter(flattenkeys_to_list(hd_dict_all_best, labels_dict, "order"))
aln_result = "|".join([":".join([str(x) for x in pair]) for pair in list(final_aln_dict.items())])

#with open(outname, "w") as outfile:


with open(outname, "w") as outfile:
    if len(aln_result) > 0:
        outfile.write(query + "\t" + aln_result + "\taln\n")
    else:
        outfile.write(query + "\tno_hits\taln\n")
    for model, value in hd_dict_all.items():
        # print highest similarity first
        value = dict(sorted(value.items(), key=lambda item: item[1], reverse=True))
        outfile.write("\n".join(model +
                        "\t" +
                        k.split(";;;")[1] +
                        "\t" +
                        k.split(";;;")[0] +
                        "\t" +
                        "|".join(labels_dict[k.split("|")[0]]) +
                        "\t" +
                        str(v) for k,v in value.items()) +
                        "\n")
