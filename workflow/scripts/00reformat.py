from Bio import SeqIO
import pandas as pd
import sys
import os

faain = sys.argv[1]
faaout = sys.argv[2]
faastats = sys.argv[3]

# header in format ><filebasename|<proteinid>

reformatted = []
protcount = 0

for seq_record in SeqIO.parse(faain, "fasta"):
    querybase = faain.split("/")[-1].split(".")[0]
    protcount += 1
    if "|" in seq_record.id:
        if seq_record.id.split("|")[0] == querybase:
            seq_record.id = seq_record.id.split()[0]
            seq_record.description = ""
            reformatted.append(seq_record)
        else:
            seq_record.id= querybase + "|" + seq_record.id.split()[0]
            seq_record.description = ""
            reformatted.append(seq_record)
    else:
        seq_record.id= querybase + "|" + seq_record.id.split()[0]
        seq_record.description = ""
        reformatted.append(seq_record)

if not os.path.isfile(faastats):
    # stats if only faa file provided, only contains gene count
    cols = ["query", "LENbp", "GCperc", "genecount", "CODINGperc", "ttable"]
    faastats_list = [[querybase, "no_fna", "no_fna", protcount, "no_fna", "no_fna"]]
    pd.DataFrame(faastats_list, columns=cols).to_csv(faastats, sep="\t", index=False)

SeqIO.write(reformatted, faaout, "fasta")
