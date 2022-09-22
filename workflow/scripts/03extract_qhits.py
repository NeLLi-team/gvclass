import os
import sys
from Bio import SeqIO
import shutil
from collections import defaultdict

"""
Extract protein sequences of hits in hmmsearch output
And merge with matching ref proteins 
"""
hmmout = sys.argv[1] # contains hits for any number of models or no hits
queryfaa = sys.argv[2] # query faa file 
queryhitsfaa = sys.argv[3] # extracted hits 

def export_hitsfaa(hits_model_dict, queryfaa, queryhitsfaa):
    for model, hits in hits_model_dict.items():
        if len(hits) > 0:
             SeqIO.write([seq_record for seq_record in SeqIO.parse(queryfaa, "fasta") \
                          if seq_record.description in hits], queryhitsfaa, "fasta")

def get_qhits(hmmout, outmodel):
    # one or more hits per model
    hits_model_dict = defaultdict(list)
    with open(hmmout, "r") as infile:
        for line in infile:
            if not line.startswith("#"):
                queryid = line.split()[0]
                model = line.split()[3]
                if queryid not in hits_model_dict[model] and model == outmodel:
                    hits_model_dict[model].append(queryid)
    return hits_model_dict


outmodel = queryhitsfaa.split("/")[-1].split(".")[0]
hits_model_dict = get_qhits(hmmout, outmodel)
export_hitsfaa(hits_model_dict, queryfaa, queryhitsfaa)
