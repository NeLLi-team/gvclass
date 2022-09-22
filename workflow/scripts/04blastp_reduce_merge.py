import sys
import subprocess
from Bio import SeqIO
from collections import defaultdict
import os

"""
Blastp of query(s) vs refs of respective marker
Extract top 100 hits from refs and merge with query(s) 
"""
queryfaa = sys.argv[1] # query sequences from hmmsearch for a single marker
reffaa = sys.argv[2] # ref seqs for a single marker
refdb = sys.argv[3] # diamond formatted db for a single marker
blastpout = sys.argv[4]
reduced_reffaa_merged = sys.argv[5] # outfile name

def run_cmd(cmd):
    sp = subprocess.Popen(cmd,
        #stderr=subprocess.PIPE,
        #stdout=subprocess.PIPE,
        shell=True)
    std_out, std_err = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


def run_blastp(queryfaa, refdb, blastpout):
    # no evalue set, what is the default?
    dblastp = ["diamond blastp --threads 8 \
                --quiet --outfmt 6 \
                --more-sensitive --query-cover 30 \
                --subject-cover 30 \
                -e 1e-10 \
                -d " + refdb + 
                " -q " + queryfaa + 
                " -o " + blastpout]
    #print ("RUNNING: " + str(dblastp))
    run_cmd(dblastp)

def parse_blastp(blastpout):
    queryids_hits_dict = defaultdict(list)
    seen = []
    # there could be multiple queries
    try:
        num_queries = len(set([line.split()[0] for line in open(blastpout)]))
        with open(blastpout, "r") as infile:
            for line in infile:
                queryname = line.split()[0]
                subjectname = line.split()[1]
                # Yield equal number of subject ids per query
                if subjectname not in seen and len(queryids_hits_dict[queryname]) <= 100/num_queries:
                    queryids_hits_dict[queryname].append(subjectname)
                    seen.append(subjectname)
    except:
        pass
    # single list with all 100 top subject ids
    besthits = [x for v in queryids_hits_dict.values() for x in v]
    return besthits


def reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout):
    run_blastp(queryfaa, refdb, blastpout)
    besthits = parse_blastp(blastpout)
    seenids = []
    outseqs = []
    for seq_record in SeqIO.parse(queryfaa, "fasta"):
        if seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)
    for seq_record in SeqIO.parse(reffaa, "fasta"):
        if seq_record.description in besthits and seq_record.id not in seenids:
            outseqs.append(seq_record)
            seenids.append(seq_record.id)
    SeqIO.write(outseqs, reduced_reffaa_merged, "fasta")

if os.path.getsize(queryfaa) > 0:
    reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout)
else: 
    print (queryfaa + " is empty, no hits, omitting blastp")
