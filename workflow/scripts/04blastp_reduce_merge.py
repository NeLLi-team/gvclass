import click
import subprocess
from Bio import SeqIO
from collections import defaultdict
import os


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
    dblastp = ["diamond blastp --threads 8 "
               "--quiet --outfmt 6 "
               "--more-sensitive --query-cover 30 "
               "--subject-cover 30 "
               "-e 1e-10 "
               "-d {} -q {} -o {}".format(refdb, queryfaa, blastpout)]
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
                pident = float(line.split()[2])
                # Yield equal number of subject ids per query
                if subjectname not in seen and len(queryids_hits_dict[queryname]) <= 100/num_queries and pident!=100:
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


@click.command()
@click.option('--queryfaa', '-q', required=True, help='Query sequences from hmmsearch for a single marker')
@click.option('--reffaa', '-r', required=True, help='Ref seqs for a single marker')
@click.option('--refdb', '-d', required=True, help='Diamond formatted db for a single marker')
@click.option('--blastpout', '-b', required=True, help='Output file name for diamond blastp')
@click.option('--reduced_reffaa_merged', '-o', required=True, help='Output file name for reduced reference sequences')
def main(queryfaa, reffaa, refdb, blastpout, reduced_reffaa_merged):
    if os.path.getsize(queryfaa) > 0:
        reduce_ref(queryfaa, reffaa, reduced_reffaa_merged, refdb, blastpout)
    else: 
        print (f"{queryfaa} is empty, no hits, omitting blast")

if __name__ == '__main__':
    main()
