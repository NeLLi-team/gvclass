import os, sys, shutil, glob
from subprocess import call
from Bio import SeqIO
import os.path
import pandas as pd

### fschulz@lbl.gov, December 2016 ###

fnafile = sys.argv[1]
gffout = sys.argv[2]
genecalling_statsout = sys.argv[3]
summary_statsout = sys.argv[4]
finalfaa = sys.argv[5]

def calc_stats(fnain, faain):
  """
  Calculate coding density, assembly size, GC%
  """
  cumGC = 0
  cumLEN = 0
  cumLENaa = 0
  genecount = 0

  for seq_record in SeqIO.parse(fnain, "fasta"):
    cumGC += seq_record.seq.upper().count('G') + seq_record.seq.upper().count('C')
    cumLEN += len(seq_record.seq)

  for seq_record in SeqIO.parse(faain, "fasta"):
    cumLENaa += len(seq_record.seq)
    genecount += 1
  return [cumLEN, format(cumGC/cumLEN*100, ".2f"), genecount, format(cumLENaa*3/cumLEN*100, ".2f")]


def check_filename(fnafile):
  """
  Input has to be .fna with no additional . in filename
  """
  if len(fnafile.split(".")) > 2:
    print ("fnafile name contains additional '.'")
    sys.exit()
  if fnafile.split(".")[1] != "fna":
    print ("fnafile ending needs to be .fna")
    sys.exit()


def rename_header(bestcode, finalfaa, gffout):
  """
  Rename sequence header to ><filename>|<proteinid>
  """
  records = []
  infilebase = bestcode.split("/")[-1].split(".faa_code")[0]
  for seq_record in SeqIO.parse(bestcode, "fasta"):
    headerbase = seq_record.description.split("|")[0]
    if headerbase != infilebase:
      seq_record.id = infilebase + "|" + str(seq_record.id).split()[0]
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
    runcommmand = str("prodigal -n -q -f gff" + 
                    " -i " + fnafile +
                    " -a " + faaout  + "_code" + str(code) +
                    " -o " + gffout + "_code" + str(code) +
                    " -g " + str(code))
    call (runcommmand, shell = True)
  else: # output for snakemake with -p meta
    runcommmand = str("prodigal -n -q -f gff -p meta" +
                    " -i " + fnafile +
                    " -a " + faaout + "_codemeta" +
                    " -o " + gffout)
    call (runcommmand, shell = True)
    shutil.copy(gffout, gffout + "_codemeta")

# check input
check_filename(fnafile)

# gene calling
# create output for snakemake
codes = [0,1,4,6,11] # pass these from snakemake config instead
for code in codes:
  run_genecalling_codes(fnafile, code, gffout)

stats_dict = {}
faaoutdir = gffout.rsplit("/", 1)[0] + "/"
for faain in glob.glob(faaoutdir + "*.faa_code*"):
  genomeid = faain.split("/")[-1]
  stats_dict[genomeid] = calc_stats(fnafile, faain)

# gene calling stats to select ttable that yields highest coding density
statsdf = pd.DataFrame.from_dict(stats_dict, columns=["LENbp", "GCperc", "genecount", "CODINGperc"], orient="index")
statsdf = statsdf.sort_values("CODINGperc", ascending=False)
statsdf["ttable"] = statsdf.index.map(lambda x : x.split("_")[-1])
statsdf.insert(0, "query", statsdf.index.map(lambda x : x.split(".faa")[0]))
bestcode = faaoutdir + list(statsdf.index)[0] # genome id with code appended that had highest coding density
statsdf.to_csv(genecalling_statsout, sep="\t", index=None)
statsdf.head(1).to_csv(summary_statsout, sep="\t", index=None)
# rename headers and cleanup
rename_header(bestcode, finalfaa, gffout)
for tempout in glob.glob(gffout+"_code*"):
    os.remove(tempout)
    if os.path.isfile(tempout.replace(".gff", ".faa")):
        os.remove(tempout.replace(".gff", ".faa"))
