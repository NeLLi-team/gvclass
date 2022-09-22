import subprocess
import sys
import pandas as pd

queryfaa = sys.argv[1] # e.g. pathtofiles + "queries/Fadolivirus.faa
modelscombined = sys.argv[2] # e.g. pathtofiles + "models/GVOGuni9.hmm"
hmmout = sys.argv[3] # e.g. pathtofiles + "Fadolivirus_GVOGuni9.out"
target = sys.argv[4] # GVOG9 (evalue) or UNI56 (cut_ga)
countout = sys.argv[5] # summarized hit counts per model from hmmout

def run_cmd(cmd):
    sp = subprocess.Popen(cmd,
        #stderr=subprocess.PIPE,
        #stdout=subprocess.PIPE,
        shell=True)

    std_out, std_err = sp.communicate()
    print('std_err: ', std_err)
    print('std_out: ', std_out)


def get_models (modelsin):
    """
    Generate list of all model names
    Some models might have no hits in hmmout, but they should be displayed in count matrix
    """
    models = []
    with open(modelsin, "r") as f:
        models = [line.split()[-1].rstrip() for line in f if line.startswith("NAME")]
    return models


def get_markercompleteness (models, hmmout, query):
    """
    Get copy numbers for each marker
    """
    # add 0s to include models that are not in hmmout
    count_dict = { x:0 for x in models }
    seen = []
    with open(hmmout, "r") as f:
        lines = [line.rstrip() for line in f if not line.startswith("#")]
        for line in lines:
            if line.split()[0] not in seen:
                count_dict[line.split()[3]] += 1
                seen.append(line.split()[0])
    count_dict = { query : count_dict }
    return count_dict


### hmmsearch ###
cutoff = "-E 1e-10"
if target == "UNI56":
    cutoff = "--cut_ga"

hmmsearch = ["hmmsearch" +
             " --noali" +
             " " + cutoff + " " +
             " --domtblout " + hmmout +
             " --cpu 4 " + modelscombined +
             " " + queryfaa ]

run_cmd(hmmsearch)

### get counts ###

query = queryfaa.split("/")[-1].split(".")[0]

models = get_models (modelscombined)
count_dict = get_markercompleteness (models, hmmout, query)
pd.DataFrame.from_dict(count_dict).T.to_csv(countout, sep="\t")
