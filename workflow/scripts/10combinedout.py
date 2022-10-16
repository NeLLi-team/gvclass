import sys
import os.path
import glob
import pandas as pd

resultsdirparent = sys.argv[1]
combinedout = sys.argv[2]


def get_combined(resultsdirparent, combinedout):
    # retain only annotation for d_NCLDV at highest available stringeny
    results = []
    for sumout in glob.glob(resultsdirparent + "*/*.summary.tab"):
        try:
            df = pd.read_csv(sumout, index_col=None, header=0, sep="\t")
            results.append(df)
        except:
            print (sumout + " is empty")
    df_results_f = pd.concat(results, axis=0, ignore_index=True)
    df_results_f = df_results_f.fillna(0)
    df_results_f_gt3 = df_results_f[ (df_results_f["domain"] == "d_NCLDV" )
            &  (df_results_f["stringency"] == "gte3" ) ]
    df_results_f_gt2 = df_results_f[ (df_results_f["domain"] == "d_NCLDV" )
            &  (df_results_f["stringency"] == "gte2" ) ]
    df_results_f_gt1 = df_results_f[ (df_results_f["domain"] == "d_NCLDV" )
            &  (df_results_f["stringency"] == "gte1" ) ]
    df_results_f_combined = pd.concat([df_results_f_gt3,df_results_f_gt2,df_results_f_gt1])
    df_results_f_combined = df_results_f_combined.drop_duplicates(subset=['query'])
    df_results_f_combined.to_csv(combinedout, sep="\t", index=None)

try:
    if not os.path.isfile(combinedout):
        get_combined(resultsdirparent, combinedout)
except:
    print ("No giant virus found")
