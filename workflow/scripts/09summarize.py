import sys
import pandas as pd
import glob
from collections import Counter

nn_tree = sys.argv[1]
gvog9_count = sys.argv[2]
uni56_count = sys.argv[3]
querystats = sys.argv[4]
summary_out = sys.argv[5]
summary_out2 = sys.argv[6]

def parse_result(tabout):
    results = []
    with open(tabout) as infile:
        for line in infile:
            line = line.strip()
            if line.startswith("GVOG"):
                # cutoff of 2 for distance to nearest neighbor in tree
                if float(line.split()[-1]) < 2:
                    results.append(line.split("\t"))
                else:
                    print (line.split()[0] + " tree distance to neighbor above threshold, hit removed")
    return results

def process_gvog9(gvog9_count):
    gvog7_models = ["GVOGm0013","GVOGm0023","GVOGm0054","GVOGm0172","GVOGm0461","GVOGm0760","GVOGm0890"]
    gvogsout_df = pd.read_csv(gvog9_count, sep="\t", index_col=0)
    gvogsout_df['GVOG7u'] = (gvogsout_df[gvog7_models] > 0).sum(axis=1)
    gvogsout_df['GVOG7t'] = (gvogsout_df[gvog7_models]).sum(axis=1)
    if list(gvogsout_df.GVOG7u)[0] > 0:
        gvogsout_df['GVOG7df'] = gvogsout_df['GVOG7t'] / gvogsout_df['GVOG7u']
    else:
        gvogsout_df['GVOG7df'] = 0
    return list(gvogsout_df.GVOG7u)[0], list(gvogsout_df.GVOG7t)[0], list(gvogsout_df.GVOG7df)[0], list(gvogsout_df.GVOGm0003)[0]


def process_uni56(uni56_count):
    UNI56out_df = pd.read_csv(uni56_count, sep="\t", index_col=0)
    UNI56out_models = list(UNI56out_df.columns)
    UNI56out_df['UNI56u'] = (UNI56out_df[UNI56out_models] > 0).sum(axis=1)
    UNI56out_df['UNI56t'] = (UNI56out_df[UNI56out_models]).sum(axis=1)
    if list(UNI56out_df.UNI56u)[0] > 0:
        UNI56out_df['UNI56df'] = UNI56out_df['UNI56t'] / UNI56out_df['UNI56u']
    else:
        UNI56out_df['UNI56df'] = 0
    return list(UNI56out_df.UNI56u)[0], list(UNI56out_df.UNI56t)[0], list(UNI56out_df.UNI56df)[0]


def most_frequent(taxstrings, taxlevel):
    frequencies = Counter(taxstrings)
    # dereplicate list to see if there is only a single count
    check = list(set([value for value in frequencies.values()]))
    # at least two different counts, one is majority
    if len(check)>1:
        return frequencies.most_common(1)[0][0]
    # single element in list is also a majority, yield it
    # but at least two markers that yielded that single element!
    elif len(check)==1 and taxlevel and len(list(set(taxstrings))) == 1 and len(taxstrings) >= 2:
        return taxstrings[0]
    # multiple elements with equal counts but domain level
    elif len(check)==1 and len(list(set(taxstrings))) > 1 and taxlevel == "domain":
        return "-".join(sorted(list(set(taxstrings))))
    else:
        return "_"


def tax_species(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[6]
    else:
        return row["taxannot"].split("|")[-1] # species

def tax_genus(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[5]
    else:
        return row["taxannot"].split("|")[1] # genus

def tax_family(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[4]
    else:
        return (row["taxannot"].split("|")[2] + "--" +  row["taxannot"].split("|")[0]).replace("---", "") # family

def tax_order(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[3]
    else:
        return row["taxannot"].split("|")[3] # order

def tax_class(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[2]
    else:
        return row["taxannot"].split("|")[4] # class
    
def tax_phylum(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0] + "__" + row["taxannot"].split("|")[1]
    else:
        return row["taxannot"].split("|")[5] # phylum
    
def tax_domains(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC"]:
        return  row["subject"].split("__")[0]
    elif row["taxannot"].split("|")[5] == "Nucleocytoviricota":
        return "NCLDV"
    else:
        return "CONFLICT"
    
def get_final_tax(df, query, stringency_s):
    # stringency defines the number of tax strings that have to be present
    # these tax strings need to strictly match at the same tax level
    # stringency majority tolerates deviations at tax levels if there is a majority
    df_q = df[ (df.queryname==query) ]
    domainlist = list(df_q["domain"])
    phylumlist= list(df_q["phylum"])
    classlist= list(df_q["class"])
    orderlist= list(df_q["order"])
    familylist= list(df_q["family"])
    genuslist= list(df_q["genus"])
    specieslist= list(df_q["species"])
    finaltax = []
    finaltax.append(query)
    if stringency_s == "majority":
        finaltax.append("s_" + most_frequent(specieslist, "species"))
        finaltax.append("g_" + most_frequent(genuslist, "genus"))
        finaltax.append("f_" + most_frequent(familylist, "family"))
        finaltax.append("o_" + most_frequent(orderlist, "order"))
        finaltax.append("c_" + most_frequent(classlist, "class"))
        finaltax.append("p_" + most_frequent(phylumlist, "phylum"))
        finaltax.append("d_" + most_frequent(domainlist, "domain"))
    else:
        stringency = int(stringency_s.replace("gte", ""))
        if stringency <= len(specieslist):
            if len(set(specieslist)) == 1:
                finaltax.append("s_" + specieslist[0])
            else:
                finaltax.append("s__")
            if len(set(genuslist)) == 1:
                finaltax.append("g_" + genuslist[0])
            else:
                finaltax.append("g__")
            if len(set(familylist)) == 1:
                finaltax.append("f_" + familylist[0])
            else:
                finaltax.append("f__")
            if len(set(orderlist)) == 1:
                finaltax.append("o_" + orderlist[0])
            else:
                finaltax.append("o__")
            if  len(set(classlist)) == 1:
                finaltax.append("c_" + classlist[0])
            else:
                finaltax.append("c__")
            if len(set(phylumlist)) == 1:
                finaltax.append("p_" + phylumlist[0])
            else:
                finaltax.append("p__")
            if len(set(domainlist)) == 1:
                finaltax.append("d_" + domainlist[0])
            elif "EUK" in domainlist and "NCLDV" in domainlist:
                finaltax.append("d_" + "EUK-NCLDV")
            elif "BAC" in domainlist and "NCLDV" in domainlist:
                finaltax.append("d_" + "NCLDV-BAC")
            elif "ARC" in domainlist and "NCLDV" in domainlist:
                finaltax.append("d_" + "NCLDV-ARC")
            elif "ARC" in domainlist and "BAC" in domainlist:
                finaltax.append("d_" + "ARC-BAC")
            elif "ARC" in domainlist and "EUK" in domainlist:
                finaltax.append("d_" + "ARC-EUK")
            elif "BAC" in domainlist and "EUK" in domainlist:
                finaltax.append("d_" + "BAC-EUK")
            else:
                finaltax.append("d__")
        else:
            finaltax.extend(["missing_markers"] * 7)
    finaltax.append(stringency_s)
    avgdist = sum(list(df_q["distance"]))/len((df_q["distance"]))
    finaltax.append(avgdist)
    return finaltax


def summarize(df):
    try:
        df["queryname"] = df["query"].apply(lambda x : x.split("|")[0])
        df["species"] =  df.apply(tax_species, axis=1)
        df["genus"] =  df.apply(tax_genus, axis=1)
        df["family"] =  df.apply(tax_family, axis=1)
        df["order"] =  df.apply(tax_order, axis=1)
        df["class"] =  df.apply(tax_class, axis=1)
        df["phylum"] =  df.apply(tax_phylum, axis=1)
        df["domain"] =  df.apply(tax_domains, axis=1)
        df = df[["queryname","species", "genus","family","order", "class", "phylum", "domain", "GVOG", "distance"]]
        df_tophit = df.drop_duplicates(subset=['GVOG'])
        query = list(df.queryname)[0]
        allresults = []
        allresults.append(get_final_tax(df_tophit, query, "gte1"))
        allresults.append(get_final_tax(df_tophit, query, "gte2"))
        allresults.append(get_final_tax(df_tophit, query, "gte3"))
        allresults.append(get_final_tax(df_tophit, query, "majority"))
        df_results  = pd.DataFrame(allresults, columns=["query", "species", "genus", "family", "order", "class", "phylum", "domain", "stringency", "avgdist"])
        return df_results
    except:
        pass


def main():
    try:
        query = summary_out.split("/")[-1].split(".")[0]
        treehits = []
        treehits.extend(parse_result(nn_tree))
        querystats_df = pd.read_csv(querystats, sep="\t")
        if len(treehits) > 0:
            df_tree = pd.DataFrame(treehits, columns=["GVOG", "query", "subject", "taxannot", "distance"])
            df_tree["distance"] = df_tree["distance"].astype(float)
            # maximum distance in tree to yield nn in the final output, parameter can be further tested and adjusted
            df_results_tree = summarize(df_tree)
            GVOG7u, GVOG7t, GVOG7df, MCP = process_gvog9(gvog9_count)
            UNI56u, UNI56t, UNI56df = process_uni56(uni56_count)
            df_results_tree["GVOG7u"] = GVOG7u
            df_results_tree["GVOG7t"] = GVOG7t
            df_results_tree["GVOG7df"] = GVOG7df
            df_results_tree["MCP"] = MCP
            df_results_tree["UNI56u"] = UNI56u
            df_results_tree["UNI56t"] = UNI56t
            df_results_tree["UNI56df"] = UNI56df
            df_results_tree = pd.merge(df_results_tree, querystats_df, on='query')
            df_results_tree.to_csv(summary_out, sep="\t", index=False)
            df_results_tree.to_csv(summary_out2, sep="\t", index=False)
        else:
            allresults = [[query, "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "no_hits", "avgdist", "0", "0", "0", "0", "0", "0", "0"]]
            cols = ["query", "species", "genus", "family", "order", "class", "phylum", "domain", "stringency", "avgdist", "GVOG9u", "GVOG9t", "GVOG9df", "MCP", "UNI56u", "UNI56t", "UNI56df"]
            df_results_tree = pd.DataFrame(allresults, columns=cols)
            df_results_tree = pd.merge(df_results_tree, querystats_df, on='query')
            df_results_tree.to_csv(summary_out, sep="\t", index=False)
            df_results_tree.to_csv(summary_out2, sep="\t", index=False)
    except:
        print ("something went wrong with " + summary_out.split("/")[-1].split(".")[0] )
        with open(summary_out, "w") as outfile:
            pass

main()
