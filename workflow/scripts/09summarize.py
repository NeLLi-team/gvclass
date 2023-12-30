import click
import pandas as pd
import glob
from collections import Counter


def parse_result(tabout):
    results = []
    with open(tabout) as infile:
        for line in infile:
            original_line = line  # Store the original line for debugging
            line = line.strip().split()
            if line[0].startswith("GVOG"):
                try:
                    distance = float(line[-1])
                    print(f"Processing: {original_line.strip()}, Distance: {distance}")  # Debug print
                except ValueError:
                    print(f"Skipping line due to ValueError: {original_line.strip()}")  # Debug print for non-float values
                    continue

                if line[0] == "GVOGm0461" and distance < 1.8:
                    print(f"Including {line[0]} with distance {distance} (Special Case)")  # Debug print
                    results.append(line)
                elif line[0] != "GVOGm0461" and distance < 2.5:
                    print(f"Including {line[0]} with distance {distance}")  # Debug print
                    results.append(line)
                else:
                    print(f"Excluding {line[0]} with distance {distance} (Above Threshold)")  # Debug print
    return results

def process_gvog9(gvog9_count):
    try:
        gvog9_models = ["GVOGm0003","GVOGm0013","GVOGm0022","GVOGm0023","GVOGm0054","GVOGm0172","GVOGm0461","GVOGm0760","GVOGm0890"]
        gvog7_models = ["GVOGm0013","GVOGm0023","GVOGm0054","GVOGm0172","GVOGm0461","GVOGm0760","GVOGm0890"]
        gvogsout_df = pd.read_csv(gvog9_count, sep="\t", index_col=0)
        gvogsout_df['GVOG7u'] = (gvogsout_df[gvog7_models] > 0).sum(axis=1)
        gvogsout_df['GVOG7t'] = (gvogsout_df[gvog7_models]).sum(axis=1)
        gvogsout_df['GVOG9u'] = (gvogsout_df[gvog9_models] > 0).sum(axis=1)
        gvogsout_df['GVOG9t'] = (gvogsout_df[gvog9_models]).sum(axis=1)
        gvogsout_df['GVOG7df'] = gvogsout_df['GVOG7t'] / gvogsout_df['GVOG7u'] if gvogsout_df['GVOG7t'].any() else 0
        gvogsout_df['GVOG9df'] = gvogsout_df['GVOG9t'] / gvogsout_df['GVOG9u'] if gvogsout_df['GVOG9t'].any() else 0
        return list(gvogsout_df.GVOG9u)[0], list(gvogsout_df.GVOG9t)[0], list(gvogsout_df.GVOG7u)[0], list(gvogsout_df.GVOG7t)[0], list(gvogsout_df.GVOG7df)[0], list(gvogsout_df.GVOG9df)[0], list(gvogsout_df.GVOGm0003)[0]
    except:
        return 0,0,0,0,0,0


def process_uni56(uni56_count):
    try:
        UNI56out_df = pd.read_csv(uni56_count, sep="\t", index_col=0)
        UNI56out_models = list(UNI56out_df.columns)
        UNI56out_df['UNI56u'] = (UNI56out_df[UNI56out_models] > 0).sum(axis=1)
        UNI56out_df['UNI56t'] = (UNI56out_df[UNI56out_models]).sum(axis=1)
        UNI56out_df['UNI56df'] = UNI56out_df['UNI56t'] / UNI56out_df['UNI56u'] if UNI56out_df['UNI56t'].any() else 0
        return list(UNI56out_df.UNI56u)[0], list(UNI56out_df.UNI56t)[0], list(UNI56out_df.UNI56df)[0]
    except:
        return 0, 0, 0


def most_frequent(taxstrings, taxlevel):
    freq = Counter(taxstrings)
    if len(freq) == 1:
        if len(taxstrings) >= 2:
            return taxstrings[0]
        else:
            return "_"
    most_common = freq.most_common(1)[0]
    if most_common[1] > len(taxstrings) // 2:
        return most_common[0]
    elif taxlevel == "domain":
        return "-".join(sorted(list(set(taxstrings))))
    else:
        return "_"


def tax_domains(row):
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC", "PHAGE"]:
        return row["subject"].split("__")[0]
    elif row["taxannot"].split("|")[5] == "Nucleocytoviricota":
        return "NCLDV"
    else:
        return "CONFLICT"

def get_final_tax(df, query, stringency_s):
    df_q = df[df.queryname == query]
    tax_levels = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    tax_lists = [list(df_q[tax]) for tax in tax_levels]
    final_tax = [query]

    if stringency_s == "majority":
        for level, lst in zip(tax_levels, tax_lists):
            most_freq = most_frequent(lst, level)
            final_tax.append(f"{level[0]}_{most_freq}")
            print(f"Level: {level}, Tax list: {lst}, Most frequent: {most_freq}")
    else:
        stringency = int(stringency_s.replace("gte", ""))
        for level, lst in zip(tax_levels, tax_lists):
            print(f"Level: {level}, List: {lst}, Stringency: {stringency}")
            if len(set(lst)) == 1:
                final_tax.append(f"{level[0]}_{lst[0]}")
                print(f"Single value for {level}: {lst[0]}")
            elif len(lst) >= stringency:
                final_tax.append(f"{level[0]}__")
                print(f"Multiple values for {level}, not meeting stringency: {lst}")
            else:
                final_tax.extend(["missing_markers"] * (7 - len(final_tax)))
                print("Insufficient data for {level}, extending with missing_markers")
                break
        else:
            final_tax.extend(["missing_markers"] * (7 - len(final_tax)))

        # Enhanced domain conflict resolution for stringencies gte1, gte2, gte3
        if "NCLDV" in tax_lists[-1]:
            conflicting_domains = [d for d in ["EUK", "BAC", "PHAGE", "ARC"] if d in tax_lists[-1]]
            if conflicting_domains:
                final_tax[-1] = f"d_NCLDV-{'-'.join(conflicting_domains)}"
                print(f"Conflict resolved to {final_tax[-1]} for {query}")
            else:
                if len(set(tax_lists[-1])) == 1:
                    final_tax[-1] = f"d_{tax_lists[-1][0]}"
                    print(f"Single domain for {query}: {final_tax[-1]}")
                else:
                    final_tax[-1] = "d__"
                    print(f"No conflict, setting domain to d__ for {query}")

    final_tax.append(stringency_s)
    final_tax.append(df_q["distance"].mean())
    print(f"Final tax for query {query}: {final_tax}")
    return final_tax


def tax_annotation(row, level_other, level_ncldv):
    taxannot_parts = row["taxannot"].split("|")
    if row["subject"].split("__")[0] in ["EUK", "ARC", "BAC", "PHAGE"]:
        tax_value = taxannot_parts[level_other] if level_other < len(taxannot_parts) else "_"
        return row["subject"].split("__")[0] + "__" + tax_value
    else:
        tax_value = taxannot_parts[level_ncldv] if level_ncldv < len(taxannot_parts) else "_"
        return tax_value

def summarize(df):
    try:
        print("Initial DataFrame:", df)  # Debug print

        # Ensure all expected columns are present in the DataFrame
        expected_columns = ["query", "subject", "taxannot", "GVOG", "distance"]
        missing_columns = [col for col in expected_columns if col not in df.columns]
        for col in missing_columns:
            df[col] = "_"  # Fill missing columns with placeholder

        df["queryname"] = df["query"].str.split("|", expand=True)[0]
        # Correctly calling tax_annotation with both arguments
        df["species"] = df.apply(lambda row: tax_annotation(row, 6, -1), axis=1)
        df["genus"] = df.apply(lambda row: tax_annotation(row, 5, 1), axis=1)
        df["family"] = df.apply(lambda row: tax_annotation(row, 4, 2), axis=1)
        df["order"] = df.apply(lambda row: tax_annotation(row, 3, 3), axis=1)
        df["class"] = df.apply(lambda row: tax_annotation(row, 2, 4), axis=1)
        df["phylum"] = df.apply(lambda row: tax_annotation(row, 1, 5), axis=1)
        df["domain"] = df.apply(lambda row: tax_domains(row), axis=1)

        df_tophit = df.drop_duplicates(subset=['GVOG'])
        num_hits = len(df_tophit)
        print(f"Number of top hits: {num_hits}")  # Debug print

        stringencies = ["gte1", "gte2", "gte3", "majority"]
        applicable_stringencies = stringencies[:min(num_hits, 3) + 1]
        allresults = [get_final_tax(df_tophit, list(df.queryname)[0], stringency) for stringency in applicable_stringencies]
        df_results = pd.DataFrame(allresults, columns=["query", "species", "genus", "family", "order", "class", "phylum", "domain", "stringency", "avgdist"])
        print(f"Results DataFrame: {df_results}")  # Debug print
        return df_results
    except Exception as e:
        print(f"An error occurred in summarize: {e}")  # Debug print for any exception
        return pd.DataFrame()  # Return an empty DataFrame if an error occurs


@click.command()
@click.option('-n', '--nn_tree', required=True, help='Path to the nearest neighbor tree file.')
@click.option('-g', '--gvog9_count', required=True, help='Path to the GVOG9 count file.')
@click.option('-u', '--uni56_count', required=True, help='Path to the UNI56 count file.')
@click.option('-q', '--querystats', required=True, help='Path to the query stats file.')
@click.option('-c', '--xgb_out', required=True, help='Path to the classifier output file.')
@click.option('-s', '--summary_out', required=True, help='Path to the output summary file.')

def main(nn_tree, gvog9_count, uni56_count, querystats, xgb_out, summary_out):
    try:
        query = summary_out.split("/")[-1].split(".")[0]
        treehits = []
        treehits.extend(parse_result(nn_tree))
        print(f"Number of tree hits: {len(treehits)}")  # Debug print

        if not treehits:
            print("No valid tree hits found.")  # Debug print

        querystats_df = pd.read_csv(querystats, sep="\t")
        with open(xgb_out) as f:
            prediction = [line.split('\t')[1].strip() for line in f if line.strip()]

        if len(treehits) > 0:
            df_tree = pd.DataFrame(treehits, columns=["GVOG", "query", "subject", "taxannot", "distance"])
            df_tree["distance"] = df_tree["distance"].astype(float)
            df_results_tree = summarize(df_tree)
            print(f"Results from summarize: {df_results_tree}")  # Debug print

            GVOG9u, GVOG9t, GVOG7u, GVOG7t, GVOG7df, GVOG9df, MCP = process_gvog9(gvog9_count)
            UNI56u, UNI56t, UNI56df = process_uni56(uni56_count)

            # Merge and save results
            if not df_results_tree.empty:
                df_results_tree["GVOG9u"] = GVOG9u
                df_results_tree["GVOG9t"] = GVOG9t
                df_results_tree["GVOG7u"] = GVOG7u
                df_results_tree["GVOG7t"] = GVOG7t
                df_results_tree["GVOG7df"] = GVOG7df
                df_results_tree["GVOG9df"] = GVOG9df
                df_results_tree["MCP"] = MCP
                df_results_tree["UNI56u"] = UNI56u
                df_results_tree["UNI56t"] = UNI56t
                df_results_tree["UNI56df"] = UNI56df

                if len(prediction) > 0:
                    df_results_tree["xgb"] = prediction[0]
                else:
                    df_results_tree["xgb"] = "no_fna"

                final_results = pd.merge(df_results_tree, querystats_df, on='query')
                print(f"Final merged results: {final_results}")  # Debug print
                final_results.to_csv(summary_out, sep="\t", index=False)
            else:
                print("No results to merge and save.")  # Debug print

        elif len(prediction) > 0:
            # Handle cases where there are predictions but no valid tree hits
            allresults = [[query, "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", prediction[0], "no_hits", "missing_markers", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]]
            cols = ["query", "species", "genus", "family", "order", "class", "phylum", "domain", "xgb", "stringency", "avgdist", "GVOG9u", "GVOG9t", "GVOG7u", "GVOG7t", "GVOG7df", "GVOG9df", "MCP", "UNI56u", "UNI56t", "UNI56df"]
            df_results_tree = pd.DataFrame(allresults, columns=cols)
            final_results = pd.merge(df_results_tree, querystats_df, on='query')
            print(f"Final merged results with predictions: {final_results}")  # Debug print
            final_results.to_csv(summary_out, sep="\t", index=False)

        else:
            # Handle cases where there are no tree hits and no predictions
            allresults = [[query, "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers", "no_fna", "no_hits", "missing_markers", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]]
            cols = ["query", "species", "genus", "family", "order", "class", "phylum", "domain", "xgb", "stringency", "avgdist", "GVOG9u", "GVOG9t", "GVOG7u", "GVOG7t", "GVOG7df", "GVOG9df", "MCP", "UNI56u", "UNI56t", "UNI56df"]
            df_results_tree = pd.DataFrame(allresults, columns=cols)
            final_results = pd.merge(df_results_tree, querystats_df, on='query')
            print(f"Final merged results without hits or predictions: {final_results}")  # Debug print
            final_results.to_csv(summary_out, sep="\t", index=False)

    except Exception as e:
        print(f"An error occurred: {e}")  # Debug print for any exception
        with open(summary_out, "w") as outfile:
            pass

if __name__ == '__main__':
    main()
