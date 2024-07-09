import click
import pandas as pd
import glob
from collections import Counter
from typing import List, Tuple, Union

def parse_result(tabout: str) -> List[List[str]]:
    results = []
    with open(tabout) as infile:
        for line in infile:
            if len(line.split()) == 5:
                results.append(line.strip().split())
    return results

def format_floats(df: pd.DataFrame) -> pd.DataFrame:
    float_cols = df.select_dtypes(include=['float64']).columns
    for col in float_cols:
        df[col] = df[col].map(lambda x: f"{x:.2f}")
    return df

def process_model_hitcounts(models_count: str) -> Tuple[int, ...]:
    try:
        mirus_models = ["mOG0000040", "mOG0000014", "mOG0000019", "mOG0000030", "mOG0000020"]
        busco_models = ["1003258at2759", "1019762at2759", "1032689at2759", "1038775at2759", "1057950at2759"]
        phage_models = ["genomad000023", "genomad000026", "genomad000172", "genomad000217", "genomad000254"]
        gvog4m_models = ["GVOGm0461", "GVOGm0022", "GVOGm0023", "GVOGm0054"]
        gvog8m_models = ["GVOGm0013", "GVOGm0022", "GVOGm0023", "GVOGm0054", "GVOGm0172", "GVOGm0461"]
        uni56_models = ["COG0013", "COG0016", "COG0018", "COG0048", "COG0049"]
        mcp_models = ["gamadvirusMCP", "yaravirusMCP", "PoxMCP", "GVOGm0003", "mOG0000014"]
        mrya_models = ["HUH", "HUHlong", "VLTF2", "VLTF3", "ATPase", "gamadvirusMCP"]

        models_count_df = pd.read_csv(models_count, sep="\t", index_col=0)

        def calculate_metrics(model_list):
            unique = models_count_df[model_list].gt(0).sum(axis=1).sum()
            total = models_count_df[model_list].sum(axis=1).sum()
            duplication = total / unique if unique > 0 else 0
            return unique, total, duplication

        gvog4_unique, gvog4_total, gvog4_dup = calculate_metrics(gvog4m_models)
        gvog8_unique, gvog8_total, gvog8_dup = calculate_metrics(gvog8m_models)
        mirus_unique, mirus_total, mirus_dup = calculate_metrics(mirus_models)
        mrya_unique, mrya_total, mrya_dup = calculate_metrics(mrya_models)
        mcp_unique, mcp_total, mcp_dup = calculate_metrics(mcp_models)
        uni56_unique, uni56_total, uni56_dup = calculate_metrics(uni56_models)
        busco_unique, busco_total, busco_dup = calculate_metrics(busco_models)
        phage_unique, phage_total, phage_dup = calculate_metrics(phage_models)

        cellular_unique = uni56_unique + busco_unique
        cellular_total = uni56_total + busco_total
        cellular_dup = cellular_total / cellular_unique if cellular_unique > 0 else 0

        return gvog4_unique, gvog8_unique, gvog8_total, gvog8_dup, mirus_unique, mirus_total, \
               mirus_dup, mrya_unique, mrya_total, mcp_unique, mcp_total, uni56_unique, uni56_total, \
               uni56_dup, busco_unique, busco_total, busco_dup, phage_unique, phage_total, \
               cellular_unique, cellular_total, cellular_dup
    except Exception as e:
        print(f"Error processing model hit counts: {e}")
        return (0,) * 23

def most_frequent(taxstrings: List[str], taxlevel: str) -> str:
    freq = Counter(taxstrings)
    if len(freq) == 1:
        return taxstrings[0] if len(taxstrings) >= 2 else "_"
    most_common = freq.most_common(1)[0]
    if most_common[1] > len(taxstrings) // 2:
        return most_common[0]
    elif taxlevel == "domain":
        return "-".join(sorted(list(set(taxstrings))))
    else:
        return "_"

def tax_domains(row: pd.Series) -> str:
    subject = row.get("subject", "")
    taxannot = row.get("taxannot", "")
    if subject.startswith(("EUK", "ARC", "BAC", "PHAGE", "MIRUS")):
        return subject.split("__")[0]
    elif "Nucleocytoviricota" in taxannot.split("|"):
        return "NCLDV"
    else:
        return "CONFLICT"

def get_final_tax(df: pd.DataFrame, query: str) -> List[Union[str, float]]:
    df_q = df[df.queryname == query]
    tax_levels = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    final_tax = [query]

    def format_counts_percent(lst: List[str]) -> str:
        counter = Counter(lst)
        total = sum(counter.values())
        sorted_counts = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        return ','.join([f"{tax}:{count}({(count / total * 100):.2f}%)" for tax, count in sorted_counts])

    taxonomy_strict = []
    taxonomy_majority = []
    for level in tax_levels:
        lst = list(df_q[level])
        formatted_counts_percent = format_counts_percent(lst)
        final_tax.append(formatted_counts_percent)
        print(f"{query}: {level}, Counts & Percent: {formatted_counts_percent}")

        if len(lst) > 0 and (lst.count(lst[0]) / len(lst)) == 1:
            taxonomy_strict.append(f"{level[0]}_{lst[0].split('__')[-1]}")
        else:
            taxonomy_strict.append(f"{level[0]}_")

        majority_tax = [tax for tax, count in Counter(lst).items() if count / len(lst) > 0.5]
        if majority_tax:
            taxonomy_majority.append(f"{level[0]}_{majority_tax[0].split('__')[-1]}")
        else:
            taxonomy_majority.append(f"{level[0]}_")

    final_tax.append(str(df_q["distance"].mean()))
    final_tax.append(";".join(taxonomy_strict))
    final_tax.append(";".join(taxonomy_majority))
    return final_tax

def tax_annotation(row: pd.Series, level_other: int, level_ncldv: int) -> str:
    try:
        taxannot_parts = row.get("taxannot", "").split("|")
        subject = row.get("subject", "").split("__")
        subject_prefix = subject[0] if subject else "_"
        tax_value = taxannot_parts[level_other] if len(taxannot_parts) > level_other else "IndexOutOfRange"
        return f"{subject_prefix}__{tax_value}"
    except IndexError as e:
        print(f"IndexError in tax_annotation: {e}, row: {row}")
        return "_"

def process_model_hitcounts_order(marker_count: str, conversiontable: str, ncldv_order: str) -> Tuple[float, ...]:
    order_models = ["OG0", "OG1023", "OG103", "OG1072", "OG107"]

    def get_relevant_ogs(ncldv_order: str, conversiontable: pd.DataFrame) -> List[str]:
        row = conversiontable[conversiontable['Order'] == ncldv_order]
        if not row.empty:
            ogs = row['Orthogroups'].values[0].split(', ')
            return ogs
        print(f"No orthogroups found for {ncldv_order}")
        return []

    def calculate_completeness_and_duplication(df: pd.DataFrame, ogs: List[str]) -> Tuple[float, float]:
        valid_ogs = [og for og in ogs if og in df.columns]
        if not valid_ogs:
            return 0, 0
        
        df_ogs = df[valid_ogs]
        completeness_percentage = (df_ogs > 0).sum().sum() / (df_ogs.shape[0] * len(valid_ogs)) * 100
        total_counts = df_ogs.sum().sum()
        non_zero_ogs_count = (df_ogs > 0).sum(axis=0).sum()
        duplication_factor = total_counts / non_zero_ogs_count if non_zero_ogs_count > 0 else 0
        
        return completeness_percentage, duplication_factor
    
    try:
        df_conversiontable = pd.read_csv(conversiontable, sep="\t")
        relevant_ogs = get_relevant_ogs(ncldv_order, df_conversiontable)
        
        if relevant_ogs:
            df_order_counts = pd.read_csv(marker_count, sep="\t")
            df_relevant_ogs = pd.DataFrame(0, columns=relevant_ogs, index=df_order_counts.index)
            df_order_counts = pd.concat([df_order_counts, df_relevant_ogs], axis=1)
            df_order_counts = df_order_counts[relevant_ogs]
            completeness, duplication = calculate_completeness_and_duplication(df_order_counts, relevant_ogs)
        else:
            completeness, duplication = 0, 0
            
        return duplication, completeness
    
    except Exception as e:
        print(f"Error processing order-level completeness and duplication: {e}")
        return 0, 0

def summarize(df: pd.DataFrame) -> pd.DataFrame:
    try:
        if df.empty:
            print("Warning: DataFrame is empty.")
            return pd.DataFrame()

        expected_columns = ["query", "subject", "taxannot", "distance"]
        for col in expected_columns:
            if col not in df.columns:
                print(f"Missing expected column: {col}")
                df[col] = "_"
        
        df["queryname"] = df["query"].str.split("|", expand=True)[0]
        df["domain"] = df.apply(lambda row: tax_domains(row), axis=1)
        df["phylum"] = df.apply(lambda row: tax_annotation(row, 0, -1), axis=1)
        df["class"] = df.apply(lambda row: tax_annotation(row, 1, 1), axis=1)
        df["order"] = df.apply(lambda row: tax_annotation(row, 2, 2), axis=1)
        df["family"] = df.apply(lambda row: tax_annotation(row, 3, 3), axis=1)
        df["genus"] = df.apply(lambda row: tax_annotation(row, 4, 4), axis=1)
        df["species"] = df.apply(lambda row: tax_annotation(row, 5, 5), axis=1)
        df_tophit = df.drop_duplicates(subset=['query'])
        num_hits = len(df_tophit)
        print(f"Number of top hits: {num_hits}")
        allresults = [get_final_tax(df_tophit, list(df.queryname)[0])]
        df_results = pd.DataFrame(allresults, columns=["query", "species", "genus", "family", "order", "class", "phylum", "domain", "avgdist", "taxonomy_strict", "taxonomy_majority"])
        return df_results
    except Exception as e:
        print(f"An error occurred in summarize: {e}")
        return pd.DataFrame()

@click.command()
@click.option('-n', '--nn_tree', required=True, help='Path to the nearest neighbor tree file.')
@click.option('-g', '--marker_count', required=True, help='Path to the marker count file.')
@click.option('-q', '--querystats', required=True, help='Path to the query stats file.')
@click.option('-f', '--conversiontable', required=False, help='Conversion table NCLDV with order completeness factor.')
@click.option('-s', '--summary_out', required=True, help='Path to the output summary file.')
def main(nn_tree: str, marker_count: str, querystats: str, conversiontable: str, summary_out: str) -> None:
    final_order = [
        "query", "taxonomy_majority","taxonomy_strict", 
        "species", "genus", "family", "order", "class", "phylum", "domain", "avgdist",
        "order_dup", "order_completeness", "gvog4_unique", "gvog8_unique", "gvog8_total", "gvog8_dup", 
        "mcp_total", "mirus_unique","mirus_total", "mirus_dup", "mrya_unique", "mrya_total", 
        "phage_unique", "phage_total", "cellular_unique", "cellular_total", "cellular_dup",
        "contigs", "LENbp", "GCperc", "genecount", "CODINGperc", "ttable", "delta"
    ]
    try:
        query = summary_out.split("/")[-1].split(".")[0]
        treehits = parse_result(nn_tree)
        print(f"Number of tree hits: {len(treehits)}")
        if not treehits:
            print("No valid tree hits found.")

        querystats_df = pd.read_csv(querystats, sep="\t")

        if treehits:
            df_tree = pd.DataFrame(treehits, columns=["model", "query", "subject", "taxannot", "distance"])
            df_tree["distance"] = df_tree["distance"].astype(float)
            df_results_tree = summarize(df_tree)

            gvog4_unique, gvog8_unique, gvog8_total, gvog8_dup, mirus_unique, mirus_total, mirus_dup, mrya_unique, \
                mrya_total, mcp_unique, mcp_total, uni56_unique, uni56_total, uni56_dup, busco_unique, busco_total, \
                busco_dup, phage_unique, phage_total, cellular_unique, cellular_total, cellular_dup = process_model_hitcounts(marker_count)

            order_dup = 0
            order_completeness = 0
            if not df_results_tree.empty and str(df_results_tree["domain"].iloc[0]).split(":")[0] == "NCLDV":
                ncldv_order = str(df_results_tree["order"].iloc[0]).split(":")[0].replace("NCLDV__", "")
                order_dup, order_completeness = process_model_hitcounts_order(marker_count, conversiontable, ncldv_order)
                print(f"completeness: {order_completeness}; duplication: {order_dup}")

            if not df_results_tree.empty:
                df_results_tree["order_dup"] = order_dup
                df_results_tree["order_completeness"] = order_completeness
                df_results_tree["gvog4_unique"] = gvog4_unique
                df_results_tree["gvog8_unique"] = gvog8_unique
                df_results_tree["gvog8_total"] = gvog8_total
                df_results_tree["gvog8_dup"] = gvog8_dup
                df_results_tree["mirus_unique"] = mirus_unique
                df_results_tree["mirus_total"] = mirus_total
                df_results_tree["mirus_dup"] = mirus_dup
                df_results_tree["mrya_unique"] = mrya_unique
                df_results_tree["mrya_total"] = mrya_total
                df_results_tree["mcp_unique"] = mcp_unique
                df_results_tree["mcp_total"] = mcp_total
                df_results_tree["phage_unique"] = phage_unique
                df_results_tree["phage_total"] = phage_total
                df_results_tree["cellular_unique"] = cellular_unique
                df_results_tree["cellular_total"] = cellular_total
                df_results_tree["cellular_dup"] = cellular_dup

                # Merge DataFrames and reorder columns
                final_results = pd.merge(df_results_tree, querystats_df, on='query')
                final_results = format_floats(final_results)
                final_results = final_results[final_order]
                final_results.to_csv(summary_out, sep="\t", index=False)
            else:
                print("No results to merge and save.")
        else:
            allresults = [[query, "missing_markers", "missing_markers", "missing_markers", "missing_markers",
                        "missing_markers", "missing_markers", "missing_markers", "missing_markers", "missing_markers",
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

            cols = [
                "query", "taxonomy_majority","taxonomy_strict", 
                "species", "genus", "family", "order", "class", "phylum", "domain", "avgdist",
                "order_dup", "order_completeness", "gvog4_unique", "gvog8_unique", "gvog8_total", "gvog8_dup", 
                "mcp_total", "mirus_unique","mirus_total", "mirus_dup", "mrya_unique", "mrya_total", 
                "phage_unique", "phage_total", "cellular_unique", "cellular_total", "cellular_dup"
            ]

            df_results_tree = pd.DataFrame(allresults, columns=cols)
            final_results = pd.merge(df_results_tree, querystats_df, on='query')
            final_results = format_floats(final_results)
            final_results = final_results[final_order]
            final_results.to_csv(summary_out, sep="\t", index=False)
    except Exception as e:
        print(f"An error occurred: {e}")
        with open(summary_out, "w") as outfile:
            pass

if __name__ == "__main__":
    main()
