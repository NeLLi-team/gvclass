import click
import pandas as pd
import glob
from collections import Counter
from typing import List, Dict, Tuple, Union

def parse_result(tabout: str) -> List[List[str]]:
    """
    Parse the result file and extract valid lines.

    Args:
        tabout (str): Path to the result file.

    Returns:
        List[List[str]]: List of parsed lines.
    """
    results = []
    with open(tabout) as infile:
        for line in infile:
            if len(line.split()) == 5:
                results.append(line.strip().split())
    return results


def format_floats(df: pd.DataFrame) -> pd.DataFrame:
    # Select columns with float dtype
    float_cols = df.select_dtypes(include=['float64']).columns
    
    # Format these columns to two decimal places
    for col in float_cols:
        df[col] = df[col].map(lambda x: f"{x:.2f}")
    
    return df


def process_model_hitcounts(models_count: str) -> Tuple[int, ...]:
    """
    Process the model hit counts and calculate various metrics.

    Args:
        models_count (str): Path to the models count file.

    Returns:
        Tuple[int, ...]: Tuple containing various calculated metrics.
    """
    try:
        mirus_models = ["mOG0000040", "mOG0000014", "mOG0000019", "mOG0000030", "mOG0000020"]
        busco_models = ["1003258at2759", "1019762at2759", "1032689at2759", "1038775at2759", "1057950at2759", "1079130at2759", "1079827at2759",
                        "1090038at2759", "1106766at2759", "1111142at2759", "1157302at2759", "1161199at2759", "1178688at2759", "1193442at2759",
                        "1194691at2759", "1197019at2759", "1200489at2759", "1220881at2759", "1222562at2759", "1236198at2759", "1249159at2759",
                        "1260807at2759", "1275837at2759", "1284731at2759", "1291729at2759", "1309031at2759", "1322299at2759", "1338131at2759",
                        "1346165at2759", "1346432at2759", "1348942at2759", "1380409at2759", "1382842at2759", "1423847at2759", "1430056at2759",
                        "1432328at2759", "1434266at2759", "1455730at2759", "1488436at2759", "1504863at2759", "1530008at2759", "1552706at2759",
                        "1567796at2759", "1576404at2759", "1588798at2759", "1617752at2759", "1623701at2759", "1645187at2759", "176625at2759",
                        "257318at2759", "296129at2759", "320059at2759", "383503at2759", "426305at2759", "543764at2759", "582756at2759", "625387at2759",
                        "674160at2759", "692986at2759", "734666at2759", "736068at2759", "777920at2759", "834694at2759", "855834at2759", "869548at2759",
                        "923657at2759", "939345at2759", "974865at2759", "976469at2759"]
        phage_models = ["genomad000023", "genomad000026", "genomad000172", "genomad000217", "genomad000254", "genomad000616", "genomad000979", "genomad001596",
                        "genomad001683", "genomad006053", "genomad008072", "genomad016854", "genomad020215", "genomad020681", "genomad029205", "genomad051630",
                        "genomad057391", "genomad066074", "genomad090963", "genomad132751"]
        gvog4m_models = ["GVOGm0461", "GVOGm0022", "GVOGm0023", "GVOGm0054"]
        gvog8m_models = ["GVOGm0013", "GVOGm0022", "GVOGm0023", "GVOGm0054", "GVOGm0172", "GVOGm0461", "GVOGm0760", "GVOGm0890"]
        uni56_models = ["COG0013", "COG0016", "COG0018", "COG0048", "COG0049", "COG0051", "COG0052", "COG0060", "COG0072", "COG0080",
                        "COG0081", "COG0085", "COG0086", "COG0087", "COG0088", "COG0089", "COG0090", "COG0091", "COG0092", "COG0093",
                        "COG0094", "COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0127", "COG0130",
                        "COG0164", "COG0172", "COG0184", "COG0185", "COG0186", "COG0193", "COG0197", "COG0198", "COG0200", "COG0201",
                        "COG0202", "COG0216", "COG0233", "COG0244", "COG0255", "COG0256", "COG0343", "COG0481", "COG0495", "COG0504",
                        "COG0519", "COG0532", "COG0533", "COG0541", "COG0691", "COG0858"]
        mcp_models = ["gamadvirusMCP", "yaravirusMCP", "PoxMCP", "GVOGm0003", "mOG0000014"]
        mrya_models = ["HUH", "HUHlong", "VLTF2", "VLTF3", "ATPase", "gamadvirusMCP", "yaravirusMCP"]

        models_count_df = pd.read_csv(models_count, sep="\t", index_col=0)

        gvog4_unique = models_count_df[gvog4m_models].gt(0).sum(axis=1).sum()
        
        gvog8_unique = models_count_df[gvog8m_models].gt(0).sum(axis=1).sum()
        gvog8_total = models_count_df[gvog8m_models].sum(axis=1).sum()
        gvog8_dup = gvog8_total / gvog8_unique if gvog8_unique > 0 else 0
        
        mirus_unique = models_count_df[mirus_models].gt(0).sum(axis=1).sum()
        mirus_total = models_count_df[mirus_models].sum(axis=1).sum()
        mirus_dup = mirus_total / mirus_unique if mirus_unique > 0 else 0
        
        mrya_unique = models_count_df[mrya_models].gt(0).sum(axis=1).sum()
        mrya_total = models_count_df[mrya_models].sum(axis=1).sum()

        mcp_unique = models_count_df[mcp_models].gt(0).sum(axis=1).sum()
        mcp_total = models_count_df[mcp_models].sum(axis=1).sum()

        uni56_unique = models_count_df[uni56_models].gt(0).sum(axis=1).sum()
        uni56_total = models_count_df[uni56_models].sum(axis=1).sum()
        uni56_dup = uni56_total / uni56_unique if uni56_unique > 0 else 0
            
        busco_unique = models_count_df[busco_models].gt(0).sum(axis=1).sum()
        busco_total =  models_count_df[busco_models].sum(axis=1).sum()
        busco_dup = busco_total / busco_unique if busco_unique > 0 else 0
        
        cellular_unique = uni56_unique + busco_unique
        cellular_total = uni56_total + busco_total
        cellular_dup =(uni56_total + busco_total) / (uni56_unique + busco_unique) if  (uni56_unique + busco_unique) > 0 else 0
            
        phage_unique = models_count_df[phage_models].gt(0).sum(axis=1).sum()
        phage_total = models_count_df[phage_models].sum(axis=1).sum()

        print(f"gvog8 count total: {gvog8_total} and gvog8 count unique: {gvog8_unique}")
        return gvog4_unique, gvog8_unique, gvog8_total, gvog8_dup, mirus_unique, mirus_total, \
            mirus_dup, mrya_unique, mrya_total, mcp_unique, mcp_total, uni56_unique, uni56_total, \
                uni56_dup, busco_unique, busco_total, busco_dup, phage_unique, phage_total, \
                    cellular_unique, cellular_total, cellular_dup
    except Exception as e:
        print(f"Error processing model hit counts: {e}")
        return (0,) * 19

def most_frequent(taxstrings: List[str], taxlevel: str) -> str:
    """
    Find the most frequent taxonomic string at a given taxonomic level.

    Args:
        taxstrings (List[str]): List of taxonomic strings.
        taxlevel (str): Taxonomic level.

    Returns:
        str: Most frequent taxonomic string or "_" if not found.
    """
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
    """
    Determine the taxonomic domain based on the subject and taxannot columns.

    Args:
        row (pd.Series): DataFrame row.

    Returns:
        str: Taxonomic domain.
    """
    subject = row.get("subject", "")
    taxannot = row.get("taxannot", "")
    
    if subject.startswith(("EUK", "ARC", "BAC", "PHAGE", "MIRUS")):
        return subject.split("__")[0]
    elif "Nucleocytoviricota" in taxannot.split("|"):
        return "NCLDV"
    else:
        return "CONFLICT"

def get_final_tax(df: pd.DataFrame, query: str) -> List[Union[str, float]]:
    """
    Get the final taxonomic annotations and counts/percentages for a query.

    Args:
        df (pd.DataFrame): DataFrame containing taxonomic annotations.
        query (str): Query sequence ID.

    Returns:
        List[Union[str, float]]: List of final taxonomic annotations and counts/percentages.
    """
    df_q = df[df.queryname == query]
    tax_levels = ["species", "genus", "family", "order", "class", "phylum", "domain"]
    final_tax = [query]

    def format_counts_percent(lst: List[str]) -> str:
        """
        Format the counts and percentages for each taxonomic level.

        Args:
            lst (List[str]): List of taxonomic annotations.

        Returns:
            str: Formatted counts and percentages.
        """
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

        # Check if the percentage is 100% for the current taxonomic level (taxonomy_strict)
        if len(lst) > 0 and (lst.count(lst[0]) / len(lst)) == 1:
            taxonomy_strict.append(f"{level[0]}_{lst[0].split('__')[-1]}")
        else:
            taxonomy_strict.append(f"{level[0]}_")

        # Check if any taxonomic string has a percentage greater than 50% (taxonomy_majority)
        majority_tax = [tax for tax, count in Counter(lst).items() if count / len(lst) > 0.5]
        if majority_tax:
            taxonomy_majority.append(f"{level[0]}_{majority_tax[0].split('__')[-1]}")
        else:
            taxonomy_majority.append(f"{level[0]}_")

    final_tax.append(str(df_q["distance"].mean()))
    final_tax.append(";".join(taxonomy_strict))  # Add the "taxonomy_strict" column
    final_tax.append(";".join(taxonomy_majority))  # Add the "taxonomy_majority" column
    return final_tax

def tax_annotation(row: pd.Series, level_other: int, level_ncldv: int) -> str:
    """
    Annotate the taxonomic level based on the subject and taxannot columns.

    Args:
        row (pd.Series): DataFrame row.
        level_other (int): Taxonomic level for non-NCLDV subjects.
        level_ncldv (int): Taxonomic level for NCLDV subjects.

    Returns:
        str: Taxonomic annotation.
    """
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
    """
    Process the model hit counts and calculate various metrics.

    Args:
        conversiontable (str): Path to dataframe that contains order level info.
        marker_count (str): Path to dataframe that contains order level marker counts.
        ncldv_order (str): Order-level prediction of the query

    Returns:
        Tuple[float, ...]: Tuple containing various calculated metrics.
    """
    
    order_models = ["OG0", "OG1023", "OG103", "OG1072", "OG107", "OG1082", "OG10", "OG11484", "OG11485", "OG1157", "OG11", "OG1212", "OG1220",
                    "OG1221", "OG1222", "OG1223", "OG1224", "OG1225", "OG1226", "OG1227", "OG1228", "OG1229", "OG1230", "OG1233", "OG1234",
                    "OG1235", "OG1252", "OG1261", "OG1262", "OG1286", "OG1289", "OG12", "OG13007", "OG13008", "OG13009", "OG1314", "OG1316",
                    "OG1322", "OG1341", "OG1347", "OG1350", "OG1382", "OG13", "OG1410", "OG1411", "OG1420", "OG1421", "OG1422", "OG1423", "OG1424",
                    "OG1425", "OG1431", "OG1446", "OG1460", "OG1461", "OG1462", "OG1463", "OG1464", "OG1465", "OG1466", "OG1467", "OG1468", "OG1469",
                    "OG1470", "OG1471", "OG1472", "OG1499", "OG14", "OG1500", "OG1501", "OG1502", "OG1506", "OG152", "OG1547", "OG1548",
                    "OG1549", "OG1550", "OG1551", "OG1554", "OG1555", "OG1569", "OG1599", "OG15", "OG1601", "OG1604", "OG1612", "OG1618", "OG1657",
                    "OG1658", "OG16616", "OG16617", "OG16619", "OG1661", "OG1685", "OG1691", "OG1699", "OG16", "OG1701", "OG1758", "OG178", "OG1793",
                    "OG1796", "OG1797", "OG179", "OG17", "OG1801", "OG1804", "OG1851", "OG1898", "OG1899", "OG18", "OG1938", "OG19487", "OG19488", "OG19489",
                    "OG19490", "OG19491", "OG19492", "OG19493", "OG19494", "OG19496", "OG19497", "OG19498", "OG19499", "OG194", "OG19500", "OG19501", "OG19504",
                    "OG1953", "OG1954", "OG1973", "OG1978", "OG19", "OG1", "OG2024", "OG2045", "OG2056", "OG2083", "OG2090", "OG2096", "OG20", "OG2104", "OG2110",
                    "OG2119", "OG215", "OG2164", "OG2178", "OG2186", "OG21", "OG2209", "OG224", "OG22", "OG23", "OG2437", "OG2473", "OG24", "OG2532", "OG255",
                    "OG258", "OG25", "OG260", "OG26", "OG273", "OG276", "OG2791", "OG27", "OG28", "OG29", "OG2", "OG302", "OG30", "OG31", "OG32", "OG33", "OG345",
                    "OG34", "OG35", "OG36", "OG375", "OG37", "OG383", "OG388", "OG38", "OG396", "OG39", "OG3", "OG40", "OG41", "OG42", "OG430", "OG43", "OG445",
                    "OG448", "OG44", "OG454", "OG45", "OG464", "OG468", "OG46", "OG47", "OG484", "OG48", "OG49", "OG4", "OG502", "OG50", "OG515", "OG51", "OG52",
                    "OG53", "OG540", "OG545", "OG54", "OG553", "OG55", "OG56", "OG57", "OG589", "OG58", "OG597", "OG59", "OG5", "OG609", "OG60", "OG61", "OG621",
                    "OG62", "OG636", "OG63", "OG646", "OG64", "OG651", "OG659", "OG65", "OG66", "OG67", "OG687", "OG68", "OG69", "OG6", "OG701", "OG703", "OG70",
                    "OG718", "OG71", "OG727", "OG728", "OG72", "OG73", "OG74", "OG751", "OG75", "OG76", "OG77", "OG78", "OG793", "OG7", "OG80", "OG82", "OG84", "OG85",
                    "OG86", "OG897", "OG89", "OG8", "OG91", "OG92", "OG93", "OG95", "OG96", "OG97", "OG98", "OG9"]

    def get_relevant_ogs(ncldv_order: str, conversiontable: pd.DataFrame) -> List[str]:
        """
        Extract the relevant orthogroups for a given order.
        
        Args:
            ncldv_order (str): NCLDV order.
            conversiontable (pd.DataFrame): DataFrame with conversion data.
        
        Returns:
            List[str]: List of relevant orthogroups.
        """
        row = conversiontable[conversiontable['Order'] == ncldv_order]
        if not row.empty:
            ogs = row['Orthogroups'].values[0].split(', ')
            return ogs
        print(f"No orthogroups found for {ncldv_order}")
        return []

    def calculate_completeness_and_duplication(df: pd.DataFrame, ogs: List[str]) -> Tuple[float, float]:
        """
        Calculate the completeness and duplication factor for given orthogroups in the DataFrame.
        
        Args:
            df (pd.DataFrame): DataFrame with order counts.
            ogs (List[str]): List of relevant orthogroups.
        
        Returns:
            Tuple[float, float]: Completeness percentage and duplication factor.
        """
        valid_ogs = [og for og in ogs if og in df.columns]
        if not valid_ogs:
            return 0, 0
        
        df_ogs = df[valid_ogs]
        completeness_percentage = (df_ogs > 0).sum().sum() / (df_ogs.shape[0] * len(valid_ogs)) * 100
        total_counts = df_ogs.sum().sum()
        non_zero_ogs_count = (df_ogs > 0).sum(axis=0).sum()  # Corrected here
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
    """
    Summarize the DataFrame by applying taxonomic annotations and counts/percentages.

    Args:
        df (pd.DataFrame): Input DataFrame.

    Returns:
        pd.DataFrame: Summarized DataFrame.
    """
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
    """
    Main function to process and summarize the data.

    Args:
        nn_tree (str): Path to the nearest neighbor tree file.
        marker_count (str): Path to the marker count file.
        querystats (str): Path to the query stats file.
        summary_out (str): Path to the output summary file.
        conversiontable (str): Path to dataframe that contains order level info.
    """
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
            print(df_results_tree["domain"].iloc[0])
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
                #df_results_tree["uni56_unique"] = uni56_unique
                #df_results_tree["uni56_total"] = uni56_total
                #df_results_tree["UNI56_dup"] = uni56_dup
                #df_results_tree["BUSCO_unique"] = busco_unique
                #df_results_tree["BUSCO_total"] = busco_total
                #df_results_tree["BUSCO_dup"] = busco_dup
                df_results_tree["phage_unique"] = phage_unique
                df_results_tree["phage_total"] = phage_total
                df_results_tree["cellular_unique"] = cellular_unique
                df_results_tree["cellular_total"] = cellular_total
                df_results_tree["cellular_dup"] = cellular_dup         
                    # Define the desired column order
                final_order = [
                    "query", "taxonomy_strict", "taxonomy_majority",
                    "species", "genus", "family", "order", "class", "phylum", "domain", "avgdist",
                    "order_dup", "order_completeness", "gvog4_unique", "gvog8_unique", "gvog8_total", "gvog8_dup", 
                    "mcp_total", "mirus_unique","mirus_total", "mirus_dup", "mrya_unique", "mrya_total", 
                    "phage_unique", "phage_total", "cellular_unique", "cellular_total", "cellular_dup"
                ]     
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
                "query", "taxonomy_strict", "taxonomy_majority",
                "species", "genus", "family", "order", "class", "phylum", "domain", "avgdist",
                "order_dup", "order_completeness", "gvog4_unique", "gvog8_unique", "gvog8_total", "gvog8_dup", 
                "mcp_total", "mirus_unique","mirus_total", "mirus_dup", "mrya_unique", "mrya_total", 
                "phage_unique", "phage_total", "cellular_unique", "cellular_total", "cellular_dup"
            ]     
            df_results_tree = pd.DataFrame(allresults, columns=cols)
            final_results = pd.merge(df_results_tree, querystats_df, on='query')
            final_results = format_floats(final_results)
            final_results.to_csv(summary_out, sep="\t", index=False)

    except Exception as e:
        print(f"An error occurred: {e}")
        with open(summary_out, "w") as outfile:
            pass

if __name__ == "__main__":
    main()