import glob
from typing import List, Dict, Tuple
import click
from ete3 import Tree
from collections import Counter, defaultdict

def flatten_keys_to_list(dict1: Dict[str, Dict[str, float] | Tuple[str, float]], labels_dict: Dict[str, List[str]], taxlevel: str) -> List[str]:
    """
    Flatten the keys of a nested dictionary and convert them to a list of taxonomic labels.

    Args:
        dict1 (Dict[str, Dict[str, float] | Tuple[str, float]]): Nested dictionary containing model names as keys and dictionaries or tuples of closest relatives as values.
        labels_dict (Dict[str, List[str]]): Dictionary mapping taxonomic IDs to lineage information.
        taxlevel (str): Taxonomic level to extract ('order', 'family', or 'genus').

    Returns:
        List[str]: List of taxonomic labels at the specified level.
    """
    level = {'order': 2, 'family': 3, 'genus': 4}[taxlevel]
    flattened = []
    for value in dict1.values():
        if isinstance(value, dict):
            flattened.extend(k.split("|")[0] for k in value.keys())
        elif isinstance(value, tuple):
            flattened.append(value[0].split("|")[0])
    return [labels_dict[tax][level] for tax in flattened]

def unpack_family(node: Tree, query: str, seen: List[str]) -> List[str]:
    """
    Recursively traverse the tree and collect all leaves under a parental node.

    Args:
        node (Tree): Current node in the tree.
        query (str): Query sequence ID.
        seen (List[str]): List of seen node names to avoid duplicates.

    Returns:
        List[str]: List of leaf names under the parental node.
    """
    leaves = []
    for child in node.get_children():
        childname = str(child).strip("/-\n")
        if child.is_leaf() and childname.split("|")[0] != query:
            if childname not in seen:
                seen.append(childname)
                leaves.append(childname)
        elif not child.is_leaf() and childname != query and childname not in seen:
            seen.append(childname)
            leaves.extend(unpack_family(child, query, seen))
    if not leaves:
        leaves = unpack_family(node.up, query, seen)
    return leaves

def get_closest_relative(query: str, leaves: List[str], tree: Tree) -> Tuple[str, float]:
    """
    Find the closest relative to the query sequence among the given leaves.

    Args:
        query (str): Query sequence ID.
        leaves (List[str]): List of leaf names.
        tree (Tree): Phylogenetic tree.

    Returns:
        Tuple[str, float]: Closest relative and its distance to the query.
    """
    min_dist = float('inf')
    closest_relative = None
    for leaf in leaves:
        if leaf != query:
            dist = tree.get_distance(leaf, query, topology_only=False)
            if dist < min_dist:
                min_dist = dist
                closest_relative = leaf
    return closest_relative, min_dist

def get_neighbors(tree: Tree, query: str) -> Dict[str, float]:
    """
    Find the nearest neighbors for all paralogs of the query sequence in the tree.

    Args:
        tree (Tree): Phylogenetic tree.
        query (str): Query sequence ID.

    Returns:
        Dict[str, float]: Dictionary of nearest neighbors and their distances for each paralog.
    """
    tree_dict_paralogs: Dict[str, float] = {}
    for node in tree.traverse():
        if node.name.split("|")[0] == query:
            queryintree = node.name
            seen: List[str] = []
            leaves = unpack_family(node, query, seen)
            closest_relative, distance = get_closest_relative(queryintree, leaves, tree)
            tree_dict_paralogs[f"{closest_relative};;;{queryintree}"] = distance
    return tree_dict_paralogs

@click.command()
@click.option('-q', '--query', required=True, help='Query to search in the tree.')
@click.option('-t', '--tree_outpath', required=True, help='Path to the tree directory.')
@click.option('-o', '--outname', required=True, help='Output file name.')
@click.option('-l', '--ncldv_labels', required=True, help='Path to the NCBI taxonomy labels file.')
def main(query: str, tree_outpath: str, outname: str, ncldv_labels: str) -> None:
    """
    Main function to find the nearest neighbors of a query sequence in phylogenetic trees.

    Args:
        query (str): Query sequence ID.
        tree_outpath (str): Path to the directory containing tree files.
        outname (str): Output file name.
        ncldv_labels (str): Path to the file containing NCBI taxonomy labels.
    """
    labels_dict: Dict[str, List[str]] = {}
    try:
        with open(ncldv_labels, "r") as infile:
            for line in infile:
                taxid, lineage = line.strip().split("\t")
                labels_dict[taxid] = lineage.split("|")
    except Exception as e:
        print(f"Error reading {ncldv_labels}: {e}")

    tree_dict_all: Dict[str, Dict[str, float]] = {}
    tree_dict_all_best: Dict[str, Tuple[str, float]] = {}
    for ftree in glob.glob(f"{tree_outpath}/*.treefile"):
        try:
            tree = Tree(ftree)
            modelname = ftree.split("/")[-1].split(".")[0]
            tree_dict_all[modelname] = get_neighbors(tree, query)
            tree_dict_all_best[modelname] = min(tree_dict_all[modelname].items(), key=lambda x: x[1])
        except Exception as e:
            print(f"Error processing {ftree}: {e}")

    final_tree_dict = Counter(flatten_keys_to_list(tree_dict_all_best, labels_dict, "order"))
    tree_result = "|".join([f"{k}:{v}" for k, v in final_tree_dict.items()])

    try:
        with open(outname, "w") as outfile:
            outfile.write(f"{query}\t{tree_result if tree_result else 'no_hits'}\ttree\n")
            for model, value in tree_dict_all.items():
                value = dict(sorted(value.items(), key=lambda x: x[1]))
                for k, v in value.items():
                    closest_rel, query_seq = k.split(";;;")
                    outfile.write(f"{model}\t{query_seq}\t{closest_rel}\t{'|'.join(labels_dict[closest_rel.split('|')[0]])}\t{v}\n")
    except Exception as e:
        print(f"Error writing to {outname}: {e}")

if __name__ == "__main__":
    main()