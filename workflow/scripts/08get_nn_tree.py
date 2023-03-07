import glob
import click
from ete3 import Tree
from collections import defaultdict
from collections import Counter
import pandas as pd


def flattenkeys_to_list(dict1, labels_dict, taxlevel):
    if taxlevel=="order":
        level = 3
    elif taxlevel=="family":
        level = 2
    elif taxlevel=="genus":
        level = 1
    dictlist = [  dict1[model] for model in dict1.keys() ]
    flattened = list({k: v for d in dictlist for k, v in d.items()}.keys())
    return [labels_dict[tax.split("|")[0]][level] for tax in flattened]


def unpack_family(node, leaves, query, seen):
    """
    get all leaves under a parental node
    """
    children = node.get_children()
    for child in children:
        childname = str(child).replace("--","").replace("/-", "").replace("\-", "").replace("\n","")
        if child.is_leaf() and childname.split("|")[0] != query:
                seen.append(childname)
                leaves.append(childname)
        elif not child.is_leaf() and childname != query and childname not in seen:
            seen.append(childname)
            unpack_family(child, leaves, query, seen) 
        elif child.is_leaf() and childname == query:
            seen.append(childname)
        else:
            pass
    if len(leaves)==0:
        parentnode = node.up
        unpack_family(parentnode, leaves, query, seen) 
    return leaves


@click.command()
@click.option('-q', '--query', required=True, help='Query to search in the tree.')
@click.option('-t', '--tree_outpath', required=True, help='Path to the tree directory.')
@click.option('-o', '--outname', required=True, help='Output file name.')
@click.option('-l', '--ncldv_labels', required=True, help='Path to the NCBI taxonomy labels file.')
def main(query, tree_outpath, outname, ncldv_labels):

    def get_closestrelative(query, leaves, node, tree):
        distance = 10
        closestrelative = "nd"
        for child in leaves:
            if child != query:
                newdist = t.get_distance(child, query, topology_only = False)
                if newdist < distance:
                    distance = newdist
                    closestrelative = child
                else:
                    pass
        return closestrelative, distance


    def get_neighbor(t, modelname, query):
        tree_dict_paralogs = {} # include nearest neighbor for all paralogs
        for node in t.traverse():
            if node.name.split("|")[0] == query:
                queryintree = node.name
                print (queryintree)
                # move up one node in the tree
                # collect all terminal leaves under the parent node
                leaves = []
                seen = []
                leaves = unpack_family(node, leaves, query, seen)
                print (leaves)
                closestrelative, distance = get_closestrelative(queryintree, leaves, node, t)
                print (closestrelative)
                tree_dict_paralogs[closestrelative + ";;;" + queryintree] = distance
        #closestrelativep = min(tree_dict_paralogs, key=tree_dict_paralogs.get) # get only ref that is most similar from all paralogs = shorted distance
        #return {closestrelativep : tree_dict_paralogs[closestrelativep]} # return only ref that is most similar
        return tree_dict_paralogs


    labels_dict = {}
    with open(ncldv_labels, "r") as infile:
        for line in infile:
            taxid = line.split("\t")[0]
            lineage = line.split("\t")[1].replace("\n","").split("|")
            labels_dict[taxid] = lineage

    # all distances including paralogs
    tree_dict_all = {}
    # dict with only single paralog that had shortest distance to reference
    tree_dict_all_best = {}
    for ftree in glob.glob(tree_outpath + "/*.FTWAG"):
        try:
            t = Tree(ftree)
            modelname = ftree.split("/")[-1].split(".")[0]
            print (f"{modelname}")
            tree_dict_all[modelname] = get_neighbor(t, modelname, query)
            print (f"{tree_dict_all}")
            tree_dict_all_best[modelname] = {min(tree_dict_all[modelname], key=tree_dict_all[modelname].get) : tree_dict_all[modelname][min(tree_dict_all[modelname], key=tree_dict_all[modelname].get)]}
        except:
            print (f"Something wrong with {ftree}")
            pass

    # this is the summary result, show number of best hits per tax string
    # cutoff for order: tree dist 2, aln 30 idperc --> needs further evaluation

    print (f"{tree_dict_all_best}@@@@@@@@{query}")
    final_tree_dict = Counter(flattenkeys_to_list(tree_dict_all_best, labels_dict, "order"))
    tree_result = "|".join([":".join([str(x) for x in pair]) for pair in list(final_tree_dict.items())])

    with open(outname, "w") as outfile:
        if len(tree_result) > 0:
            outfile.write(f"{query}\t{tree_result}\ttree\n")
        else:
            outfile.write(f"{query}\tno_hits\ttree\n")
        for model, value in tree_dict_all.items():
            # print lowest distance first
            value = dict(sorted(value.items(), key=lambda item: item[1], reverse=False))
            #outfile.write("\n".join(f"{model}\t{k.split(';;;')[1]}\t{k.split(';;;')[0]}\t{'|'.join(labels_dict[k.split('|')[0]])}\t{v}" for k,v in value.items()) + "\n")
            outfile.write("\n".join([
                f"{model}\t{k.split(';;;')[1]}\t{k.split(';;;')[0]}\t{'|'.join(labels_dict[k.split('|')[0]])}\t{v}"
                for k, v in value.items()]) + "\n")


if __name__ == "__main__":
    main()
