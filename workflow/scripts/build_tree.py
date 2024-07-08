import os
import glob
import subprocess
from typing import List
import click


def run_cmd(cmd: str) -> None:
    """
    Run a command using subprocess and print the standard output and error.

    Args:
        cmd (str): The command to run.
    """
    try:
        sp = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = sp.communicate()
        if std_err:
            print('std_err:', std_err.decode('utf-8'))
        if std_out:
            print('std_out:', std_out.decode('utf-8'))
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}\n{e}")


def cleanup_iqtree_output(prefix: str) -> None:
    """
    Remove all IQ-TREE generated files except for the .treefile.

    Args:
        prefix (str): The prefix of the IQ-TREE output files.
    """
    for filepath in glob.glob(prefix + "*"):
        if not filepath.endswith(".treefile"):
            os.remove(filepath)
    print("Cleaned up IQ-TREE output files, keeping only .treefile")


def build_tree_fasttree(aln: str, treeout: str) -> None:
    """
    Build a tree using FastTree.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
    """
    cmd = f"FastTree -lg -spr 4 -mlacc 2 -slownni < {aln} > {treeout}"
    run_cmd(cmd)


def build_tree_iqtree(aln: str, treeout: str) -> None:
    """
    Build a tree using IQ-TREE.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
    """
    output_prefix = treeout.rsplit('.', 1)[0]
    cmd = f"iqtree -s {aln} -m LG4X -fast -nt 4 -quiet -pre {output_prefix}"
    run_cmd(cmd)

    # Move the .treefile to the desired output
    os.rename(f"{output_prefix}.treefile", treeout)

    # Clean up IQ-TREE output files except for the .treefile
    cleanup_iqtree_output(output_prefix)


@click.command()
@click.option("--aln", "-a", type=click.Path(exists=True), required=True, help="Alignment file in FASTA format.")
@click.option("--treeout", "-t", type=click.Path(), required=True, help="Output file for the tree.")
@click.option("--method", "-m", type=click.Choice(['fasttree', 'iqtree']), default="fasttree", help="Method for tree building.")
def main(aln: str, treeout: str, method: str) -> None:
    """
    Build a tree using the specified method.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
        method (str): The method for tree building (fasttree or iqtree).
    """
    # Check if alignment file is not empty
    if os.path.getsize(aln) > 0:
        try:
            if method == "fasttree":
                build_tree_fasttree(aln, treeout)
            elif method == "iqtree":
                build_tree_iqtree(aln, treeout)
            
            if os.path.exists(treeout):
                print(f"Tree built successfully with {method} and saved to {treeout}.")
            else:
                print(f"Error: {treeout} not found after running {method}.")
        except subprocess.CalledProcessError as e:
            print(f"Error occurred while running {method}: {e}")
    else:
        print(f"The alignment file {aln} is empty. Skipping tree building.")


if __name__ == "__main__":
    main()