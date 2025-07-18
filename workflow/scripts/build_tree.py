import os
import glob
import subprocess
import click

# Import our IQ-TREE Python wrapper
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils import run_command_safe_legacy, validate_file_path, GVClassError

try:
    from src.core.iqtree_python import IQTree
    IQTREE_PYTHON_AVAILABLE = True
except ImportError:
    IQTREE_PYTHON_AVAILABLE = False



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
    Build a tree using VeryFastTree (Python implementation of FastTree).

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
    """
    aln_file = validate_file_path(aln, must_exist=True)
    out_file = validate_file_path(treeout, must_exist=False)
    
    try:
        # First try to use VeryFastTree Python package
        # VeryFastTree can be run as a module or through its API
        # Using command-line style through Python
        import subprocess
        
        cmd = [sys.executable, "-m", "veryfasttree", "-lg", "-spr", "4", "-mlacc", "2", "-slownni", str(aln_file)]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        with open(out_file, 'w') as output_file:
            output_file.write(result.stdout)
            
        if result.stderr:
            print(f'VeryFastTree stderr: {result.stderr}')
            
    except (ImportError, subprocess.CalledProcessError) as e:
        print(f"VeryFastTree failed: {e}")
        print("Falling back to FastTree command")
        
        # Fallback to original FastTree if available
        cmd = ["FastTree", "-lg", "-spr", "4", "-mlacc", "2", "-slownni"]
        
        with open(aln_file, 'r') as input_file:
            result = subprocess.run(
                cmd,
                stdin=input_file,
                capture_output=True,
                text=True,
                check=True
            )
            
            with open(out_file, 'w') as output_file:
                output_file.write(result.stdout)
                
            if result.stderr:
                print(f'FastTree stderr: {result.stderr}')



def build_tree_python_iqtree(aln: str, treeout: str, threads: int = 4) -> None:
    """
    Build a tree using our Python IQ-TREE wrapper.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
        threads (int): Number of threads to use.
    """
    if not IQTREE_PYTHON_AVAILABLE:
        print(f"Python IQ-TREE wrapper not available (IQTREE_PYTHON_AVAILABLE={IQTREE_PYTHON_AVAILABLE}), falling back to command-line")
        return build_tree_iqtree(aln, treeout, threads)
    
    aln_file = validate_file_path(aln, must_exist=True)
    out_file = validate_file_path(treeout, must_exist=False)
    
    output_prefix = str(out_file).rsplit('.', 1)[0]
    
    try:
        # Use our Python wrapper
        iqtree = IQTree()
        outputs = iqtree.run(
            alignment=str(aln_file),
            model="LG4X",
            threads=threads,
            fast=True,
            prefix=output_prefix,
            quiet=True
        )
        
        # Move the .treefile to the desired output
        treefile_path = outputs["tree"]
        if os.path.exists(treefile_path):
            import shutil
            shutil.move(treefile_path, str(out_file))
        else:
            raise GVClassError(f"IQ-TREE did not produce expected treefile: {treefile_path}")
        
        # Clean up IQ-TREE output files
        iqtree.cleanup(output_prefix)
        
    except Exception as e:
        print(f"Error running Python IQ-TREE wrapper: {e}")
        print("Falling back to command-line IQ-TREE")
        return build_tree_iqtree(aln, treeout, threads)


def build_tree_iqtree(aln: str, treeout: str, threads: int = 4) -> None:
    """
    Build a tree using IQ-TREE.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
        threads (int): Number of threads to use.
    """
    aln_file = validate_file_path(aln, must_exist=True)
    out_file = validate_file_path(treeout, must_exist=False)
    
    output_prefix = str(out_file).rsplit('.', 1)[0]
    
    # Check if we're on NixOS - for now, just use iqtree directly
    # The conda environment should handle library paths
    iqtree_cmd = "iqtree"
    
    cmd = [
        iqtree_cmd,
        "-s", str(aln_file),
        "-m", "LG4X",
        "-fast",
        "-nt", str(threads),
        "-quiet",
        "-pre", output_prefix
    ]
    
    run_command_safe_legacy(cmd)

    # Move the .treefile to the desired output
    treefile_path = f"{output_prefix}.treefile"
    if os.path.exists(treefile_path):
        os.rename(treefile_path, str(out_file))
    else:
        raise GVClassError(f"IQ-TREE did not produce expected treefile: {treefile_path}")

    # Clean up IQ-TREE output files except for the .treefile
    cleanup_iqtree_output(output_prefix)


@click.command()
@click.option("--aln", "-a", type=click.Path(exists=True), required=True, help="Alignment file in FASTA format.")
@click.option("--treeout", "-t", type=click.Path(), required=True, help="Output file for the tree.")
@click.option("--method", "-m", type=click.Choice(['fasttree', 'iqtree']), default="fasttree", help="Method for tree building.")
@click.option("--threads", "-j", default=4, help="Number of threads to use for IQ-TREE.")
def main(aln: str, treeout: str, method: str, threads: int) -> None:
    """
    Build a tree using the specified method.

    Args:
        aln (str): The path to the alignment file in FASTA format.
        treeout (str): The path to the output tree file.
        method (str): The method for tree building (fasttree or iqtree).
        threads (int): Number of threads to use for IQ-TREE.
    """
    try:
        # Validate input
        aln_file = validate_file_path(aln, must_exist=True)
        
        # Validate threads parameter
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("Threads must be a positive integer")
        
        # Check if alignment file is not empty
        if os.path.getsize(aln_file) > 0:
            if method == "fasttree":
                build_tree_fasttree(str(aln_file), treeout)
            elif method == "iqtree":
                # Use our Python IQ-TREE wrapper
                build_tree_python_iqtree(str(aln_file), treeout, threads)
            
            if os.path.exists(treeout):
                print(f"Tree built successfully with {method} and saved to {treeout}.")
            else:
                print(f"Error: {treeout} not found after running {method}.")
        else:
            print(f"The alignment file {aln} is empty. Skipping tree building.")
            
    except (ValueError, FileNotFoundError, GVClassError, subprocess.CalledProcessError) as e:
        print(f"Error occurred while running {method}: {e}")
        raise click.ClickException(str(e))


if __name__ == "__main__":
    main()