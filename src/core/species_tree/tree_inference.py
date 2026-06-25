"""Species-tree inference (Section 5): run VeryFastTree on the supermatrix.

Thin wrapper over the shared :func:`src.core.alignment.run_veryfasttree` helper so
the supermatrix tree is built with the same inference settings as the per-marker
gene trees.
"""

from __future__ import annotations

from pathlib import Path

from src.core.alignment import run_iqtree, run_veryfasttree


def infer_species_tree(
    supermatrix_faa: Path, out_tree: Path, threads: int = 4, tree_method: str = "iqtree"
) -> Path:
    """Infer the concatenated species tree from the supermatrix; write Newick.

    Uses IQ-TREE (``--fast``, Q.pfam+R10+F) by default, the same engine/model as
    the per-marker gene trees; set ``tree_method="fasttree"`` for VeryFastTree.
    Returns the tree path; raises if no tree is produced.
    """
    if tree_method == "fasttree":
        newick = run_veryfasttree(supermatrix_faa, threads)
        if not newick or not newick.strip():
            raise RuntimeError(f"VeryFastTree produced no tree for {supermatrix_faa}")
        out_tree.write_text(newick)
        return out_tree
    newick = run_iqtree(supermatrix_faa, out_tree, threads)
    if not newick or not newick.strip():
        raise RuntimeError(f"IQ-TREE produced no tree for {supermatrix_faa}")
    return out_tree
