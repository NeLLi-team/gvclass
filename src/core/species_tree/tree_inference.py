"""Species-tree inference (Section 5): run VeryFastTree on the supermatrix.

Thin wrapper over the shared :func:`src.core.alignment.run_veryfasttree` helper so
the supermatrix tree is built with the same inference settings as the per-marker
gene trees.
"""

from __future__ import annotations

from pathlib import Path

from src.core.alignment import run_iqtree, run_veryfasttree


def infer_species_tree(
    supermatrix_faa: Path, out_tree: Path, threads: int = 4, tree_method: str = "veryfasttree"
) -> Path:
    """Infer the concatenated species tree from the supermatrix; write Newick.

    Uses VeryFastTree by default (fast); set ``tree_method="iqtree"`` to build the
    species tree with IQ-TREE (``Q.pfam+R10+F``; slower, optionally with ultrafast
    bootstrap via ``--iqtree-mode ufboot``). Returns the tree path; raises if none.
    """
    if tree_method in ("veryfasttree", "fasttree"):
        newick = run_veryfasttree(supermatrix_faa, threads)
        if not newick or not newick.strip():
            raise RuntimeError(f"VeryFastTree produced no tree for {supermatrix_faa}")
        out_tree.write_text(newick)
        return out_tree
    newick = run_iqtree(supermatrix_faa, out_tree, threads)
    if not newick or not newick.strip():
        raise RuntimeError(f"IQ-TREE produced no tree for {supermatrix_faa}")
    return out_tree
