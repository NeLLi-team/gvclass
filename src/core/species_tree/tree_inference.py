"""Species-tree inference (Section 5): run VeryFastTree on the supermatrix.

Thin wrapper over the shared :func:`src.core.alignment.run_veryfasttree` helper so
the supermatrix tree is built with the same inference settings as the per-marker
gene trees.
"""

from __future__ import annotations

from pathlib import Path

from src.core.alignment import run_veryfasttree


def infer_species_tree(
    supermatrix_faa: Path, out_tree: Path, threads: int = 4
) -> Path:
    """Infer the concatenated species tree from the supermatrix; write Newick.

    Returns the tree path. Raises if VeryFastTree produces an empty result.
    """
    newick = run_veryfasttree(supermatrix_faa, threads)
    if not newick or not newick.strip():
        raise RuntimeError(f"VeryFastTree produced no tree for {supermatrix_faa}")
    out_tree.write_text(newick)
    return out_tree
