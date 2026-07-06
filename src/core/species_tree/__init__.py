"""Opt-in NCLDV (and, via config, PPV/MIRUS) species-tree feature for gvclass.

This package builds a concatenated GVOG8 supermatrix species tree and assigns a
root-invariant tree-placement taxonomy per query. It is gated on gvclass calling
a genome NCLDV and is fully isolated from the production per-marker gene-tree
path, so enabling ``--species-tree`` never perturbs the standard classification
(``taxonomy_majority``) or the golden regression output.
"""
