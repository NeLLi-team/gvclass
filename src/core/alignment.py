"""Shared multiple-sequence-alignment and tree-inference primitives.

These two helpers are the reusable cores extracted **verbatim** (behaviour-
preserving) from :class:`src.core.marker_processing.MarkerProcessor`. They are
shared by two callers:

* the per-marker gene-tree path (``MarkerProcessor.align_and_trim`` /
  ``MarkerProcessor.build_tree``), and
* the ``--species-tree`` supermatrix path (``src.core.species_tree`` —
  Sections 4 & 5),

so the species-tree feature reuses the exact FAMSA parameters and VeryFastTree
invocation as the validated gene-tree pipeline instead of duplicating them.

Both functions are intentionally thin: callers own sequence selection, the
"too few sequences" guard, file I/O naming, stop-codon cleaning, and trimming.
Keeping these helpers free of those concerns is what lets MarkerProcessor stay
byte-for-byte identical while the supermatrix builder applies its own policy.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Tuple, Union

import pyfamsa
import veryfasttree


def align_sequences_pyfamsa(
    records: Iterable, threads: int = 4
) -> List[Tuple[str, str]]:
    """Align protein records with FAMSA (pyfamsa) and return the aligned MSA.

    This is the alignment core lifted from ``MarkerProcessor.align_and_trim``;
    it uses the identical ``pyfamsa.Aligner(threads=...)`` call so output is
    unchanged from the validated gene-tree path.

    Args:
        records: Iterable of sequence records, each exposing ``.id`` (str) and
            ``.seq`` (str-convertible) — e.g. ``Bio.SeqRecord.SeqRecord``.
        threads: Worker threads passed straight to ``pyfamsa.Aligner`` (kept raw
            to match the historical per-marker call; FAMSA output is
            thread-count independent).

    Returns:
        List of ``(sequence_id, aligned_sequence)`` tuples in FAMSA's output
        order, with ids/sequences decoded to ``str``. Callers write these out
        however they need (the gene-tree path writes ``>{id}\\n{seq}\\n``).
    """
    famsa_seqs = [
        pyfamsa.Sequence(str(record.id).encode(), str(record.seq).encode())
        for record in records
    ]

    aligner = pyfamsa.Aligner(threads=threads)
    msa = aligner.align(famsa_seqs)

    return [(seq.id.decode(), seq.sequence.decode()) for seq in msa]


def run_veryfasttree(alignment_path: Union[str, Path], threads: int = 4) -> str:
    """Infer an approximately-maximum-likelihood tree with VeryFastTree.

    This is the inference core lifted from ``MarkerProcessor.build_tree``
    (``fasttree`` branch); it uses the identical ``veryfasttree.run`` call so
    output matches the validated gene-tree path. Caller owns the
    "too few sequences" guard, writing the newick, and result validation.

    Args:
        alignment_path: Path to the (trimmed) alignment in a format VeryFastTree
            reads (FASTA here).
        threads: Worker threads; clamped to ``>= 1`` exactly as the per-marker
            call did.

    Returns:
        The Newick tree string returned by ``veryfasttree.run``.
    """
    return veryfasttree.run(
        alignment=str(alignment_path),
        quiet=True,
        threads=max(1, int(threads)),
    )


#: Default IQ-TREE model for the per-marker gene trees and the species tree.
IQTREE_MODEL = "Q.pfam+R10+F"

#: IQ-TREE search mode for the SPECIES tree, set once per run via
#: ``configure_iqtree()`` (the pipeline is a single-process ThreadPoolExecutor so
#: all worker threads share it). Per-marker gene trees always use ``--fast``.
_IQTREE_MODE = "fast"


def configure_iqtree(mode: str) -> None:
    """Set the species-tree IQ-TREE search mode for this run.

    ``"fast"`` (default) or ``"ufboot"``; any other value falls back to ``"fast"``.
    """
    global _IQTREE_MODE
    _IQTREE_MODE = mode if mode in ("fast", "ufboot") else "fast"


def iqtree_search_args() -> List[str]:
    """IQ-TREE search args for the configured species-tree mode.

    ``fast``   -> ``--fast`` (quick FastTree-like search; default)
    ``ufboot`` -> ``-B 1000 -bnni`` (ultrafast bootstrap, 1000 replicates; IQ-TREE
    also writes a ``<prefix>.contree`` consensus tree with branch support).
    """
    if _IQTREE_MODE == "ufboot":
        return ["-B", "1000", "-bnni"]
    return ["--fast"]


def run_iqtree(
    alignment_path: Union[str, Path],
    out_tree: Union[str, Path],
    threads: int = 4,
    model: str = IQTREE_MODEL,
) -> str:
    """Infer a tree with IQ-TREE in ``--fast`` mode under ``model``.

    Runs ``iqtree --fast -m <model> -s <aln> -pre <prefix>`` (the prefix is
    derived from ``out_tree`` so IQ-TREE writes ``<prefix>.treefile == out_tree``)
    and returns the Newick string. Caller owns the too-few-sequences guard and
    result validation.
    """
    import subprocess

    out_tree = Path(out_tree)
    prefix = out_tree.with_suffix("")
    cmd = [
        "iqtree",
        *iqtree_search_args(),
        "-s",
        str(alignment_path),
        "-m",
        model,
        "-nt",
        str(max(1, int(threads))),
        "-pre",
        str(prefix),
        "-quiet",
    ]
    # Ultrafast bootstrap (-B 1000 -bnni) is far heavier than --fast; give it a
    # much larger ceiling so a real species-tree supermatrix can finish.
    timeout = 14400 if _IQTREE_MODE == "ufboot" else 1800
    subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=timeout)
    treefile = prefix.with_suffix(".treefile")
    return treefile.read_text() if treefile.exists() else ""
