"""Unit tests for the shared alignment/tree helpers (``src.core.alignment``).

Section 1 of the ``--species-tree`` feature extracts the FAMSA alignment and
VeryFastTree inference cores out of ``MarkerProcessor`` into ``src.core.alignment``
so the supermatrix path can reuse them. These tests pin the helper contract;
behaviour-preservation of the per-marker path itself is covered by the golden
example regression test.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from src.core.alignment import align_sequences_pyfamsa, run_veryfasttree


def _rec(seq_id: str, seq: str) -> SimpleNamespace:
    """A minimal SeqRecord-like object exposing ``.id`` and ``.seq``."""
    return SimpleNamespace(id=seq_id, seq=seq)


_SEQS = [
    _rec("alpha", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG"),
    _rec("beta", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSG"),
    _rec("gamma", "MKQAYIAKQRQISFVKSHFSRQLEDRLGLIEVQAPILSRVGDGTQDNLSG"),
    _rec("delta", "MKTAYIWKQRQISFVKSHFSRQLEERLGLIEVQACILSRVGDGTQDNLSG"),
]


def test_align_returns_equal_width_msa_preserving_taxa() -> None:
    aligned = align_sequences_pyfamsa(_SEQS, threads=1)

    # One row per input, ids preserved.
    assert {seq_id for seq_id, _ in aligned} == {"alpha", "beta", "gamma", "delta"}
    # All rows are the same (aligned) width.
    widths = {len(seq) for _, seq in aligned}
    assert len(widths) == 1
    # Alignment is at least as wide as the longest input (gaps only add width).
    assert widths.pop() >= max(len(r.seq) for r in _SEQS)


def test_align_accepts_str_convertible_seq() -> None:
    """``.seq`` need not be a str — anything str-convertible (e.g. Bio.Seq) works."""

    class _SeqLike:
        def __init__(self, s: str) -> None:
            self._s = s

        def __str__(self) -> str:
            return self._s

    recs = [
        SimpleNamespace(id="x", seq=_SeqLike("MKTAYIAK")),
        SimpleNamespace(id="y", seq=_SeqLike("MKTAYIAK")),
        SimpleNamespace(id="z", seq=_SeqLike("MKQAYIAK")),
    ]
    aligned = align_sequences_pyfamsa(recs, threads=1)
    # FAMSA preserves input order for these identical-length inputs.
    assert [seq_id for seq_id, _ in aligned] == ["x", "y", "z"]


def test_run_veryfasttree_emits_parseable_newick(tmp_path: Path) -> None:
    aligned = align_sequences_pyfamsa(_SEQS, threads=1)
    aln_path = tmp_path / "aln.fasta"
    with open(aln_path, "w") as handle:
        for seq_id, seq in aligned:
            handle.write(f">{seq_id}\n{seq}\n")

    newick = run_veryfasttree(aln_path, threads=1)

    assert isinstance(newick, str) and newick.strip().endswith(";")

    ete3 = pytest.importorskip("ete3")
    tree = ete3.Tree(newick)
    assert {leaf.name for leaf in tree.get_leaves()} == {
        "alpha",
        "beta",
        "gamma",
        "delta",
    }
