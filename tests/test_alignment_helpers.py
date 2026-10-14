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
from src.core.marker_processing import MarkerProcessor


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


def test_build_tree_skips_two_unaligned_sequences(tmp_path: Path) -> None:
    processor = MarkerProcessor("OG312", tmp_path / "db", tmp_path / "query")
    alignment = tmp_path / "OG312.mafft01"
    alignment.write_text(">query\nMKT\n>reference\nMKTAY\n", encoding="utf-8")
    tree = processor.tree_dir / "OG312.treefile"
    error = processor.tree_dir / "OG312.ERROR"
    tree.write_text("stale tree", encoding="utf-8")
    error.write_text("stale error", encoding="utf-8")

    result = processor.build_tree(alignment, threads=1)

    assert result is None
    assert not tree.exists()
    assert not error.exists()
    skipped = processor.tree_dir / "OG312.SKIPPED"
    assert "Sequences: 2" in skipped.read_text(encoding="utf-8")


def test_build_tree_rejects_three_unequal_length_sequences(tmp_path: Path) -> None:
    processor = MarkerProcessor("marker", tmp_path / "db", tmp_path / "query")
    alignment = tmp_path / "marker.mafft01"
    alignment.write_text(
        ">query\nMKT\n>reference_1\nMKTAY\n>reference_2\nMKTA\n",
        encoding="utf-8",
    )
    skipped = processor.tree_dir / "marker.SKIPPED"
    skipped.write_text("stale skip", encoding="utf-8")

    with pytest.raises(ValueError, match="same length"):
        processor.build_tree(alignment, threads=1)

    assert not skipped.exists()
    assert not (processor.tree_dir / "marker.treefile").exists()
