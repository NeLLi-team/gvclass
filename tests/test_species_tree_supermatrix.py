"""Unit tests for species-tree Section 4: supermatrix construction.

Covers reversible leaf namespacing, gap-filling of missing (taxon, group) cells,
partition offset correctness, group/leaf dropping, and the witchi -> pytrimal ->
unpruned trimming fallback chain.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

import src.core.species_tree.supermatrix as sm
from src.core.species_tree.supermatrix import (
    build_supermatrix,
    denamespace,
    is_query_leaf,
    namespace_query,
)


def _read_fasta(path: Path):
    out = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def test_namespacing_is_reversible() -> None:
    assert namespace_query("PkV-RF01") == "QUERY__PkV-RF01"
    assert is_query_leaf("QUERY__PkV-RF01")
    assert not is_query_leaf("NCLDV__G1")
    assert denamespace("QUERY__PkV-RF01") == "PkV-RF01"
    # A reference leaf keeps its id; a query stem that itself starts with a domain
    # prefix still round-trips because the QUERY__ namespace is stripped exactly once.
    assert denamespace("NCLDV__G1") == "NCLDV__G1"
    assert denamespace(namespace_query("NCLDV__weird")) == "NCLDV__weird"


def _identity_patches(monkeypatch) -> None:
    """Make align + trim identity so block widths equal the raw input lengths."""
    monkeypatch.setattr(
        sm,
        "align_sequences_pyfamsa",
        lambda records, threads: [(r.id, str(r.seq)) for r in records],
    )
    monkeypatch.setattr(
        sm,
        "_trim_group_alignment",
        lambda aln_path, threads, expected, method="witchi": sm._read_fasta_dict(aln_path),
    )


def test_build_supermatrix_gaps_and_partitions(tmp_path, monkeypatch) -> None:
    _identity_patches(monkeypatch)
    # grpA width 4, grpB width 6. r2 is missing grpB -> gap-filled.
    taxa = {
        "QUERY__q1": {"grpA": "AAAA", "grpB": "BBBBBB"},
        "NCLDV__r1": {"grpA": "CCCC", "grpB": "DDDDDD"},
        "NCLDV__r2": {"grpA": "EEEE"},
        "NCLDV__r3": {"grpA": "FFFF", "grpB": "GGGGGG"},
    }
    result = build_supermatrix(
        taxa, ["grpA", "grpB"], tmp_path / "work", tmp_path / "out", min_taxa_per_group=3
    )
    assert result is not None
    assert result.included_groups == ["grpA", "grpB"]
    assert result.group_widths == {"grpA": 4, "grpB": 6}

    rows = _read_fasta(result.supermatrix_faa)
    assert set(rows) == {"QUERY__q1", "NCLDV__r1", "NCLDV__r2", "NCLDV__r3"}
    # All rows are the concatenated width (4 + 6 = 10).
    assert {len(s) for s in rows.values()} == {10}
    # r2 missing grpB -> its second block is all gaps.
    assert rows["NCLDV__r2"] == "EEEE" + "-" * 6
    # Query block present in both partitions.
    assert rows["QUERY__q1"] == "AAAA" + "BBBBBB"

    partitions = result.partitions.read_text().strip().splitlines()
    assert partitions == ["AA, grpA = 1-4", "AA, grpB = 5-10"]


def test_group_with_too_few_taxa_is_dropped(tmp_path, monkeypatch) -> None:
    _identity_patches(monkeypatch)
    # grpB has only 2 taxa -> dropped; grpA has 3 -> kept.
    taxa = {
        "NCLDV__r1": {"grpA": "AAAA", "grpB": "BBBB"},
        "NCLDV__r2": {"grpA": "CCCC", "grpB": "DDDD"},
        "NCLDV__r3": {"grpA": "EEEE"},
    }
    result = build_supermatrix(
        taxa, ["grpA", "grpB"], tmp_path / "work", tmp_path / "out", min_taxa_per_group=3
    )
    assert result is not None
    assert result.included_groups == ["grpA"]
    assert "grpB" in result.dropped_groups
    rows = _read_fasta(result.supermatrix_faa)
    assert {len(s) for s in rows.values()} == {4}


def test_taxon_absent_from_all_surviving_groups_excluded(tmp_path, monkeypatch) -> None:
    _identity_patches(monkeypatch)
    # r9 appears only in grpZ which is dropped (1 taxon) -> r9 excluded entirely.
    taxa = {
        "NCLDV__r1": {"grpA": "AAAA"},
        "NCLDV__r2": {"grpA": "CCCC"},
        "NCLDV__r3": {"grpA": "EEEE"},
        "NCLDV__r9": {"grpZ": "ZZZZ"},
    }
    result = build_supermatrix(
        taxa, ["grpA", "grpZ"], tmp_path / "work", tmp_path / "out", min_taxa_per_group=3
    )
    assert result is not None
    assert result.included_groups == ["grpA"]
    rows = _read_fasta(result.supermatrix_faa)
    assert "NCLDV__r9" not in rows
    assert "NCLDV__r9" in result.dropped_leaves


def test_returns_none_when_no_group_survives(tmp_path, monkeypatch) -> None:
    _identity_patches(monkeypatch)
    taxa = {"NCLDV__r1": {"grpA": "AAAA"}, "NCLDV__r2": {"grpA": "CCCC"}}
    result = build_supermatrix(
        taxa, ["grpA"], tmp_path / "work", tmp_path / "out", min_taxa_per_group=3
    )
    assert result is None


# ---------------------------------------------------------------------------
# Trimming fallback chain
# ---------------------------------------------------------------------------

_ALN = (
    ">QUERY__q1\nMKTAYIAKQR-ISFVKSHFSR\n"
    ">NCLDV__r1\nMKTAYIAKQRQISFVKSHFSR\n"
    ">NCLDV__r2\nMKQAYIAKQR-ISFVKSHFSR\n"
)
_ALN_TAXA = {"QUERY__q1", "NCLDV__r1", "NCLDV__r2"}


def test_build_supermatrix_uses_trimmed_width(tmp_path, monkeypatch) -> None:
    """Gap-fill and partition offsets use the TRIMMED block width, not raw."""
    monkeypatch.setattr(
        sm,
        "align_sequences_pyfamsa",
        lambda records, threads: [(r.id, str(r.seq)) for r in records],
    )

    def narrow_trim(aln_path, threads, expected, method="witchi"):
        # Drop the last two columns from every row (raw 6 -> trimmed 4).
        return {k: v[:-2] for k, v in sm._read_fasta_dict(aln_path).items()}

    monkeypatch.setattr(sm, "_trim_group_alignment", narrow_trim)
    taxa = {
        "NCLDV__r1": {"grpA": "AAAAAA", "grpB": "TTTTTT"},
        "NCLDV__r2": {"grpA": "CCCCCC", "grpB": "GGGGGG"},
        "NCLDV__r3": {"grpA": "EEEEEE", "grpB": "CCCCCC"},
        "NCLDV__r4": {"grpA": "DDDDDD"},  # missing grpB
    }
    result = build_supermatrix(
        taxa, ["grpA", "grpB"], tmp_path / "work", tmp_path / "out", min_taxa_per_group=3
    )
    assert result is not None
    assert result.group_widths == {"grpA": 4, "grpB": 4}  # trimmed, not raw 6
    rows = _read_fasta(result.supermatrix_faa)
    assert {len(s) for s in rows.values()} == {8}
    assert rows["NCLDV__r4"] == "DDDD" + "-" * 4  # gap-fill at trimmed width
    assert result.partitions.read_text().strip().splitlines() == [
        "AA, grpA = 1-4",
        "AA, grpB = 5-8",
    ]


@pytest.mark.skipif(shutil.which("witchi") is None, reason="witchi not on PATH")
def test_trim_real_witchi_path(tmp_path) -> None:
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    # Real witchi run -> uniform, non-empty block over exactly the three taxa.
    assert set(block) == _ALN_TAXA
    assert sm._uniform_nonempty(block)


def test_trim_falls_back_to_pytrimal(tmp_path, monkeypatch) -> None:
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    monkeypatch.setattr(sm, "_witchi_prune", lambda aln_path, threads: None)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    assert set(block) == _ALN_TAXA
    assert sm._uniform_nonempty(block)


def test_trim_falls_back_to_unpruned(tmp_path, monkeypatch) -> None:
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    monkeypatch.setattr(sm, "_witchi_prune", lambda aln_path, threads: None)

    def _boom(aln_path, out_path):
        raise RuntimeError("pytrimal unavailable")

    monkeypatch.setattr(sm, "_pytrimal_trim", _boom)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    # Unpruned: identical to the raw alignment.
    assert block == sm._read_fasta_dict(aln)


def test_trim_rejects_empty_witchi_output(tmp_path, monkeypatch) -> None:
    """A 0-column witchi result is rejected so pruning never empties a block."""
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)

    def fake_witchi(aln_path, threads):  # writes after the stale-clear step
        out = aln_path.parent / "grp.aln_empty_pruned.fasta"
        out.write_text(">QUERY__q1\n\n>NCLDV__r1\n\n>NCLDV__r2\n\n")
        return out

    monkeypatch.setattr(sm, "_witchi_prune", fake_witchi)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    assert sm._uniform_nonempty(block)  # pytrimal fallback used
    assert set(block) == _ALN_TAXA


def test_trim_rejects_block_with_wrong_taxa(tmp_path, monkeypatch) -> None:
    """A pruned file whose taxa don't match the input (e.g. stale) is rejected."""
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)

    def fake_witchi(aln_path, threads):
        out = aln_path.parent / "grp.aln_wrong_pruned.fasta"
        out.write_text(">STALE_x\nMKTA\n>STALE_y\nMKTA\n>STALE_z\nMKTA\n")
        return out

    monkeypatch.setattr(sm, "_witchi_prune", fake_witchi)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    # Stale taxa rejected -> pytrimal/unpruned with the correct taxa.
    assert set(block) == _ALN_TAXA


def test_trim_rejects_pytrimal_dropping_taxon(tmp_path, monkeypatch) -> None:
    """pytrimal automated1 can REMOVE an all-gap sequence; the taxon-set check
    rejects such a block so rows are never silently corrupted."""
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    monkeypatch.setattr(sm, "_witchi_prune", lambda aln_path, threads: None)

    def drop_taxon(aln_path, out_path):
        out_path.write_text(">QUERY__q1\nMKTA\n>NCLDV__r1\nMKTA\n")  # NCLDV__r2 dropped
        return out_path

    monkeypatch.setattr(sm, "_pytrimal_trim", drop_taxon)
    block = sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    # Dropped-taxon block rejected -> unpruned, all three taxa intact.
    assert set(block) == _ALN_TAXA


def test_trim_method_none_returns_unpruned(tmp_path) -> None:
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    block = sm._trim_group_alignment(
        aln, threads=1, expected_leaves=set(_ALN_TAXA), method="none"
    )
    assert block == sm._read_fasta_dict(aln)  # raw, untrimmed


def test_trim_method_pytrimal_skips_witchi(tmp_path, monkeypatch) -> None:
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    called = {"witchi": False}

    def _spy(aln_path, threads):
        called["witchi"] = True
        return None

    monkeypatch.setattr(sm, "_witchi_prune", _spy)
    block = sm._trim_group_alignment(
        aln, threads=1, expected_leaves=set(_ALN_TAXA), method="pytrimal"
    )
    assert called["witchi"] is False  # witchi not invoked
    assert set(block) == _ALN_TAXA


def test_build_supermatrix_trim_method_none(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        sm,
        "align_sequences_pyfamsa",
        lambda records, threads: [(r.id, str(r.seq)) for r in records],
    )
    taxa = {
        "NCLDV__r1": {"grpA": "AAAAAA"},
        "NCLDV__r2": {"grpA": "CCCCCC"},
        "NCLDV__r3": {"grpA": "EEEEEE"},
    }
    result = build_supermatrix(
        taxa, ["grpA"], tmp_path / "w", tmp_path / "o", min_taxa_per_group=3, trim_method="none"
    )
    assert result is not None
    assert result.group_widths == {"grpA": 6}  # untrimmed


def test_trim_clears_stale_pruned_artifacts(tmp_path, monkeypatch) -> None:
    """Stale <stem>_*_pruned.fasta is removed before trimming so the glob in the
    real witchi path cannot pick it up."""
    aln = tmp_path / "grp.aln.fasta"
    aln.write_text(_ALN)
    stale = tmp_path / "grp.aln_stale_pruned.fasta"
    stale.write_text(">STALE_x\nMKTA\n")
    monkeypatch.setattr(sm, "_witchi_prune", lambda aln_path, threads: None)
    sm._trim_group_alignment(aln, threads=1, expected_leaves=set(_ALN_TAXA))
    assert not stale.exists()
