"""Unit tests for species-tree Section 3: representative-protein gathering.

Covers ``ReferenceProteinIndex`` (genome resolution, longest-protein, fetch),
query-representative selection from the dedicated trees (ref-closest paralog +
tie-break), reference-representative gathering (placed paralog preserved, longest
used only to fill additional groups), and the >= min_markers (3/8) filter.
"""

from __future__ import annotations

from pathlib import Path

from src.core.tree_analysis import TreeAnalyzer
from src.core.species_tree.neighbor_selection import NeighborHit
from src.core.species_tree.protein_gathering import (
    ReferenceProteinIndex,
    gather_query_representatives,
    gather_reference_representatives,
    select_query_representative_from_tree,
)


def _labels(tmp_path: Path) -> Path:
    path = tmp_path / "labels.tsv"
    path.write_text(
        "\n".join(
            f"NCLDV__{g}\tNCLDV|Nucleocytoviricota|Imitervirales|f|g|s_{g}"
            for g in ("G1", "G2", "G3")
        )
        + "\n"
    )
    return path


def _seq(n: int) -> str:
    return "M" + "A" * (n - 1)


def _write_tree(tmp_path: Path, name: str, newick: str) -> Path:
    path = tmp_path / name
    path.write_text(newick.strip() + "\n")
    return path


# ---------------------------------------------------------------------------
# ReferenceProteinIndex
# ---------------------------------------------------------------------------


def test_reference_index_genome_resolution_and_longest(tmp_path: Path) -> None:
    analyzer = TreeAnalyzer(_labels(tmp_path))
    faa = tmp_path / "grpA.faa"
    faa.write_text(
        f">NCLDV__G1|cA_1\n{_seq(10)}\n"
        f">NCLDV__G1|cA_2\n{_seq(20)}\n"
        f">NCLDV__G3|cA_1\n{_seq(15)}\n"
    )
    index = ReferenceProteinIndex(faa, analyzer)

    assert set(index.proteins_for("NCLDV__G1")) == {"NCLDV__G1|cA_1", "NCLDV__G1|cA_2"}
    # Longest by length -> cA_2 (20 > 10).
    assert index.longest_protein("NCLDV__G1") == "NCLDV__G1|cA_2"
    assert index.longest_protein("NCLDV__G3") == "NCLDV__G3|cA_1"
    assert index.longest_protein("NCLDV__ABSENT") is None
    assert index.sequence("NCLDV__G1|cA_1") == _seq(10)
    assert "NCLDV__G1|cA_1" in index
    index.close()


def test_reference_index_longest_tie_breaks_lexicographically(tmp_path: Path) -> None:
    analyzer = TreeAnalyzer(_labels(tmp_path))
    faa = tmp_path / "grp.faa"
    # Two equal-length proteins -> lexicographically smaller id wins.
    faa.write_text(f">NCLDV__G1|cB_2\n{_seq(20)}\n>NCLDV__G1|cB_1\n{_seq(20)}\n")
    index = ReferenceProteinIndex(faa, analyzer)
    assert index.longest_protein("NCLDV__G1") == "NCLDV__G1|cB_1"
    index.close()


# ---------------------------------------------------------------------------
# Query representative selection
# ---------------------------------------------------------------------------


def test_query_rep_picks_ref_closest_paralog(tmp_path: Path) -> None:
    # p_1 sits next to an NCLDV ref (dist 0.2); p_2 is far (min 0.6).
    tree = _write_tree(
        tmp_path,
        "grp.treefile",
        "((PkV-RF01|p_1:0.1,NCLDV__G1_1:0.1):0.1,"
        "(PkV-RF01|p_2:0.5,NCLDV__G2_1:0.1):0.1);",
    )
    chosen = select_query_representative_from_tree(tree, "PkV-RF01", keep_prefix="NCLDV__")
    assert chosen is not None
    protein, dist = chosen
    assert protein == "PkV-RF01|p_1"
    assert dist == 0.2


def test_query_rep_tie_breaks_by_length_then_id(tmp_path: Path) -> None:
    # Both paralogs equidistant (0.4) to the single NCLDV ref.
    tree = _write_tree(
        tmp_path,
        "tie.treefile",
        "(NCLDV__G1_1:0.2,(PkV-RF01|p_1:0.1,PkV-RF01|p_2:0.1):0.1);",
    )
    # Length tie-break: p_2 longer -> p_2 wins despite higher id.
    chosen = select_query_representative_from_tree(
        tree,
        "PkV-RF01",
        keep_prefix="NCLDV__",
        query_lengths={"PkV-RF01|p_1": 100, "PkV-RF01|p_2": 200},
    )
    assert chosen is not None and chosen[0] == "PkV-RF01|p_2"

    # No lengths -> deterministic id tie-break (p_1 < p_2).
    chosen_id = select_query_representative_from_tree(tree, "PkV-RF01", keep_prefix="NCLDV__")
    assert chosen_id is not None and chosen_id[0] == "PkV-RF01|p_1"


def test_query_rep_none_without_query_or_ref(tmp_path: Path) -> None:
    no_ref = _write_tree(tmp_path, "noref.treefile", "(PkV-RF01|p_1:0.1,BAC__B_1:0.2);")
    assert select_query_representative_from_tree(no_ref, "PkV-RF01") is None
    no_q = _write_tree(tmp_path, "noq.treefile", "(NCLDV__G1_1:0.1,NCLDV__G2_1:0.2);")
    assert select_query_representative_from_tree(no_q, "PkV-RF01") is None


def test_gather_query_representatives_min_markers(tmp_path: Path) -> None:
    trees = {
        g: _write_tree(
            tmp_path,
            f"{g}.treefile",
            f"((PkV-RF01|{g}_1:0.1,NCLDV__G1_1:0.1):0.1,NCLDV__G2_1:0.3);",
        )
        for g in ("grpA", "grpB", "grpC")
    }
    reps = gather_query_representatives(trees, "PkV-RF01", min_markers=3)
    assert reps is not None
    assert set(reps) == {"grpA", "grpB", "grpC"}
    assert reps["grpA"] == "PkV-RF01|grpA_1"

    # Only two markers -> below the 3/8 floor -> excluded.
    two = {k: trees[k] for k in ("grpA", "grpB")}
    assert gather_query_representatives(two, "PkV-RF01", min_markers=3) is None


# ---------------------------------------------------------------------------
# Reference representative gathering + >= min_markers filter
# ---------------------------------------------------------------------------


def _ref_indexes(tmp_path: Path, analyzer: TreeAnalyzer):
    (tmp_path / "grpA.faa").write_text(
        f">NCLDV__G1|cA_1\n{_seq(10)}\n>NCLDV__G1|cA_2\n{_seq(20)}\n"
        f">NCLDV__G3|cA_1\n{_seq(15)}\n"
    )
    (tmp_path / "grpB.faa").write_text(
        f">NCLDV__G1|cB_1\n{_seq(12)}\n>NCLDV__G1|cB_2\n{_seq(25)}\n"
        f">NCLDV__G3|cB_1\n{_seq(12)}\n"
    )
    (tmp_path / "grpC.faa").write_text(
        f">NCLDV__G1|cC_1\n{_seq(12)}\n>NCLDV__G3|cC_1\n{_seq(12)}\n"
    )
    (tmp_path / "grpD.faa").write_text(f">NCLDV__G2|cD_1\n{_seq(12)}\n")
    return {
        g: ReferenceProteinIndex(tmp_path / f"{g}.faa", analyzer)
        for g in ("grpA", "grpB", "grpC", "grpD")
    }


def test_reference_rep_preserves_placed_paralog_and_fills_longest(tmp_path: Path) -> None:
    analyzer = TreeAnalyzer(_labels(tmp_path))
    indexes = _ref_indexes(tmp_path, analyzer)
    groups = ["grpA", "grpB", "grpC", "grpD"]

    # G1 placed in grpA as the SHORT protein cA_1 (not the longer cA_2).
    hits = [
        NeighborHit("NCLDV__G1", "NCLDV__G1|cA_1", "grpA", 0.2),
        NeighborHit("NCLDV__G2", "NCLDV__G2|cD_1", "grpD", 0.3),
    ]
    kept, dropped = gather_reference_representatives(hits, indexes, groups, min_markers=3)

    # G1 kept (grpA placed + grpB/grpC filled = 3 markers); G2 dropped (1 marker).
    assert set(kept) == {"NCLDV__G1"}
    assert "NCLDV__G2" in dropped

    g1 = kept["NCLDV__G1"]
    # Placed paralog preserved exactly (NOT replaced by the longer cA_2).
    assert g1["grpA"] == "NCLDV__G1|cA_1"
    # Additional groups filled with the genome's LONGEST protein.
    assert g1["grpB"] == "NCLDV__G1|cB_2"  # 25 > 12
    assert g1["grpC"] == "NCLDV__G1|cC_1"
    assert "grpD" not in g1  # G1 has no protein in grpD
    assert len(g1) == 3

    for index in indexes.values():
        index.close()


def test_reference_rep_genome_with_two_markers_dropped(tmp_path: Path) -> None:
    analyzer = TreeAnalyzer(_labels(tmp_path))
    indexes = _ref_indexes(tmp_path, analyzer)
    groups = ["grpA", "grpB", "grpC", "grpD"]
    # G3 present in grpA, grpB, grpC (3 markers) but only ever a neighbor in grpA.
    hits = [NeighborHit("NCLDV__G3", "NCLDV__G3|cA_1", "grpA", 0.2)]
    kept, dropped = gather_reference_representatives(hits, indexes, groups, min_markers=3)
    assert set(kept) == {"NCLDV__G3"}  # filled grpB+grpC -> 3/8 kept
    assert kept["NCLDV__G3"]["grpA"] == "NCLDV__G3|cA_1"
    assert len(kept["NCLDV__G3"]) == 3

    # Raise the floor to 4 -> now G3 (3 markers) is dropped.
    kept4, dropped4 = gather_reference_representatives(hits, indexes, groups, min_markers=4)
    assert kept4 == {}
    assert "NCLDV__G3" in dropped4

    for index in indexes.values():
        index.close()


def test_reference_rep_duplicate_hit_picks_nearest_order_independent(tmp_path: Path) -> None:
    """Combined mode: the same (genome, group) placed by two queries resolves to
    the NEAREST paralog regardless of hit order (no first-query-wins)."""
    analyzer = TreeAnalyzer(_labels(tmp_path))
    indexes = _ref_indexes(tmp_path, analyzer)
    groups = ["grpA", "grpB", "grpC", "grpD"]

    far = NeighborHit("NCLDV__G1", "NCLDV__G1|cA_2", "grpA", 0.9)
    near = NeighborHit("NCLDV__G1", "NCLDV__G1|cA_1", "grpA", 0.2)

    kept_fwd, _ = gather_reference_representatives([far, near], indexes, groups, min_markers=3)
    kept_rev, _ = gather_reference_representatives([near, far], indexes, groups, min_markers=3)

    # Nearest placement (cA_1 @ 0.2) wins both ways — order independent.
    assert kept_fwd["NCLDV__G1"]["grpA"] == "NCLDV__G1|cA_1"
    assert kept_rev["NCLDV__G1"]["grpA"] == "NCLDV__G1|cA_1"

    for index in indexes.values():
        index.close()


def test_reference_rep_skips_unindexed_fill_group(tmp_path: Path) -> None:
    """A group in `groups` but absent from `ref_indexes` is silently skipped."""
    analyzer = TreeAnalyzer(_labels(tmp_path))
    indexes = _ref_indexes(tmp_path, analyzer)
    groups = ["grpA", "grpB", "grpC", "grpZ"]  # grpZ has no index
    hits = [NeighborHit("NCLDV__G1", "NCLDV__G1|cA_1", "grpA", 0.2)]

    kept, _ = gather_reference_representatives(hits, indexes, groups, min_markers=3)
    assert set(kept) == {"NCLDV__G1"}
    assert "grpZ" not in kept["NCLDV__G1"]  # unindexed group skipped, no error
    assert len(kept["NCLDV__G1"]) == 3

    for index in indexes.values():
        index.close()


def test_reference_index_context_manager_and_idempotent_close(tmp_path: Path) -> None:
    analyzer = TreeAnalyzer(_labels(tmp_path))
    faa = tmp_path / "grp.faa"
    faa.write_text(f">NCLDV__G1|c_1\n{_seq(10)}\n")
    with ReferenceProteinIndex(faa, analyzer) as index:
        assert index.proteins_for("NCLDV__G1") == ["NCLDV__G1|c_1"]
    # Exited -> closed; a second close() is a no-op (idempotent), no exception.
    index.close()
