"""Query leaf matching, nearest neighbours and rank indices of marker trees."""

from __future__ import annotations

from pathlib import Path

from src.core.species_tree.protein_gathering import (
    select_query_representative_from_tree,
)
from src.core.tree_analysis import TreeAnalyzer, _is_query_leaf


def test_is_query_leaf_exact_and_delimited_matching() -> None:
    """_is_query_leaf should match exact, pipe-delimited, or single-underscore suffixes,
    without matching unrelated reference leaves that happen to share a prefix."""
    assert _is_query_leaf("G1", "G1") is True
    assert _is_query_leaf("G1|c1_1", "G1") is True
    assert _is_query_leaf("G1_1", "G1") is True
    assert _is_query_leaf("G10|c1_1", "G1") is False
    assert _is_query_leaf("G10_1", "G1") is False
    assert _is_query_leaf("NCLDV__G1_1", "NCLDV") is False
    assert _is_query_leaf("NCLDV__QUERY_1", "NCLDV__QUERY") is True


def test_find_nearest_neighbors_and_query_reps_ignore_prefix_colliding_refs(
    tmp_path: Path,
) -> None:
    """Query 'G1' must not match reference 'G10|c_1' as a query leaf."""
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "G10\td_Vir|p_NCLDV|c_Megaviricetes|o_Imitervirales|f_Mimiviridae|g_Mimi|s_Mimi_sp\n"
    )
    tree_file = tmp_path / "OG1.treefile"
    tree_file.write_text("((G1|q_1:0.1,G10|c_1:0.2):0.1,G10|c_2:0.5);\n")

    analyzer = TreeAnalyzer(labels)
    neighbors = analyzer.get_neighbors_from_tree(tree_file, query_id="G1")
    assert list(neighbors.keys()) == ["G1|q_1"]
    assert list(neighbors["G1|q_1"].keys()) == ["G10|c_1"]
    assert abs(neighbors["G1|q_1"]["G10|c_1"] - 0.3) < 1e-6

    rep = select_query_representative_from_tree(
        tree_file,
        "G1",
        keep_prefix="G10",
        query_lengths={"G1|q_1": 5, "G10|c_1": 12},
    )
    assert rep is not None
    assert rep[0] == "G1|q_1"


def test_marker_trees_use_7_rank_indices(
    tmp_path: Path,
) -> None:
    """7-rank labels (domain|phylum|class|order|family|genus|species) must map
    indices 3..6 to order, family, genus, species in the tree_nn file."""
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "NCLDV__G1\td_Viruses|p_Nucleocytoviricota|c_Megaviricetes|"
        "o_Imitervirales|f_Mimiviridae|g_Mimivirus|s_Mimivirus_bradfordmassiliense\n"
    )
    tree_dir = tmp_path / "trees"
    tree_dir.mkdir()
    tree_file = tree_dir / "GVOGm0003.treefile"
    tree_file.write_text("(QUERY|c1_1:0.1,NCLDV__G1|p1:0.2);\n")

    analyzer = TreeAnalyzer(labels)
    output_tsv = tmp_path / "nearest_neighbors.tsv"
    all_neighbors = analyzer.process_marker_trees(tree_dir, "QUERY", output_tsv)
    assert "GVOGm0003" in all_neighbors

    content = output_tsv.read_text()
    assert "NCLDV;o_Imitervirales;f_Mimiviridae;g_Mimivirus" in content
