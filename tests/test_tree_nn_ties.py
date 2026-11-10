"""Tied nearest neighbours: the tree reader keeps all, the vote uses their lineage."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest

from src.core.labels import LabelIndex
from src.core.summarize import labels_taxonomy
from src.core.tree_analysis import TreeAnalyzer

IM_19 = ("NCLDV", "Nucleocytoviricota", "Megaviricetes", "Imitervirales", "IM_19")
TAXONOMY = {
    "NCLDV__ref_a": [*IM_19, "g1787", "singleton"],
    "NCLDV__ref_b": [*IM_19, "g_other", "s_other"],
    "NCLDV__ref_twin": [*IM_19, "g1787", "singleton"],
    "NCLDV__ref_far": [
        *IM_19[:3],
        "Algavirales",
        "Phycodnaviridae",
        "g_far",
        "s_far",
    ],
    "BAC__ref_bac": ["BAC", "Pseudomonadota", "Gamma", "Entero", "Entero_f", "g", "s"],
    "PPV__ref_ppv": ["PPV", "Preplasmiviricota", "c", "o", "f", "g", "s"],
}
LABELS = LabelIndex.from_taxonomy(TAXONOMY)
FAMILY_LINEAGE = (*IM_19, "", "")
FULL_LINEAGE = (*IM_19, "g1787", "singleton")

TIED_NEWICKS = [
    "(query_1:0.1,ref_b_1:0.2,ref_a_1:0.2,ref_c_1:0.25);",
    "(query_1:0.1,ref_a_1:0.2,ref_b_1:0.2,ref_c_1:0.25);",
]


def _write_tree(tmp_path: Path, newick: str) -> tuple[Path, Path]:
    labels = tmp_path / "labels.tsv"
    labels.write_text(
        "ref_a\td_Vir|p_N|c_M|o_Imitervirales|f_Mimiviridae|g_G|s_S\n",
        encoding="utf-8",
    )
    tree_file = tmp_path / "OG1.treefile"
    tree_file.write_text(newick + "\n", encoding="utf-8")
    return labels, tree_file


@pytest.mark.parametrize("newick", TIED_NEWICKS)
def test_tree_reader_returns_every_nearest_neighbour_sorted_by_name(
    tmp_path: Path, newick: str
) -> None:
    labels, tree_file = _write_tree(tmp_path, newick)

    neighbors = TreeAnalyzer(labels).get_neighbors_from_tree(tree_file, "query")

    assert neighbors == {
        "query_1": {"ref_a_1": pytest.approx(0.3), "ref_b_1": pytest.approx(0.3)}
    }
    assert list(neighbors["query_1"]) == ["ref_a_1", "ref_b_1"]


@pytest.mark.parametrize("newick", TIED_NEWICKS)
def test_tree_nn_file_has_one_row_per_tied_neighbour(
    tmp_path: Path, newick: str
) -> None:
    labels, tree_file = _write_tree(tmp_path, newick)
    output_file = tmp_path / "query.tree_nn"

    TreeAnalyzer(labels).process_marker_trees(
        tmp_path, "query", output_file, tree_paths=[tree_file]
    )

    rows = [
        line.split("\t")[:4]
        for line in output_file.read_text(encoding="utf-8").splitlines()[2:]
    ]
    assert rows == [
        ["OG1", "query_1", "ref_a_1", "0.300000"],
        ["OG1", "query_1", "ref_b_1", "0.300000"],
    ]


def _collect(
    tree_nn_results: dict[str, dict[str, dict[str, float]]],
) -> tuple[dict[str, Counter], dict[str, Counter], list[float]]:
    return labels_taxonomy.collect_taxonomy_counts(LABELS, tree_nn_results)


def test_tied_references_vote_with_their_common_lineage() -> None:
    tree_nn = {"m1": {"q_1": {"NCLDV__ref_a_1": 0.2, "NCLDV__ref_b_1": 0.2}}}

    tax_counters, lineage_counters, distances = _collect(tree_nn)

    assert lineage_counters == {"m1": Counter({FAMILY_LINEAGE: 1})}
    assert tax_counters["order"] == Counter({"NCLDV__Imitervirales": 1})
    assert tax_counters["family"] == Counter({"NCLDV__IM_19": 1})
    assert tax_counters["genus"] == Counter()
    assert tax_counters["species"] == Counter()
    assert distances == [0.2]


def test_tied_references_with_the_same_lineage_change_nothing() -> None:
    tied = {"m1": {"q_1": {"NCLDV__ref_a_1": 0.2, "NCLDV__ref_twin_1": 0.2}}}
    single = {"m1": {"q_1": {"NCLDV__ref_a_1": 0.2}}}

    tax_counters, lineage_counters, distances = _collect(tied)

    assert (tax_counters, lineage_counters, distances) == _collect(single)
    assert lineage_counters == {"m1": Counter({FULL_LINEAGE: 1})}
    assert tax_counters["species"] == Counter({"NCLDV__singleton": 1})


def test_tied_references_from_different_domains_cast_no_lineage_vote() -> None:
    tree_nn = {"m1": {"q_1": {"NCLDV__ref_a_1": 0.2, "BAC__ref_bac_1": 0.2}}}

    tax_counters, lineage_counters, distances = _collect(tree_nn)

    assert lineage_counters == {}
    assert all(not counter for counter in tax_counters.values())
    assert distances == [0.2]


def test_unlabelled_references_do_not_take_part_in_a_tie() -> None:
    tree_nn = {
        "m1": {
            "q_1": {
                "UNKNOWN__x_1": 0.1,
                "NCLDV__ref_a_1": 0.2,
                "NCLDV__ref_b_1": 0.2,
                "NCLDV__ref_far_1": 0.4,
            }
        }
    }

    _, lineage_counters, distances = _collect(tree_nn)

    assert lineage_counters == {"m1": Counter({FAMILY_LINEAGE: 1})}
    assert distances == [0.2]


def test_ties_across_markers_credit_the_smallest_marker_name() -> None:
    tree_nn = {
        "m_c": {"q_1": {"NCLDV__ref_far_1": 0.3}},
        "m_b": {"q_1": {"NCLDV__ref_a_1": 0.2}},
        "m_a": {"q_1": {"NCLDV__ref_b_1": 0.2}},
    }

    _, lineage_counters, _ = _collect(tree_nn)

    assert lineage_counters == {"m_a": Counter({FAMILY_LINEAGE: 1})}


def test_a_closer_placement_in_another_marker_decides_alone() -> None:
    tree_nn = {
        "m_a": {"q_1": {"NCLDV__ref_a_1": 0.3, "NCLDV__ref_b_1": 0.3}},
        "m_b": {"q_1": {"NCLDV__ref_b_1": 0.1}},
    }

    _, lineage_counters, distances = _collect(tree_nn)

    assert lineage_counters == {"m_b": Counter({(*IM_19, "g_other", "s_other"): 1})}
    assert distances == [0.1]


@pytest.mark.parametrize("reverse", [False, True])
def test_counts_do_not_depend_on_insertion_order(reverse: bool) -> None:
    placements = [
        ("m_a", "q_1", "NCLDV__ref_a_1", 0.2),
        ("m_a", "q_1", "NCLDV__ref_b_1", 0.2),
        ("m_b", "q_1", "NCLDV__ref_twin_1", 0.2),
        ("m_b", "q_2", "NCLDV__ref_far_1", 0.5),
        ("m_b", "q_2", "NCLDV__ref_a_1", 0.5),
    ]
    if reverse:
        placements.reverse()
    tree_nn: dict[str, dict[str, dict[str, float]]] = {}
    for marker, protein, neighbor, distance in placements:
        tree_nn.setdefault(marker, {}).setdefault(protein, {})[neighbor] = distance

    tax_counters, lineage_counters, distances = _collect(tree_nn)

    assert lineage_counters == {
        "m_a": Counter({FAMILY_LINEAGE: 1}),
        "m_b": Counter(
            {("NCLDV", "Nucleocytoviricota", "Megaviricetes", "", "", "", ""): 1}
        ),
    }
    assert tax_counters["family"] == Counter({"NCLDV__IM_19": 1})
    assert sorted(distances) == [0.2, 0.5]


@pytest.mark.parametrize("reverse", [False, True])
def test_nearest_domain_of_a_tie_is_the_shared_domain_or_none(reverse: bool) -> None:
    same_domain = [("NCLDV__ref_a_1", 0.2), ("NCLDV__ref_far_1", 0.2)]
    mixed_domain = [("NCLDV__ref_a_1", 0.2), ("PPV__ref_ppv_1", 0.2)]
    if reverse:
        same_domain, mixed_domain = same_domain[::-1], mixed_domain[::-1]
    tree_nn = {"m1": {"q_same": dict(same_domain), "q_mixed": dict(mixed_domain)}}

    domains = labels_taxonomy.protein_nearest_domain(LABELS, tree_nn)

    assert domains == {"q_same": "NCLDV"}
