"""M11: per-contig classifier must resolve pipeless PHAGE reference genome IDs.

The subject resolver once kept the trailing _<gene> suffix on pipeless reference
leaves (e.g. PHAGE__...__GCA-003814125-1_7), so the label lookup missed them and
their viral votes were silently dropped.
"""

from __future__ import annotations

from src.core.contamination_scoring import ContaminationScorer
from src.core.labels import LabelIndex


def _scorer(taxonomy: dict[str, list[str]]) -> ContaminationScorer:
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.labels = LabelIndex.from_taxonomy(taxonomy)
    return scorer


def test_subject_taxonomy_resolves_pipeless_phage_and_pipe_forms():
    scorer = _scorer(
        {
            "PHAGE__VARDNA__GCA-003814125-1": ["PHAGE", "p", "c", "o", "f", "g", "s"],
            "NCLDV__Mimi_1": ["NCLDV", "p", "c", "Imitervirales", "Mimiviridae"],
        }
    )

    # Pipeless PHAGE protein -> strip trailing _7 to reach the genome label.
    assert scorer._subject_taxonomy("PHAGE__VARDNA__GCA-003814125-1_7") == (
        "PHAGE",
        "o",
        "f",
    )
    # A genome id that itself ends in _1 AND is a label keeps its label.
    assert scorer._subject_taxonomy("NCLDV__Mimi_1") == (
        "NCLDV",
        "Imitervirales",
        "Mimiviridae",
    )
    # Pipe form resolves to the pre-pipe genome.
    assert scorer._subject_taxonomy("NCLDV__Mimi_1|Ga0_5") == (
        "NCLDV",
        "Imitervirales",
        "Mimiviridae",
    )


def test_subject_taxonomy_of_an_unlabelled_subject_keeps_the_id_domain():
    scorer = _scorer({"NCLDV__a": ["NCLDV"]})

    assert scorer._subject_taxonomy("BAC__unknown|c_3") == ("BAC", "", "")
    assert scorer._subject_taxonomy("plain_7") == ("NA", "", "")


def test_per_protein_dominant_lineage_counts_pipeless_phage_vote():
    scorer = _scorer(
        {
            "PHAGE__VARDNA__GCA-003814125-1": [
                "PHAGE",
                "p",
                "c",
                "Caudovirales",
                "f",
                "gg",
                "s",
            ]
        }
    )
    tree_nn = {"markerA": {"prot1": {"PHAGE__VARDNA__GCA-003814125-1_7": 0.1}}}

    result = scorer._per_protein_dominant_lineage("prot1", tree_nn)

    assert result == {"domain": "PHAGE", "order": "Caudovirales", "genus": "gg"}


def test_existing_pipe_neighbor_unchanged():
    scorer = _scorer(
        {"NCLDV__Mimi_1": ["NCLDV", "p", "c", "Imitervirales", "f", "g335", "s"]}
    )
    tree_nn = {"markerA": {"prot1": {"NCLDV__Mimi_1": 0.1}}}

    result = scorer._per_protein_dominant_lineage("prot1", tree_nn)

    assert result["order"] == "Imitervirales"
