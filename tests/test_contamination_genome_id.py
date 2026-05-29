"""M11: per-contig classifier must resolve pipeless PHAGE reference genome IDs.

_parse_subject_genome_id previously kept the trailing _<gene> suffix on pipeless
reference leaves (e.g. PHAGE__...__GCA-003814125-1_7), so labels.get() missed
them and their viral votes were silently dropped.
"""

from __future__ import annotations

from src.core.contamination_scoring import ContaminationScorer


def _scorer(labels):
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.labels = labels
    return scorer


def test_genome_id_resolution_handles_pipeless_phage():
    labels = {
        "PHAGE__VARDNA__GCA-003814125-1": {"domain": "PHAGE", "order": "o", "genus": "g"},
        "NCLDV__Mimi_1": {"domain": "NCLDV", "order": "Imitervirales", "genus": "g335"},
    }
    scorer = _scorer(labels)
    # Pipeless PHAGE protein -> strip trailing _7 to reach the genome label.
    assert (
        scorer._parse_subject_genome_id("PHAGE__VARDNA__GCA-003814125-1_7")
        == "PHAGE__VARDNA__GCA-003814125-1"
    )
    # A genome id that itself ends in _1 AND is a label must be kept as-is.
    assert scorer._parse_subject_genome_id("NCLDV__Mimi_1") == "NCLDV__Mimi_1"
    # Pipe form resolves to the pre-pipe genome.
    assert scorer._parse_subject_genome_id("NCLDV__Mimi_1|Ga0_5") == "NCLDV__Mimi_1"


def test_per_protein_dominant_lineage_counts_pipeless_phage_vote():
    labels = {
        "PHAGE__VARDNA__GCA-003814125-1": {
            "domain": "PHAGE",
            "order": "Caudovirales",
            "genus": "gg",
        }
    }
    scorer = _scorer(labels)
    tree_nn = {"markerA": {"prot1": ["PHAGE__VARDNA__GCA-003814125-1_7"]}}
    result = scorer._per_protein_dominant_lineage("prot1", tree_nn)
    assert result == {"domain": "PHAGE", "order": "Caudovirales", "genus": "gg"}


def test_existing_pipe_neighbor_unchanged():
    labels = {
        "NCLDV__Mimi_1": {"domain": "NCLDV", "order": "Imitervirales", "genus": "g335"}
    }
    scorer = _scorer(labels)
    tree_nn = {"markerA": {"prot1": ["NCLDV__Mimi_1"]}}
    result = scorer._per_protein_dominant_lineage("prot1", tree_nn)
    assert result["order"] == "Imitervirales"
