"""Unit tests for the v1.4.3 Phase 2 per-contig taxonomic-purity classifier.

Covers the six classifier-level cases listed in
docs/plans/active/v1_4_3_contig_purity_contamination.md under the
"Tests" section.
"""

from __future__ import annotations

from pathlib import Path

import pytest


def _make_scorer(tmp_path: Path):
    """Build a ContaminationScorer with a synthetic labels table covering
    the cellular lineages the tests reference."""
    from src.core.contamination_scoring import ContaminationScorer
    from tests.conftest import stage_db_resources

    db = tmp_path / "db"
    stage_db_resources(db, markers=False)
    labels = db / "labels.tsv"
    # _load_labels expects genome_id<TAB>pipe_delimited_tax.
    # Pipe columns (7 fields): [0]=phylum, [1]=class, [2]=order, [3]=family, [4]=genus, [5]=species.
    # ContaminationScorer reads "order" from parts[3] (index after phylum|class|order);
    # see _load_labels: tax = parts[1].split("|"); labels[genome_id] = {..., "order": tax[3]}.
    # So for a (phylum, class, order, family, genus, species) sequence we pad the leading slot.
    labels.write_text(
        "EUK__Nannochloropsis_A\tEuk|Stramenopiles|Eustigmatophyceae|Eustigmatales|Monodopsidaceae|Nannochloropsis|Nannochloropsis_gaditana\n"
        "EUK__Nannochloropsis_B\tEuk|Stramenopiles|Eustigmatophyceae|Eustigmatales|Monodopsidaceae|Nannochloropsis|Nannochloropsis_oceanica\n"
        "EUK__Acanthocystis_A\tEuk|Centroheliozoa|Haptista|Acanthocystidae|Acanthocystidae|Acanthocystis|NA\n"
        "EUK__NA_A\tEuk|NA|NA|NA|NA|NA|NA\n"
        "EUK__Other_A\tEuk|NA|NA|OtherOrder|NA|NA|NA\n"
        "NCLDV__Mimi_1\tNCLDV|Nucleocytoviricota|Megaviricetes|Imitervirales|Mimiviridae|Mimivirus|Acanthamoeba_polyphaga_mimivirus\n"
        "BAC__Eco_A\tBac|Pseudomonadota|Gammaproteobacteria|Enterobacterales|Enterobacteriaceae|Escherichia|Escherichia_coli\n"
    )
    # Minimal completeness table (not used by the classifier).
    (db / "order_completeness.tab").write_text(
        "Order\tOrthogroups\tAverage_Percent\tStd_Percent\n"
    )
    return ContaminationScorer(db)


def _tree_nn(marker: str, protein: str, neighbor: str, distance: float = 0.5):
    """Build a single-row tree_nn_results slice."""
    return {marker: {protein: {neighbor: distance}}}


def _merge(*slices):
    out: dict = {}
    for slice_ in slices:
        for marker, proteins in slice_.items():
            out.setdefault(marker, {}).update(proteins)
    return out


def _blast_hit(
    query_id: str,
    subject_domain: str,
    subject_order: str = "",
    identity: float = 80.0,
    score: float = 120.0,
):
    from src.core.contamination_scoring import BlastHit

    return BlastHit(
        query_id=query_id,
        subject_id=f"{subject_domain}__synthetic",
        identity=identity,
        score=score,
        subject_domain=subject_domain,
        subject_order=subject_order,
        subject_family="",
    )


# ---------------------------------------------------------------------------
# Six classifier cases from the plan
# ---------------------------------------------------------------------------


def test_viral_bearing_contig_keeps_euk_markers_unflagged(tmp_path: Path) -> None:
    """Contig with one NCLDV marker + two EUK-resolving markers must be
    classified as viral_bearing, not cellular_coherent, regardless of
    how coherent the EUK markers look."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M_NCLDV", "p_viral", "NCLDV__Mimi_1"),
        _tree_nn("M_EUK_1", "p_euk_1", "EUK__Nannochloropsis_A"),
        _tree_nn("M_EUK_2", "p_euk_2", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        "p_viral": _blast_hit("p_viral", "NCLDV", "Imitervirales"),
        "p_euk_1": _blast_hit("p_euk_1", "EUK", "Eustigmatales", identity=90.0),
        "p_euk_2": _blast_hit("p_euk_2", "EUK", "Eustigmatales", identity=90.0),
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_1",
        contig_proteins=["p_viral", "p_euk_1", "p_euk_2"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "viral_bearing"
    assert result["n_viral_proteins"] == 1


def test_cellular_coherent_contig_with_consistent_order(tmp_path: Path) -> None:
    """Contig with 4 markers all resolving to EUK__Nannochloropsis passes."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_B"),
        _tree_nn("M3", "p3", "EUK__Nannochloropsis_A"),
        _tree_nn("M4", "p4", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=85.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_cell",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "cellular_coherent"
    assert result["n_cellular_proteins"] == 4
    # All four proteins agree on (EUK, Eustigmatales).
    assert result["cellular_lineage_purity"] == pytest.approx(1.0)


def test_cellular_scatter_fails_purity(tmp_path: Path) -> None:
    """Four EUK markers scattered across Nannochloropsis / Acanthocystis /
    NA / Other -> purity 0.25, classifier returns ambiguous."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Acanthocystis_A"),
        _tree_nn("M3", "p3", "EUK__NA_A"),
        _tree_nn("M4", "p4", "EUK__Other_A"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "various", identity=85.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_scatter",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"
    assert result["cellular_lineage_purity"] < 0.6


def test_below_marker_floor_is_ambiguous(tmp_path: Path) -> None:
    """Two cellular markers fall under the 3-marker floor -> ambiguous."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=90.0)
        for pid in ("p1", "p2")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_short",
        contig_proteins=["p1", "p2"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"


def test_low_identity_rejects_cellular_coherent(tmp_path: Path) -> None:
    """Four markers all tree-resolving to one EUK order but median identity
    45% -> HGT divergent, classifier must refuse the cellular call."""
    scorer = _make_scorer(tmp_path)
    tree = _merge(
        _tree_nn("M1", "p1", "EUK__Nannochloropsis_A"),
        _tree_nn("M2", "p2", "EUK__Nannochloropsis_A"),
        _tree_nn("M3", "p3", "EUK__Nannochloropsis_B"),
        _tree_nn("M4", "p4", "EUK__Nannochloropsis_B"),
    )
    best_hits = {
        pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=45.0)
        for pid in ("p1", "p2", "p3", "p4")
    }
    result = scorer._classify_contig_purity(
        contig_id="contig_divergent",
        contig_proteins=["p1", "p2", "p3", "p4"],
        tree_nn_results=tree,
        best_hits=best_hits,
    )
    assert result["classification"] == "ambiguous"


def test_protein_fraction_aggregates_correctly(tmp_path: Path) -> None:
    """Full bin-level aggregation: 3 contigs, one cellular_coherent, two
    viral_bearing. Check protein fraction over marker-bearing proteins."""
    scorer = _make_scorer(tmp_path)

    protein_to_contig = {
        # cellular contig (4 proteins, all EUK Nannochloropsis)
        "pC1": "contig_cell",
        "pC2": "contig_cell",
        "pC3": "contig_cell",
        "pC4": "contig_cell",
        # viral contig 1 (3 proteins)
        "pV1_a": "contig_viral_1",
        "pV1_b": "contig_viral_1",
        "pV1_c": "contig_viral_1",
        # viral contig 2 (3 proteins)
        "pV2_a": "contig_viral_2",
        "pV2_b": "contig_viral_2",
        "pV2_c": "contig_viral_2",
    }
    contig_lengths = {
        "contig_cell": 12000,
        "contig_viral_1": 30000,
        "contig_viral_2": 20000,
    }
    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", pid, "EUK__Nannochloropsis_A")
            for i, pid in enumerate(("pC1", "pC2", "pC3", "pC4"))
        ],
        *[
            _tree_nn(f"M_viral_{i}", pid, "NCLDV__Mimi_1")
            for i, pid in enumerate(
                ("pV1_a", "pV1_b", "pV1_c", "pV2_a", "pV2_b", "pV2_c")
            )
        ],
    )
    best_hits = {
        **{
            pid: _blast_hit(pid, "EUK", "Eustigmatales", identity=85.0)
            for pid in ("pC1", "pC2", "pC3", "pC4")
        },
        **{
            pid: _blast_hit(pid, "NCLDV", "Imitervirales", identity=90.0)
            for pid in (
                "pV1_a",
                "pV1_b",
                "pV1_c",
                "pV2_a",
                "pV2_b",
                "pV2_c",
            )
        },
    }

    out = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths=contig_lengths,
        tree_nn_results=tree,
        best_hits=best_hits,
        attribution_mode="fna_gene_calling",
    )
    assert out["per_contig_classifier_active"] is True
    assert out["cellular_coherent_contig_count"] == 1
    assert out["viral_bearing_contig_count"] == 2
    # 4 of 10 marker-bearing proteins on the cellular coherent contig.
    assert out["cellular_coherent_protein_fraction"] == pytest.approx(40.0)


def test_bp_fraction_aggregates_correctly(tmp_path: Path) -> None:
    """Diagnostic bp fraction: 12kb cellular of 62kb total -> ~19.35%."""
    scorer = _make_scorer(tmp_path)
    protein_to_contig = {
        "pC1": "contig_cell",
        "pC2": "contig_cell",
        "pC3": "contig_cell",
        "pV1_a": "contig_viral_1",
        "pV1_b": "contig_viral_1",
    }
    contig_lengths = {"contig_cell": 12000, "contig_viral_1": 50000}
    tree = _merge(
        *[
            _tree_nn(f"M_cell_{i}", pid, "EUK__Nannochloropsis_A")
            for i, pid in enumerate(("pC1", "pC2", "pC3"))
        ],
        _tree_nn("M_viral_0", "pV1_a", "NCLDV__Mimi_1"),
        _tree_nn("M_viral_1", "pV1_b", "NCLDV__Mimi_1"),
    )
    best_hits = {
        **{pid: _blast_hit(pid, "EUK", "Eustigmatales", 85.0)
           for pid in ("pC1", "pC2", "pC3")},
        **{pid: _blast_hit(pid, "NCLDV", "Imitervirales", 90.0)
           for pid in ("pV1_a", "pV1_b")},
    }
    out = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths=contig_lengths,
        tree_nn_results=tree,
        best_hits=best_hits,
        attribution_mode="fna_gene_calling",
    )
    expected_fraction = round((12000 / 62000) * 100.0, 2)
    assert out["cellular_coherent_bp_fraction"] == pytest.approx(expected_fraction, abs=0.01)
