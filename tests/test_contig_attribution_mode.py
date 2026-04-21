"""Unit tests for the v1.4.3 Phase 2 ``contig_attribution_mode`` fallback paths.

Covers the six attribution-mode cases listed in
docs/plans/active/v1_4_3_contig_purity_contamination.md under the
"Weak-attribution fallback" test matrix.
"""

from __future__ import annotations

import logging
from pathlib import Path


def _make_scorer(tmp_path: Path):
    from src.core.contamination_scoring import ContaminationScorer
    from tests.conftest import stage_db_resources

    db = tmp_path / "db"
    stage_db_resources(db, markers=False)
    (db / "labels.tsv").write_text("")
    (db / "order_completeness.tab").write_text(
        "Order\tOrthogroups\tAverage_Percent\tStd_Percent\n"
    )
    return ContaminationScorer(db)


def test_fna_gene_calling_mode_runs_per_contig_classifier(tmp_path: Path) -> None:
    """`.fna` runs trust contig attribution unconditionally."""
    scorer = _make_scorer(tmp_path)
    protein_to_contig = {f"p_{i}": "contig_A" for i in range(10)}
    mode = scorer._detect_contig_attribution_mode("fna", protein_to_contig)
    assert mode == "fna_gene_calling"


def test_faa_contig_encoded_mode_runs_per_contig_classifier(tmp_path: Path) -> None:
    """`.faa` runs with IDs encoding the contig ID keep per-contig mode."""
    scorer = _make_scorer(tmp_path)
    # 10 proteins on 2 contigs -> 2 / 10 = 0.2 << 0.8 threshold
    protein_to_contig = {f"p_{i}": "contig_X" for i in range(5)}
    protein_to_contig.update({f"q_{i}": "contig_Y" for i in range(5)})
    mode = scorer._detect_contig_attribution_mode("faa", protein_to_contig)
    assert mode == "faa_contig_encoded"


def test_faa_weak_per_protein_fallback_emits_neutral_features(tmp_path: Path) -> None:
    """20 unique pseudo-contigs (one-protein-each) triggers the weak regime."""
    scorer = _make_scorer(tmp_path)
    protein_to_contig = {f"p_{i}": f"pseudo_{i}" for i in range(20)}
    mode = scorer._detect_contig_attribution_mode("faa", protein_to_contig)
    assert mode == "faa_weak_per_protein"

    features = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths={cid: 500 for cid in protein_to_contig.values()},
        tree_nn_results={},
        best_hits={},
        attribution_mode=mode,
    )
    assert features["per_contig_classifier_active"] is False
    assert features["cellular_coherent_contig_count"] == 0
    assert features["cellular_coherent_protein_fraction"] == 0.0
    assert features["cellular_coherent_bp_fraction"] == 0.0
    assert features["cellular_lineage_purity_median"] == 0.0
    assert features["cellular_hit_identity_median"] == 0.0
    assert features["viral_bearing_contig_count"] == 0


def test_faa_weak_per_protein_fallback_preserves_legacy_cellular_trigger(
    tmp_path: Path,
) -> None:
    """With neutral per-contig features AND a non-zero rule-based cellular
    source, the legacy `cellular` label still fires from
    ``_classify_contamination_type``."""
    from src.core.summarize_full import FullSummarizer
    from tests.conftest import stage_db_resources

    stage_db_resources(tmp_path, markers=False)
    summarizer = FullSummarizer(tmp_path)
    summarizer.contamination_scorer.ml_threshold = 5.0

    result = {
        "estimated_contamination": 18.0,
        "contamination_source_v1": "cellular",
        "cellular_coherent_contig_count": 0,
    }
    label = summarizer._classify_contamination_type(result)
    assert label == "cellular"


def test_low_data_regime_below_5_proteins_uses_legacy_rule(tmp_path: Path) -> None:
    """< 5 marker-bearing proteins -> low-data regime -> neutral features
    regardless of attribution mode."""
    scorer = _make_scorer(tmp_path)
    protein_to_contig = {"p_1": "c1", "p_2": "c1", "p_3": "c2"}  # only 3
    mode = scorer._detect_contig_attribution_mode("fna", protein_to_contig)
    assert mode == "fna_gene_calling"
    features = scorer.compute_contig_purity_features(
        protein_to_contig=protein_to_contig,
        contig_lengths={"c1": 10000, "c2": 5000},
        tree_nn_results={},
        best_hits={},
        attribution_mode=mode,
    )
    assert features["per_contig_classifier_active"] is False
    assert features["cellular_coherent_contig_count"] == 0


def test_novel_virus_downgrade_to_uncertain(tmp_path: Path) -> None:
    """Per v1.4.3 Phase 2: when the per-contig classifier finds zero
    cellular_coherent contigs AND the bin carries ≥3 viral_bearing
    contigs AND the rule-based source is viral_mixture, the label
    downgrades from ``mixed_viral`` to ``uncertain`` — this is the
    novel-virus-with-HGT-scatter pattern, not real multi-viral mixing."""
    from src.core.summarize_full import FullSummarizer
    from tests.conftest import stage_db_resources

    stage_db_resources(tmp_path, markers=False)
    summarizer = FullSummarizer(tmp_path)
    summarizer.contamination_scorer.ml_threshold = 5.0

    result = {
        "estimated_contamination": 27.0,
        "contamination_source_v1": "viral_mixture",
        "cellular_coherent_contig_count": 0,
        "viral_bearing_contig_count": 3,
    }
    assert summarizer._classify_contamination_type(result) == "uncertain"


def test_mixed_viral_preserved_when_multi_giant_virus_with_no_viral_bearing_contigs(
    tmp_path: Path,
) -> None:
    """The downgrade ONLY fires when viral_bearing >= 3. A bin that looks
    multi-viral-order but has too few viral_bearing contigs keeps the
    legacy ``mixed_viral`` label."""
    from src.core.summarize_full import FullSummarizer
    from tests.conftest import stage_db_resources

    stage_db_resources(tmp_path, markers=False)
    summarizer = FullSummarizer(tmp_path)
    summarizer.contamination_scorer.ml_threshold = 5.0

    result = {
        "estimated_contamination": 27.0,
        "contamination_source_v1": "viral_mixture",
        "cellular_coherent_contig_count": 0,
        "viral_bearing_contig_count": 0,
    }
    assert summarizer._classify_contamination_type(result) == "mixed_viral"


def test_weak_attribution_emits_info_log(tmp_path: Path, caplog) -> None:
    """The scorer must log an INFO record naming the fallback path when
    contig_attribution_mode resolves to faa_weak_per_protein."""
    scorer = _make_scorer(tmp_path)

    query_faa = tmp_path / "weak.faa"
    # `protein_to_contig_id` strips the trailing "_<N>" token, so we use
    # IDs where the portion before the trailing "_X" is unique per
    # protein — that way every protein maps to its own pseudo-contig and
    # the weak-attribution heuristic fires.
    query_faa.write_text(
        "\n".join(f">stateless{i}_A\nMQ" for i in range(10)) + "\n"
    )
    blast_dir = tmp_path / "blastp_out"
    blast_dir.mkdir()

    caplog.set_level(logging.INFO, logger="src.core.contamination_scoring")

    scorer.collect_contig_features(
        query_fna=Path("/nonexistent.fna"),
        query_faa=query_faa,
        blast_dir=blast_dir,
        primary_order="",
        primary_family="",
        tree_nn_results={},
    )
    messages = [rec.getMessage() for rec in caplog.records if rec.levelno == logging.INFO]
    assert any(
        "faa_weak_per_protein" in msg.lower() or "weak contig attribution" in msg.lower()
        for msg in messages
    ), messages
