"""Regression tests for contamination runtime behavior and summary output."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from unittest.mock import Mock

import pytest

from src.bin.gvclass_cli import CliOutput, combine_summary_files
from src.core.summarize_full import FullSummarizer
from src.pipeline.summary_writer import FINAL_SUMMARY_COLUMNS


def _write_minimal_database(db_path: Path) -> None:
    (db_path / "gvclassFeb26_labels.tsv").write_text("")
    (db_path / "order_completeness.tab").write_text(
        "Order\tOrthogroups\tAverage_Percent\tStd_Percent\n"
    )


@pytest.fixture
def temp_database_dir(tmp_path: Path) -> Path:
    db_path = tmp_path / "db"
    db_path.mkdir()
    _write_minimal_database(db_path)
    return db_path


def test_add_contamination_metrics_requires_trained_model(
    temp_database_dir: Path, tmp_path: Path
) -> None:
    summarizer = FullSummarizer(temp_database_dir)
    summarizer.contamination_scorer.score_rule_based = Mock(
        return_value={
            "contamination_score_v1": 12.5,
            "contamination_flag_v1": "low",
            "contamination_source_v1": "cellular",
            "contamination_cellular_signal_v1": 25.0,
            "contamination_phage_signal_v1": 0.0,
            "contamination_duplication_signal_v1": 0.0,
            "contamination_viral_mixture_signal_v1": 0.0,
            "contamination_nonviral_hit_fraction_v1": 10.0,
            "estimated_contamination": 12.5,
            "estimated_contamination_strategy": "rule_based_v1",
        }
    )
    summarizer.contamination_scorer.collect_contig_features = Mock(
        return_value={"suspicious_bp_fraction": 0.0, "suspicious_contig_count": 0}
    )
    summarizer.contamination_scorer.ml_available = False

    with pytest.raises(RuntimeError, match="trained contamination model"):
        summarizer._add_contamination_metrics(
            result={},
            query_id="query1",
            query_output_dir=tmp_path / "query1",
            marker_counts={},
            tax_counters={"order": Counter(), "family": Counter()},
        )


def test_add_contamination_metrics_uses_ml_estimate_and_query_fna_path(
    temp_database_dir: Path, tmp_path: Path
) -> None:
    summarizer = FullSummarizer(temp_database_dir)
    query_output_dir = tmp_path / "query1"
    (query_output_dir / "query_fna").mkdir(parents=True)
    (query_output_dir / "query_faa").mkdir(parents=True)
    query_fna = query_output_dir / "query_fna" / "query1.fna"
    query_faa = query_output_dir / "query_faa" / "query1.faa"
    query_fna.write_text(">contig1\nATGC\n")
    query_faa.write_text(">query1|contig1_1\nM\n")

    summarizer.contamination_scorer.score_rule_based = Mock(
        return_value={
            "contamination_score_v1": 18.0,
            "contamination_flag_v1": "low",
            "contamination_source_v1": "cellular",
            "contamination_cellular_signal_v1": 42.0,
            "contamination_phage_signal_v1": 0.0,
            "contamination_duplication_signal_v1": 0.0,
            "contamination_viral_mixture_signal_v1": 0.0,
            "contamination_nonviral_hit_fraction_v1": 12.0,
            "estimated_contamination": 18.0,
            "estimated_contamination_strategy": "rule_based_v1",
        }
    )
    summarizer.contamination_scorer.collect_contig_features = Mock(
        return_value={"suspicious_bp_fraction": 33.3, "suspicious_contig_count": 2}
    )
    summarizer.contamination_scorer.predict_contamination = Mock(
        return_value={
            "estimated_contamination": 4.25,
            "estimated_contamination_strategy": "hist_gbm_v1",
        }
    )
    summarizer.contamination_scorer.ml_available = True

    result: dict[str, object] = {}
    summarizer._add_contamination_metrics(
        result=result,
        query_id="query1",
        query_output_dir=query_output_dir,
        marker_counts={},
        tax_counters={"order": Counter(), "family": Counter()},
    )

    collect_args = summarizer.contamination_scorer.collect_contig_features.call_args[0]
    assert collect_args[0] == query_fna
    assert collect_args[1] == query_faa
    assert result["contamination_score_v1"] == 18.0
    assert result["estimated_contamination"] == 4.25
    assert result["estimated_contamination_strategy"] == "hist_gbm_v1"
    assert result["suspicious_bp_fraction_v2"] == 33.3
    assert result["suspicious_contig_count_v2"] == 2


def test_combine_summary_files_preserves_existing_final_summary(tmp_path: Path) -> None:
    summary_tsv = tmp_path / "gvclass_summary.tsv"
    summary_csv = tmp_path / "gvclass_summary.csv"
    expected_tsv = "query\ttaxonomy_strict\testimated_contamination\nq1\td_NCLDV\t4.25\n"
    expected_csv = "query,taxonomy_strict,estimated_contamination\nq1,d_NCLDV,4.25\n"
    summary_tsv.write_text(expected_tsv)
    summary_csv.write_text(expected_csv)
    (tmp_path / "q1.summary.tab").write_text(
        "query\ttaxonomy_majority\nq1\td_NCLDV\n"
    )

    combine_summary_files(tmp_path, CliOutput(plain_output=True))

    assert summary_tsv.read_text() == expected_tsv
    assert summary_csv.read_text() == expected_csv


def test_final_summary_columns_surface_single_qc_estimates() -> None:
    assert "estimated_completeness" in FINAL_SUMMARY_COLUMNS
    assert "estimated_contamination" in FINAL_SUMMARY_COLUMNS

    for column in [
        "order_completeness",
        "order_completeness_raw",
        "order_completeness_v2",
        "estimated_completeness_strategy",
        "order_weighted_completeness",
        "order_weighted_completeness_raw",
        "order_confidence_score",
        "contamination_score_v1",
        "contamination_flag_v1",
        "contamination_source_v1",
        "estimated_contamination_strategy",
        "suspicious_bp_fraction_v2",
        "suspicious_contig_count_v2",
        "cellular_unique",
        "cellular_total",
    ]:
        assert column not in FINAL_SUMMARY_COLUMNS

    for column in ["order_dup", "gvog8_dup", "vp_df", "mirus_df", "cellular_dup"]:
        assert column in FINAL_SUMMARY_COLUMNS


def test_classify_contamination_type_uses_threshold_and_source(
    temp_database_dir: Path,
) -> None:
    summarizer = FullSummarizer(temp_database_dir)
    summarizer.contamination_scorer.ml_threshold = 5.0

    assert summarizer._classify_contamination_type(
        {"estimated_contamination": 2.0, "contamination_source_v1": "cellular"}
    ) == "clean"
    assert summarizer._classify_contamination_type(
        {"estimated_contamination": 12.0, "contamination_source_v1": "cellular"}
    ) == "cellular"
    assert summarizer._classify_contamination_type(
        {"estimated_contamination": 12.0, "contamination_source_v1": "viral_mixture"}
    ) == "mixed_viral"


def test_write_contamination_candidates_file_gated_by_contamination_type(
    tmp_path: Path,
) -> None:
    from src.pipeline.query_processing_engine import _write_contamination_candidates_file

    query_output_dir = tmp_path / "query1"
    summary_data = {
        "contamination_type": "cellular",
        "estimated_contamination": 18.0,
        "_contamination_candidates": [
            {
                "contig_id": "contig_7",
                "candidate_type": "cellular",
                "reason": "cellular_hits",
                "length_bp": 12000,
                "cellular_marker_count": 2,
                "phage_marker_count": 0,
                "viral_marker_count": 1,
                "nonviral_fraction": 87.5,
                "foreign_viral_fraction": 0.0,
            }
        ],
    }

    output = _write_contamination_candidates_file(
        query_output_dir=query_output_dir,
        query_name="query1",
        summary_data=summary_data,
        logger=Mock(),
    )

    assert output is not None
    assert output.exists()
    text = output.read_text()
    assert "contig_7" in text
    assert "cellular" in text
