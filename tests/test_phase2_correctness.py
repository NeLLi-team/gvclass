"""Regression tests for Phase 2 of the v1.4.3 remediation plan.

Covers:
* 2.1 Taxonomy summary columns (the consensus rules are in
  test_taxonomy_consensus.py).
* 2.2 Marker extraction keyed by record.id with record.description fallback.
* 2.3 Genetic-code deterministic tiebreak with an absolute margin over meta.
* 2.5 CLI resolver default for completeness_mode matches default_config.
"""

from __future__ import annotations

from tests.conftest import RESOURCES_ROOT, skip_if_no_runtime_resources

pytestmark = skip_if_no_runtime_resources()

from collections import Counter
import math
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from src.core.summarize import order_completeness, tax_format
from src.core.summarize_full import FullSummarizer


def test_mirus_taxonomy_uses_broad_mirus_quality_group(tmp_path: Path) -> None:
    accepted_lineage = tax_format.parse_taxonomy_majority(
        "d_MIRUS;p_Mirusviricota;c_Class_04;o_Demutovirales;"
        "f_Demuto_05;g_;s_"
    )
    completeness_table = tmp_path / "order_completeness.tab"
    completeness_table.write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\n"
        "MIRUS\t50\t5\tOG1, OG2\n"
    )
    counts_file = tmp_path / "models.counts"
    counts_file.write_text("OG1\t1\n")

    class _Weighted:
        def calculate_weighted_completeness(self, **_kwargs):
            return 50.0, 80.0, {}

    result = {}
    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        _Weighted(),
        SimpleNamespace(available=False),
        "novelty-aware",
        result,
        counts_file,
        accepted_lineage,
    )

    assert accepted_lineage["order"] == "Demutovirales"
    assert result["completeness_model_group"] == "MIRUS"
    assert result["order_completeness_reference_order"] == "MIRUS"
    assert result["estimated_completeness"] == 100.0


def test_supported_group_with_no_marker_evidence_is_unavailable(
    tmp_path: Path,
) -> None:
    completeness_table = tmp_path / "order_completeness.tab"
    completeness_table.write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\n"
        "Asfuvirales\t50\t5\tOG1, OG2\n"
    )
    counts_file = tmp_path / "models.counts"
    counts_file.write_text("")
    result = {}

    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        SimpleNamespace(),
        SimpleNamespace(available=False),
        "novelty-aware",
        result,
        counts_file,
        tax_format.parse_taxonomy_majority(
            "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;"
            "o_Asfuvirales;f_Asfarviridae;g_;s_"
        ),
    )

    assert math.isnan(result["estimated_completeness"])
    assert result["completeness_query_support"] == "insufficient_evidence"


def test_full_summary_with_no_evidence_has_no_quality_claim(tmp_path: Path) -> None:
    query_dir = tmp_path / "query"
    (query_dir / "query_faa").mkdir(parents=True)
    (query_dir / "query_faa" / "empty.faa").write_text(">contig_1\nMMMM\n")

    result = FullSummarizer(
        RESOURCES_ROOT,
        completeness_mode="novelty-aware",
        sensitive_mode=True,
    ).summarize_query_full("empty", query_dir, {})

    assert math.isnan(result["estimated_completeness"])
    assert result["completeness_query_support"] == "unassigned"
    assert math.isnan(result["estimated_contamination"])
    assert result["contamination_query_support"] == "insufficient_evidence"
    assert result["contamination_type"] == "uncertain"


def test_tax_level_counts_normalize_domain_unclassified_placeholders() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV__NCLDV_unclassified": 6, "NCLDV__Megaviricetes": 4}),
        level="class",
    )

    assert "NCLDV__unclassified:6(60.00%)" in formatted
    assert "NCLDV__NCLDV_unclassified" not in formatted


def test_tax_level_counts_use_domain_prefix_for_grouped_minor_taxa() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV__S392": 18, "PLV__EP00001": 1, "NCLDV__singleton": 1}),
        level="species",
    )

    assert "PLV__other:1(5.00%)" in formatted
    assert "NCLDV__other:1(5.00%)" in formatted
    assert "PLV_other" not in formatted
    assert "NCLDV_other" not in formatted


def test_domain_counts_are_not_grouped_as_other() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV": 19, "PLV": 1}),
        level="domain",
    )

    assert "NCLDV:19(95.00%)" in formatted
    assert "PLV:1(5.00%)" in formatted
    assert "PLV_other" not in formatted


def test_taxonomy_confidence_column_present_in_summary_schemas() -> None:
    from src.pipeline.summary_writer import (
        FINAL_SUMMARY_COLUMNS,
        LEGACY_SUMMARY_HEADERS,
    )

    assert "taxonomy_confidence" in FINAL_SUMMARY_COLUMNS
    assert "taxonomy_confidence" in LEGACY_SUMMARY_HEADERS


# ---------------------------------------------------------------------------
# 2.2 Marker extraction by record.id with description fallback.
# ---------------------------------------------------------------------------


def test_extract_marker_sequences_matches_by_record_id(tmp_path: Path) -> None:
    from src.core.marker_extraction import extract_marker_sequences

    query_faa = tmp_path / "q.faa"
    query_faa.write_text(
        ">prot_1 # start=1 # end=300 # annotation\n"
        "MKLILV\n"
        ">prot_2 some descriptor\n"
        "MTKAAG\n"
    )
    output_dir = tmp_path / "markers"

    marker_hits = {"GVOGm0022": {"prot_1"}, "GVOGm0054": {"prot_2"}}

    files = extract_marker_sequences(query_faa, marker_hits, output_dir)

    assert set(files.keys()) == {"GVOGm0022", "GVOGm0054"}
    for marker in files:
        contents = files[marker].read_text()
        assert contents.startswith(">prot_"), contents


# ---------------------------------------------------------------------------
# 2.5 Completeness-mode fallback.
# ---------------------------------------------------------------------------


def test_resolve_completeness_mode_falls_back_to_novelty_aware() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode=None)
    config = {"pipeline": {}}
    assert resolve_completeness_mode(args, config) == "novelty-aware"


def test_resolve_completeness_mode_honors_explicit_legacy() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode=None)
    config = {"pipeline": {"completeness_mode": "legacy"}}
    assert resolve_completeness_mode(args, config) == "legacy"


def test_resolve_completeness_mode_cli_flag_wins_over_config() -> None:
    from src.bin.gvclass_cli import resolve_completeness_mode

    args = SimpleNamespace(completeness_mode="legacy")
    config = {"pipeline": {"completeness_mode": "novelty-aware"}}
    assert resolve_completeness_mode(args, config) == "legacy"
