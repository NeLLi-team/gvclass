"""Order-level completeness metrics and the summary of a query without evidence."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from src.core.novelty_completeness import NoveltyAwareCompletenessScorer
from src.core.summarize import order_completeness, tax_format
from src.core.summarize_full import FullSummarizer
from tests.conftest import RESOURCES_ROOT, skip_if_no_runtime_resources


def _unavailable_novelty_scorer() -> NoveltyAwareCompletenessScorer:
    scorer = NoveltyAwareCompletenessScorer.__new__(NoveltyAwareCompletenessScorer)
    scorer.available = False
    return scorer


def test_mirus_taxonomy_uses_broad_mirus_quality_group(tmp_path: Path) -> None:
    accepted_lineage = tax_format.parse_taxonomy_majority(
        "d_MIRUS;p_Mirusviricota;c_Class_04;o_Demutovirales;f_Demuto_05;g_;s_"
    )
    completeness_table = tmp_path / "order_completeness.tab"
    completeness_table.write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\nMIRUS\t50\t5\tOG1, OG2\n"
    )

    class _Weighted:
        def calculate_weighted_completeness(self, **_kwargs):
            return 50.0, 80.0, {}

    result = {}
    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        _Weighted(),
        _unavailable_novelty_scorer(),
        "novelty-aware",
        result,
        {"OG1": 1},
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
    result = {}

    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        SimpleNamespace(),
        _unavailable_novelty_scorer(),
        "novelty-aware",
        result,
        {},
        tax_format.parse_taxonomy_majority(
            "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;"
            "o_Asfuvirales;f_Asfarviridae;g_;s_"
        ),
    )

    assert math.isnan(result["estimated_completeness"])
    assert result["completeness_query_support"] == "insufficient_evidence"


class _Weighted:
    """Weighted calculator that finds half of the expected markers."""

    def calculate_weighted_completeness(self, **_kwargs: object) -> tuple:
        return 50.0, 80.0, {}


def _order_metrics(tmp_path: Path, baseline: str) -> dict:
    """Return the metrics of a query with OG1 twice and OG2 missing."""
    completeness_table = tmp_path / "order_completeness.tab"
    completeness_table.write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\n"
        f"Asfuvirales\t{baseline}\t5\tOG1, OG2\n"
    )
    result: dict = {}
    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        _Weighted(),
        SimpleNamespace(),
        "novelty-aware",
        result,
        {"OG1": 2},
        tax_format.parse_taxonomy_majority("d_NCLDV;o_Asfuvirales"),
    )
    return result


@pytest.mark.parametrize("baseline", ["nan", "inf"])
def test_non_finite_baseline_makes_normalized_estimates_unavailable(
    tmp_path: Path, baseline: str
) -> None:
    result = _order_metrics(tmp_path, baseline)

    assert math.isnan(result["order_completeness"])
    assert math.isnan(result["order_weighted_completeness"])
    assert math.isnan(result["estimated_completeness"])
    assert result["completeness_query_support"] == "unavailable"


@pytest.mark.parametrize("baseline", ["nan", "inf"])
def test_non_finite_baseline_keeps_raw_recovery_duplication_and_weighted_value(
    tmp_path: Path, baseline: str
) -> None:
    result = _order_metrics(tmp_path, baseline)

    assert result["order_completeness_raw"] == 50.0
    assert result["order_dup"] == 2.0
    assert result["order_weighted_completeness_raw"] == 50.0
    assert result["order_confidence_score"] == 80.0
    assert result["order_completeness_reference_order"] == "Asfuvirales"


@pytest.mark.parametrize("mode", ["novelty-aware", "legacy"])
def test_unavailable_novelty_resources_still_fill_every_v2_key(
    tmp_path: Path, mode: str
) -> None:
    completeness_table = tmp_path / "order_completeness.tab"
    completeness_table.write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\n"
        "Asfuvirales\t50\t5\tOG1, OG2\n"
    )

    class _Weighted:
        def calculate_weighted_completeness(self, **_kwargs):
            return 50.0, 80.0, {}

    result = {}

    order_completeness._add_order_metrics(
        completeness_table,
        pd.read_csv(completeness_table, sep="\t"),
        _Weighted(),
        _unavailable_novelty_scorer(),
        mode,
        result,
        {"OG1": 1},
        tax_format.parse_taxonomy_majority("d_NCLDV;o_Asfuvirales"),
    )

    v2_keys = {key for key in result if key.startswith("order_completeness_v2")}
    assert v2_keys == {
        "order_completeness_v2",
        "order_completeness_v2_strategy",
        "order_completeness_v2_strategy2_raw",
        "order_completeness_v2_strategy2_normalized",
        "order_completeness_v2_support_score",
        "order_completeness_v2_ood_flag",
        "order_completeness_v2_reference_group",
        "order_completeness_v2_validation_mode",
        "order_completeness_v2_informative_fraction",
    }
    assert math.isnan(result["order_completeness_v2"])
    assert result["order_completeness_v2_strategy"] == "unavailable"
    assert result["estimated_completeness"] == 100.0


@skip_if_no_runtime_resources()
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
