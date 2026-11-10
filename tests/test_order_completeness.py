"""Order-level completeness metrics when the reference baseline is not finite."""

from __future__ import annotations

import math
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from src.core.summarize import order_completeness, tax_format


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
