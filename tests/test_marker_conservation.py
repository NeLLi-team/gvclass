"""M12: weighted completeness must load real per-order marker conservation.

Previously create_weighted_calculator pointed at markers/stats.tsv (an HMM size
table), so marker_stats_df was always None and every weight defaulted to 1.0,
making weighted_order_completeness identical to the unweighted value.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from tests.conftest import skip_if_no_runtime_resources

REPO_ROOT = Path(__file__).resolve().parents[1]
RESOURCES = REPO_ROOT / "resources"
CONSERVATION = RESOURCES / "markers" / "marker_conservation.tsv"
TIERS = RESOURCES / "completeness" / "tiers.tsv"


@skip_if_no_runtime_resources()
def test_marker_conservation_file_present_and_schema():
    assert CONSERVATION.exists(), "marker_conservation.tsv not generated"
    with CONSERVATION.open() as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    for col in ("marker_name", "order_name", "percent_genomes_with_marker"):
        assert col in header


@skip_if_no_runtime_resources()
def test_calculator_loads_conservation_stats():
    from src.core.weighted_completeness import create_weighted_calculator

    calc = create_weighted_calculator(RESOURCES)
    assert calc.marker_stats_df is not None
    for col in ("marker_name", "order_name", "percent_genomes_with_marker"):
        assert col in calc.marker_stats_df.columns


@skip_if_no_runtime_resources()
def test_weights_are_not_uniform():
    from src.core.weighted_completeness import create_weighted_calculator

    calc = create_weighted_calculator(RESOURCES)
    df = calc.marker_stats_df
    assert df is not None
    order = "Algavirales"
    markers = tuple(df[df["order_name"] == order]["marker_name"].tolist())
    assert markers, f"no markers for {order} in conservation table"
    weights = calc.get_marker_weights(order, markers)
    # The whole point: not the uniform 1.0 fallback.
    distinct = {round(w, 6) for w in weights.values()}
    assert distinct != {1.0}, "weights are uniform -> conservation not applied"


@skip_if_no_runtime_resources()
def test_order_core_markers_have_high_prevalence():
    # Every ORDER-level core-tier marker must have prevalence >= core_prevalence
    # (0.9 -> 90%) in the conservation table they are both derived from.
    conservation: dict[tuple[str, str], float] = {}
    with CONSERVATION.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            conservation[(row["order_name"], row["marker_name"])] = float(
                row["percent_genomes_with_marker"]
            )

    checked = 0
    with TIERS.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["group_level"] != "order" or row["tier"] != "core":
                continue
            key = (row["order_name"], row["marker_name"])
            assert key in conservation, f"core marker missing from conservation: {key}"
            assert conservation[key] >= 90.0 - 1e-6, (
                f"core marker {key} prevalence {conservation[key]} < 90%"
            )
            checked += 1
    assert checked > 0, "no order-level core markers found to check"
