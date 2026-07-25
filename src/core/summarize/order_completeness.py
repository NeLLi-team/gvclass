"""
Order-level completeness metrics for the GVClass summary.

Stateless counterparts of the ``FullSummarizer`` methods that score a query
against its assigned taxonomic order; the class keeps thin delegating wrappers
for every name here.

The reference tables and scorers these functions need
(``completeness_table``, ``order_stats_df``, ``weighted_calculator``,
``novelty_scorer``, ``completeness_mode``) are passed in explicitly, so nothing
here reads instance state.
"""

from collections import Counter
import logging
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd

from src.core.novelty_completeness import NoveltyAwareCompletenessScorer
from src.core.summarize.paths_stats import _load_marker_counts
from src.core.weighted_completeness import WeightedCompletenessCalculator

logger = logging.getLogger(__name__)


def calculate_order_metrics(
    completeness_table: Path,
    order_stats_df: pd.DataFrame,
    weighted_calculator: WeightedCompletenessCalculator,
    novelty_scorer: NoveltyAwareCompletenessScorer,
    completeness_mode: str,
    counts_file: Path,
    order_tax: str,
    family_tax: str = "",
) -> Dict[str, Any]:
    """Calculate order-specific completeness and duplication metrics.

    Returns:
        Dictionary of raw and normalized order completeness metrics.
    """
    try:
        order_ogs = _load_order_orthogroups(completeness_table, order_tax)
        if not order_ogs:
            return _default_order_metrics(order_tax)

        marker_counts = _load_marker_counts(counts_file)
        raw_completeness, duplication = _calculate_traditional_order_metrics(
            marker_counts, order_ogs
        )
        raw_weighted_completeness, confidence_score = _calculate_weighted_order_metrics(
            weighted_calculator, marker_counts, order_tax, order_ogs, raw_completeness
        )
        normalized_completeness, baseline_mean, baseline_std = (
            _normalize_order_completeness(order_stats_df, raw_completeness, order_tax)
        )
        normalized_weighted_completeness, _, _ = _normalize_order_completeness(
            order_stats_df, raw_weighted_completeness, order_tax
        )
        result = {
            "order_completeness": normalized_completeness,
            "order_completeness_raw": raw_completeness,
            "order_dup": duplication,
            "order_weighted_completeness": normalized_weighted_completeness,
            "order_weighted_completeness_raw": raw_weighted_completeness,
            "order_confidence_score": confidence_score,
            "order_completeness_baseline_mean": baseline_mean,
            "order_completeness_baseline_std": baseline_std,
            "order_completeness_reference_order": order_tax,
            "order_completeness_strategy": "order_baseline_ratio_v1",
            "estimated_completeness": normalized_completeness,
            "estimated_completeness_strategy": "order_baseline_ratio_v1",
        }
        if novelty_scorer.available:
            result.update(
                novelty_scorer.calculate(
                    marker_counts=marker_counts,
                    order_tax=order_tax,
                    family_tax=family_tax,
                    strategy1_raw=raw_completeness,
                    strategy1_score=normalized_completeness,
                    selected_mode=completeness_mode,
                )
            )
        return result

    except Exception as e:
        logger.error(f"Error calculating order metrics: {e}")
        return _default_order_metrics(order_tax)


def _default_order_metrics(order_tax: str = "") -> Dict[str, Any]:
    """Return default order completeness metrics."""
    return {
        "order_completeness": 0.0,
        "order_completeness_raw": 0.0,
        "order_dup": 0.0,
        "order_weighted_completeness": 0.0,
        "order_weighted_completeness_raw": 0.0,
        "order_confidence_score": 0.0,
        "order_completeness_baseline_mean": 0.0,
        "order_completeness_baseline_std": 0.0,
        "order_completeness_reference_order": order_tax,
        "order_completeness_strategy": "order_baseline_ratio_v1",
        "order_completeness_v2": 0.0,
        "order_completeness_v2_strategy": "novelty_aware_v1",
        "order_completeness_v2_strategy2_raw": 0.0,
        "order_completeness_v2_strategy2_normalized": 0.0,
        "order_completeness_v2_support_score": 0.0,
        "order_completeness_v2_ood_flag": "unassigned",
        "order_completeness_v2_reference_group": "unavailable",
        "order_completeness_v2_validation_mode": "unavailable",
        "order_completeness_v2_informative_fraction": 0.0,
        "estimated_completeness": 0.0,
        "estimated_completeness_strategy": "order_baseline_ratio_v1",
    }


def _load_order_orthogroups(completeness_table: Path, order_tax: str) -> List[str]:
    """Load expected orthogroups for a specific taxonomic order."""
    comp_df = pd.read_csv(completeness_table, sep="\t")
    order_row = comp_df[comp_df["Order"] == order_tax]
    if order_row.empty:
        return []
    ogs_str = str(order_row["Orthogroups"].values[0])
    return [og.strip() for og in ogs_str.split(",") if og.strip()]


def _calculate_traditional_order_metrics(
    marker_counts: Dict[str, int], order_ogs: List[str]
) -> Tuple[float, float]:
    """Calculate completeness and duplication from expected order markers."""
    present_ogs = sum(1 for og in order_ogs if marker_counts.get(og, 0) > 0)
    total_hits = sum(marker_counts.get(og, 0) for og in order_ogs)

    completeness = (present_ogs / len(order_ogs)) * 100 if order_ogs else 0.0
    duplication = total_hits / present_ogs if present_ogs > 0 else 0.0
    return completeness, duplication


def _calculate_weighted_order_metrics(
    weighted_calculator: WeightedCompletenessCalculator,
    marker_counts: Dict[str, int],
    order_tax: str,
    order_ogs: List[str],
    fallback_completeness: float,
) -> Tuple[float, float]:
    """Calculate ML-weighted completeness with fallback to traditional metrics."""
    try:
        weighted_completeness, confidence_score, _ = (
            weighted_calculator.calculate_weighted_completeness(
                marker_counts=marker_counts,
                taxonomic_order=order_tax,
                expected_markers=order_ogs,
            )
        )
        logger.debug(
            f"Weighted completeness for {order_tax}: {weighted_completeness:.2f}% "
            f"(traditional: {fallback_completeness:.2f}%)"
        )
        return weighted_completeness, confidence_score
    except Exception as e:
        logger.warning(f"Weighted completeness calculation failed: {e}")
        return fallback_completeness, 50.0


def _normalize_order_completeness(
    order_stats_df: pd.DataFrame, raw_completeness: float, order_tax: str
) -> Tuple[float, float, float]:
    """Normalize marker recovery against complete-reference order baselines."""
    if order_stats_df.empty:
        return raw_completeness, 0.0, 0.0

    order_row = order_stats_df[order_stats_df["Order"] == order_tax]
    if order_row.empty:
        return raw_completeness, 0.0, 0.0

    baseline_mean = float(order_row["Average_Percent"].iloc[0])
    baseline_std = float(order_row["Std_Percent"].iloc[0])
    if baseline_mean <= 0:
        return raw_completeness, baseline_mean, baseline_std

    normalized = (raw_completeness / baseline_mean) * 100.0
    normalized = max(0.0, min(100.0, normalized))
    return normalized, baseline_mean, baseline_std


def _extract_order_taxonomy(tax_counters: Dict[str, Counter]) -> str:
    """Extract most-supported order taxonomy token from counters."""
    if not tax_counters["order"]:
        return ""

    most_common_order = tax_counters["order"].most_common(1)[0][0]
    if "__" not in most_common_order:
        return ""
    parts = most_common_order.split("__")
    return parts[1] if len(parts) > 1 else most_common_order


def _extract_family_taxonomy(tax_counters: Dict[str, Counter]) -> str:
    """Extract most-supported family taxonomy token from counters."""
    if not tax_counters["family"]:
        return ""

    most_common_family = tax_counters["family"].most_common(1)[0][0]
    if "__" not in most_common_family:
        return ""
    parts = most_common_family.split("__")
    return parts[1] if len(parts) > 1 else most_common_family


def _add_order_metrics(
    completeness_table: Path,
    order_stats_df: pd.DataFrame,
    weighted_calculator: WeightedCompletenessCalculator,
    novelty_scorer: NoveltyAwareCompletenessScorer,
    completeness_mode: str,
    result: Dict[str, Any],
    counts_file: Path,
    tax_counters: Dict[str, Counter],
) -> None:
    """Populate order-level completeness metrics from order assignment."""
    order_tax = _extract_order_taxonomy(tax_counters)
    family_tax = _extract_family_taxonomy(tax_counters)
    if order_tax:
        result.update(
            calculate_order_metrics(
                completeness_table,
                order_stats_df,
                weighted_calculator,
                novelty_scorer,
                completeness_mode,
                counts_file,
                order_tax,
                family_tax,
            )
        )
        return

    result.update(_default_order_metrics(order_tax))
