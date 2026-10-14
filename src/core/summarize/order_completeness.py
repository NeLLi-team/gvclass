"""
Order-level completeness metrics for the GVClass summary.

Stateless helpers that score a query against its assigned taxonomic order.

The reference tables and scorers these functions need
(``completeness_table``, ``order_stats_df``, ``weighted_calculator``,
``novelty_scorer``, ``completeness_mode``) are passed in explicitly, so nothing
here reads instance state.
"""

import logging
import math
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
    quality_group: str,
    family_tax: str = "",
) -> Dict[str, Any]:
    """Calculate order-specific completeness and duplication metrics.

    Returns:
        Dictionary of raw and normalized order completeness metrics.
    """
    try:
        order_ogs = _load_order_orthogroups(
            completeness_table, quality_group, order_stats_df=order_stats_df
        )
        if not order_ogs:
            return _default_order_metrics(
                quality_group,
                query_support="unsupported_target",
                model_reliability="unavailable",
            )

        marker_counts = _load_marker_counts(counts_file)
        if not any(marker_counts.get(og, 0) > 0 for og in order_ogs):
            return _default_order_metrics(
                quality_group,
                query_support="insufficient_evidence",
                model_reliability="unvalidated_legacy",
            )
        raw_completeness, duplication = _calculate_traditional_order_metrics(
            marker_counts, order_ogs
        )
        raw_weighted_completeness, confidence_score = _calculate_weighted_order_metrics(
            weighted_calculator,
            marker_counts,
            quality_group,
            order_ogs,
            raw_completeness,
        )
        normalized_completeness, baseline_mean, baseline_std = (
            _normalize_order_completeness(
                order_stats_df, raw_completeness, quality_group
            )
        )
        normalized_weighted_completeness, _, _ = _normalize_order_completeness(
            order_stats_df, raw_weighted_completeness, quality_group
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
            "order_completeness_reference_order": quality_group,
            "order_completeness_strategy": "order_baseline_ratio_v1",
            "estimated_completeness": normalized_completeness,
            "estimated_completeness_strategy": "order_baseline_ratio_v1",
            "completeness_query_support": "supported",
            "completeness_model_group": quality_group,
            "completeness_model_reliability": "unvalidated_legacy",
        }
        if novelty_scorer.available:
            result.update(
                novelty_scorer.calculate(
                    marker_counts=marker_counts,
                    order_tax=quality_group,
                    family_tax=family_tax,
                    strategy1_raw=raw_completeness,
                    strategy1_score=normalized_completeness,
                    selected_mode=completeness_mode,
                )
            )
        return result

    except Exception as e:
        logger.error(f"Error calculating order metrics: {e}")
        return _default_order_metrics(
            quality_group,
            query_support="unavailable",
            model_reliability="unavailable",
        )


def _default_order_metrics(
    quality_group: str = "",
    query_support: str = "unassigned",
    model_reliability: str = "unavailable",
) -> Dict[str, Any]:
    """Return default order completeness metrics."""
    resource_group = (
        quality_group
        if query_support not in {"unassigned", "unsupported_target"}
        else ""
    )
    return {
        "order_completeness": math.nan,
        "order_completeness_raw": math.nan,
        "order_dup": math.nan,
        "order_weighted_completeness": math.nan,
        "order_weighted_completeness_raw": math.nan,
        "order_confidence_score": 0.0,
        "order_completeness_baseline_mean": math.nan,
        "order_completeness_baseline_std": math.nan,
        "order_completeness_reference_order": quality_group or "unavailable",
        "order_completeness_strategy": "unavailable",
        "order_completeness_v2": math.nan,
        "order_completeness_v2_strategy": "unavailable",
        "order_completeness_v2_strategy2_raw": math.nan,
        "order_completeness_v2_strategy2_normalized": math.nan,
        "order_completeness_v2_support_score": 0.0,
        "order_completeness_v2_ood_flag": query_support,
        "order_completeness_v2_reference_group": (
            f"order:{resource_group}" if resource_group else "unavailable"
        ),
        "order_completeness_v2_validation_mode": model_reliability,
        "order_completeness_v2_informative_fraction": 0.0,
        "estimated_completeness": math.nan,
        "estimated_completeness_strategy": "unavailable",
        "completeness_query_support": query_support,
        "completeness_model_group": resource_group or "unavailable",
        "completeness_model_reliability": model_reliability,
    }


def _load_order_orthogroups(
    completeness_table: Path,
    order_tax: str,
    order_stats_df: pd.DataFrame | None = None,
) -> List[str]:
    """Load expected orthogroups for a specific taxonomic order."""
    if (
        order_stats_df is not None
        and not order_stats_df.empty
        and "Order" in order_stats_df.columns
        and "Orthogroups" in order_stats_df.columns
    ):
        comp_df = order_stats_df
    else:
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


def _select_quality_group(accepted_lineage: Dict[str, str]) -> str:
    """Map an accepted lineage to the matching completeness resource group."""
    if (
        accepted_lineage.get("domain") == "MIRUS"
        or accepted_lineage.get("phylum") == "Mirusviricota"
    ):
        return "MIRUS"
    return accepted_lineage.get("order", "")


def _add_order_metrics(
    completeness_table: Path,
    order_stats_df: pd.DataFrame,
    weighted_calculator: WeightedCompletenessCalculator,
    novelty_scorer: NoveltyAwareCompletenessScorer,
    completeness_mode: str,
    result: Dict[str, Any],
    counts_file: Path,
    accepted_lineage: Dict[str, str],
) -> None:
    """Populate order-level completeness metrics from order assignment."""
    quality_group = _select_quality_group(accepted_lineage)
    family_tax = accepted_lineage.get("family", "")
    if quality_group:
        result.update(
            calculate_order_metrics(
                completeness_table,
                order_stats_df,
                weighted_calculator,
                novelty_scorer,
                completeness_mode,
                counts_file,
                quality_group,
                family_tax,
            )
        )
        return

    result.update(_default_order_metrics())
