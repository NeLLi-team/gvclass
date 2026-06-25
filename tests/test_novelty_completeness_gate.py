"""Unit tests for novelty-aware completeness gate (H3) and fallback (M13).

These build a scorer with controlled in-memory state (bypassing the on-disk
resource bundle) so the decision branches are exercised without the DB.
"""

from __future__ import annotations

import pytest

from src.core.novelty_completeness import (
    NoveltyAwareCompletenessScorer,
    Strategy2Config,
    TierSet,
)


class _FakeModel:
    def __init__(self, value: float) -> None:
        self._value = value

    def predict(self, _X):
        return [self._value]


def _config() -> Strategy2Config:
    return Strategy2Config(
        name="cfg_test",
        core_prevalence=0.9,
        shared_prevalence=0.65,
        accessory_prevalence=0.5,
        core_max_copy=1.2,
        shared_max_copy=1.5,
        family_min_ref_genomes=3,
        family_min_informative_markers=2,
        weight_core=0.7,
        weight_shared=0.2,
        weight_accessory=0.1,
    )


def _make_scorer(*, models, metadata, tiers, baselines, available=True):
    scorer = NoveltyAwareCompletenessScorer.__new__(NoveltyAwareCompletenessScorer)
    scorer.available = available
    scorer.config = _config()
    scorer.feature_names = ["strategy1_score", "strategy2_score"]
    scorer.order_expected_marker_counts = {"Algavirales": 71, "Imitervirales": 44}
    scorer.tiers = tiers
    scorer.baselines = baselines
    scorer.model_metadata = metadata
    scorer.models = models
    return scorer


def _tiered_scorer(test_r2):
    order = "Algavirales"
    key = ("order", order, order)
    tiers = {key: TierSet(core=("OG1", "OG2"), shared=(), accessory=())}
    baselines = {key: {"baseline_mean": 40.0, "baseline_median": 40.0, "n_refs": 10}}
    metadata = {
        order: {
            "validation_type": "leave_one_family_out",
            "test_r2": test_r2,
            "test_mae": 14.0,
        }
    }
    return _make_scorer(
        models={order: _FakeModel(80.0)},
        metadata=metadata,
        tiers=tiers,
        baselines=baselines,
    )


# --- H3: R2 gate must key on the shipped `test_r2` column -------------------


def test_r2_gate_fires_for_low_r2_order():
    scorer = _tiered_scorer(test_r2=0.4898)  # < R2_ADVISORY_FLOOR (0.5)
    result = scorer.calculate(
        {"OG1": 1, "OG2": 1}, "Algavirales", "", 50.0, 50.0, selected_mode="novelty-aware"
    )
    assert result["estimated_completeness_strategy"] == "strategy2_r2_below_gate"
    assert "r2_below_gate" in result["order_completeness_v2_ood_flag"]
    assert result["completeness_model_reliability"] == "advisory_only"


def test_r2_gate_high_quality_above_floor():
    scorer = _tiered_scorer(test_r2=0.9921)  # >= R2_HIGH_QUALITY (0.7)
    result = scorer.calculate(
        {"OG1": 1, "OG2": 1}, "Algavirales", "", 50.0, 50.0, selected_mode="novelty-aware"
    )
    assert result["completeness_model_reliability"] == "high"
    assert result["estimated_completeness_strategy"] == "novelty_aware_v1"


# --- M13: non-legacy fallback must use strategy1_score, not 0.0 -------------


def test_untiered_order_uses_strategy1_not_zero():
    scorer = _make_scorer(models={}, metadata={}, tiers={}, baselines={})
    result = scorer.calculate(
        {}, "Imitervirales", "", 63.5, 63.5, selected_mode="novelty-aware"
    )
    assert result["estimated_completeness"] == 63.5
    assert result["estimated_completeness_strategy"] == "order_baseline_ratio_v1_no_tier_set"


def test_scorer_unavailable_uses_strategy1():
    scorer = _make_scorer(
        models={}, metadata={}, tiers={}, baselines={}, available=False
    )
    result = scorer.calculate(
        {}, "Imitervirales", "", 63.5, 63.5, selected_mode="novelty-aware"
    )
    assert result["estimated_completeness"] == 63.5
    assert result["estimated_completeness_strategy"].startswith("order_baseline_ratio_v1")


def test_empty_order_tax_uses_strategy1():
    scorer = _make_scorer(models={}, metadata={}, tiers={}, baselines={})
    result = scorer.calculate({}, "", "", 63.5, 63.5, selected_mode="novelty-aware")
    assert result["estimated_completeness"] == 63.5


def test_legacy_mode_fallback_unchanged():
    scorer = _make_scorer(models={}, metadata={}, tiers={}, baselines={})
    result = scorer.calculate(
        {}, "Imitervirales", "", 63.5, 63.5, selected_mode="legacy"
    )
    assert result["estimated_completeness"] == 63.5
    assert result["estimated_completeness_strategy"] == "order_baseline_ratio_v1"
