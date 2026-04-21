"""Novelty-aware completeness scoring for GVClass.

This module performs runtime inference only. All heavy training and reference
profile derivation are expected to be done offline and shipped as database
resources.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Tuple

import pandas as pd
from joblib import load

logger = logging.getLogger(__name__)

CONFIG_FILE = "completeness/config.json"
TIERS_FILE = "completeness/tiers.tsv"
BASELINES_FILE = "completeness/baselines.tsv"
MODEL_METADATA_FILE = "completeness/model_metadata.tsv"
MODEL_BUNDLE_FILE = "completeness/model.joblib"

OOD_STRICT_FLAGS = {
    "unassigned",
    "no_informative_markers",
    "zero_informative_hits",
    "r2_below_gate",
}

#: R² thresholds for the Phase 3.3 gate on the novelty-aware regressor.
#:
#: * ``r2_holdout < R2_ADVISORY_FLOOR`` — the regressor is worse than a
#:   trivial baseline; surface the strategy-2 tier score as the primary
#:   estimate, emit the ML score as ``estimated_completeness_advisory``
#:   only, and tag the record with OOD flag ``r2_below_gate``.
#: * ``R2_ADVISORY_FLOOR <= r2_holdout < R2_HIGH_QUALITY`` — surface the ML
#:   estimate but label quality as ``moderate``.
#: * ``r2_holdout >= R2_HIGH_QUALITY`` — quality ``high``.
#: * missing metadata — fall back to strategy-2 tiers with quality
#:   ``advisory_only`` and an OOD flag.
R2_ADVISORY_FLOOR = 0.5
R2_HIGH_QUALITY = 0.7


@dataclass(frozen=True)
class Strategy2Config:
    name: str
    core_prevalence: float
    shared_prevalence: float
    accessory_prevalence: float
    core_max_copy: float
    shared_max_copy: float
    family_min_ref_genomes: int
    family_min_informative_markers: int
    weight_core: float
    weight_shared: float
    weight_accessory: float


@dataclass(frozen=True)
class GroupDefinition:
    level: str
    name: str
    order_name: str


@dataclass(frozen=True)
class TierSet:
    core: Tuple[str, ...]
    shared: Tuple[str, ...]
    accessory: Tuple[str, ...]


class NoveltyAwareCompletenessScorer:
    """Runtime novelty-aware completeness scorer.

    This scorer expects resources produced by the offline builder script.
    """

    def __init__(self, database_path: Path):
        self.database_path = database_path
        self.config_path = database_path / CONFIG_FILE
        self.tiers_path = database_path / TIERS_FILE
        self.baselines_path = database_path / BASELINES_FILE
        self.model_metadata_path = database_path / MODEL_METADATA_FILE
        self.model_bundle_path = database_path / MODEL_BUNDLE_FILE
        self.available = False
        self.config: Strategy2Config | None = None
        self.feature_names: list[str] = []
        self.order_expected_marker_counts: Dict[str, int] = {}
        self.tiers: Dict[Tuple[str, str, str], TierSet] = {}
        self.baselines: Dict[Tuple[str, str, str], Dict[str, float]] = {}
        self.model_metadata: Dict[str, Dict[str, object]] = {}
        self.models: Dict[str, object] = {}
        self._load_resources()

    def _load_resources(self) -> None:
        required = [
            self.config_path,
            self.tiers_path,
            self.baselines_path,
            self.model_metadata_path,
            self.model_bundle_path,
        ]
        missing = [path.name for path in required if not path.exists()]
        if missing:
            logger.warning(
                "Novelty-aware completeness resources missing: %s",
                ", ".join(missing),
            )
            return

        try:
            config_data = json.loads(self.config_path.read_text())
            self.config = Strategy2Config(**config_data["strategy2_config"])
            self.order_expected_marker_counts = {
                str(key): int(value)
                for key, value in config_data["order_expected_marker_counts"].items()
            }

            tiers_df = pd.read_csv(self.tiers_path, sep="\t")
            tiers: Dict[Tuple[str, str, str], Dict[str, list[str]]] = {}
            for _, row in tiers_df.iterrows():
                key = (str(row["group_level"]), str(row["group_name"]), str(row["order_name"]))
                entry = tiers.setdefault(key, {"core": [], "shared": [], "accessory": []})
                tier_name = str(row["tier"])
                if tier_name in entry:
                    entry[tier_name].append(str(row["marker_name"]))
            self.tiers = {
                key: TierSet(
                    core=tuple(sorted(value["core"])),
                    shared=tuple(sorted(value["shared"])),
                    accessory=tuple(sorted(value["accessory"])),
                )
                for key, value in tiers.items()
            }

            baselines_df = pd.read_csv(self.baselines_path, sep="\t")
            self.baselines = {}
            for _, row in baselines_df.iterrows():
                key = (str(row["group_level"]), str(row["group_name"]), str(row["order_name"]))
                self.baselines[key] = {
                    "baseline_mean": float(row["baseline_mean"]),
                    "baseline_median": float(row["baseline_median"]),
                    "n_refs": int(row["n_refs"]),
                }

            model_meta_df = pd.read_csv(self.model_metadata_path, sep="\t")
            self.model_metadata = {
                str(row["order_name"]): row.to_dict() for _, row in model_meta_df.iterrows()
            }

            bundle = load(self.model_bundle_path)
            self.models = bundle["models"]
            self.feature_names = list(bundle["feature_names"])
            self.available = True
        except Exception as exc:
            logger.error("Failed to load novelty-aware completeness resources: %s", exc)
            self.available = False

    def select_group(self, order_tax: str, family_tax: str) -> Tuple[GroupDefinition | None, TierSet]:
        family_key = ("family", family_tax, order_tax)
        if family_tax and family_key in self.tiers:
            return GroupDefinition("family", family_tax, order_tax), self.tiers[family_key]
        order_key = ("order", order_tax, order_tax)
        if order_key in self.tiers:
            return GroupDefinition("order", order_tax, order_tax), self.tiers[order_key]
        return None, TierSet((), (), ())

    @staticmethod
    def _presence_fraction(marker_counts: Dict[str, int], markers: Tuple[str, ...]) -> float:
        if not markers:
            return 0.0
        present = sum(1 for marker in markers if marker_counts.get(marker, 0) > 0)
        return present / len(markers)

    def score_tier_set(self, marker_counts: Dict[str, int], tier_set: TierSet) -> Dict[str, float]:
        if self.config is None:
            return {
                "strategy2_raw": 0.0,
                "core_fraction": 0.0,
                "shared_fraction": 0.0,
                "accessory_fraction": 0.0,
                "core_marker_count": 0,
                "shared_marker_count": 0,
                "accessory_marker_count": 0,
                "informative_marker_count": 0,
            }

        core_fraction = self._presence_fraction(marker_counts, tier_set.core)
        shared_fraction = self._presence_fraction(marker_counts, tier_set.shared)
        accessory_fraction = self._presence_fraction(marker_counts, tier_set.accessory)

        used_weights = []
        weighted_values = []
        if tier_set.core:
            used_weights.append(self.config.weight_core)
            weighted_values.append(self.config.weight_core * core_fraction)
        if tier_set.shared:
            used_weights.append(self.config.weight_shared)
            weighted_values.append(self.config.weight_shared * shared_fraction)
        if tier_set.accessory:
            used_weights.append(self.config.weight_accessory)
            weighted_values.append(self.config.weight_accessory * accessory_fraction)

        if not used_weights:
            strategy2_raw = 0.0
        else:
            strategy2_raw = 100.0 * (sum(weighted_values) / sum(used_weights))

        return {
            "strategy2_raw": round(strategy2_raw, 2),
            "core_fraction": round(core_fraction * 100.0, 2),
            "shared_fraction": round(shared_fraction * 100.0, 2),
            "accessory_fraction": round(accessory_fraction * 100.0, 2),
            "core_marker_count": len(tier_set.core),
            "shared_marker_count": len(tier_set.shared),
            "accessory_marker_count": len(tier_set.accessory),
            "informative_marker_count": len(tier_set.core) + len(tier_set.shared),
        }

    def normalize_strategy2_score(
        self,
        raw_score: float,
        group_def: GroupDefinition,
    ) -> Tuple[float, float, int]:
        key = (group_def.level, group_def.name, group_def.order_name)
        baseline = self.baselines.get(key, {})
        baseline_mean = float(baseline.get("baseline_mean", 0.0))
        n_refs = int(baseline.get("n_refs", 0))
        if baseline_mean <= 0:
            return raw_score, baseline_mean, n_refs
        normalized = min(100.0, (raw_score / baseline_mean) * 100.0)
        return round(normalized, 2), baseline_mean, n_refs

    @staticmethod
    def compute_support_metrics(
        informative_marker_count: int,
        core_marker_count: int,
        shared_marker_count: int,
        core_fraction: float,
        shared_fraction: float,
        baseline_n_refs: int,
        group_level: str,
    ) -> Dict[str, float | int | str]:
        core_present = round((core_fraction / 100.0) * core_marker_count)
        shared_present = round((shared_fraction / 100.0) * shared_marker_count)
        informative_present = int(core_present + shared_present)
        informative_fraction = (
            0.0 if informative_marker_count == 0 else informative_present / informative_marker_count
        )
        ref_support = min(1.0, baseline_n_refs / 10.0)
        support_score = round((0.85 * informative_fraction) + (0.15 * ref_support), 4)

        if group_level == "unavailable":
            ood_flag = "unassigned"
        elif informative_marker_count == 0:
            ood_flag = "no_informative_markers"
        elif informative_present == 0:
            ood_flag = "zero_informative_hits"
        elif support_score < 0.2:
            ood_flag = "low_support"
        elif baseline_n_refs < 5:
            ood_flag = "low_reference_support"
        else:
            ood_flag = "ok"

        return {
            "informative_marker_count": informative_marker_count,
            "informative_present_count": informative_present,
            "informative_fraction": round(informative_fraction * 100.0, 2),
            "support_score": round(support_score * 100.0, 2),
            "ood_flag": ood_flag,
        }

    @staticmethod
    def apply_ood_cap(
        prediction: float,
        strategy1_score: float,
        strategy2_score: float,
        ood_flag: str,
        support_score: float,
        baseline_n_refs: int,
    ) -> float:
        conservative_floor = max(strategy1_score, strategy2_score)
        support = max(0.0, min(1.0, support_score / 100.0))
        ref_support = max(0.0, min(1.0, baseline_n_refs / 20.0))
        if ood_flag in OOD_STRICT_FLAGS:
            if ood_flag == "unassigned":
                return conservative_floor
            shrink = 0.20 + (0.40 * ref_support) + (0.20 * support)
            adjusted = conservative_floor + ((prediction - conservative_floor) * shrink)
            return max(conservative_floor, min(prediction, adjusted))
        if ood_flag == "low_support":
            shrink = 0.45 + (0.30 * ref_support) + (0.20 * support)
            adjusted = conservative_floor + ((prediction - conservative_floor) * shrink)
            return max(conservative_floor, min(prediction, adjusted))
        if ood_flag == "low_reference_support":
            shrink = 0.60 + (0.15 * ref_support) + (0.10 * support)
            adjusted = conservative_floor + ((prediction - conservative_floor) * shrink)
            return max(conservative_floor, min(prediction, adjusted))
        return prediction

    def default_metrics(self, order_tax: str = "", family_tax: str = "") -> Dict[str, object]:
        return {
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
            "estimated_completeness_strategy": "novelty_aware_v1",
        }

    def calculate(
        self,
        marker_counts: Dict[str, int],
        order_tax: str,
        family_tax: str,
        strategy1_raw: float,
        strategy1_score: float,
        selected_mode: str = "legacy",
    ) -> Dict[str, object]:
        if not self.available or not order_tax:
            metrics = self.default_metrics(order_tax, family_tax)
            metrics["estimated_completeness"] = (
                strategy1_score if selected_mode == "legacy" else metrics["order_completeness_v2"]
            )
            metrics["estimated_completeness_strategy"] = (
                "order_baseline_ratio_v1" if selected_mode == "legacy" else metrics["order_completeness_v2_strategy"]
            )
            return metrics

        group_def, tier_set = self.select_group(order_tax, family_tax)
        if group_def is None:
            metrics = self.default_metrics(order_tax, family_tax)
            metrics["estimated_completeness"] = (
                strategy1_score if selected_mode == "legacy" else metrics["order_completeness_v2"]
            )
            metrics["estimated_completeness_strategy"] = (
                "order_baseline_ratio_v1" if selected_mode == "legacy" else metrics["order_completeness_v2_strategy"]
            )
            return metrics

        tier_scores = self.score_tier_set(marker_counts, tier_set)
        strategy2_norm, baseline_mean, baseline_n_refs = self.normalize_strategy2_score(
            float(tier_scores["strategy2_raw"]),
            group_def,
        )
        support_metrics = self.compute_support_metrics(
            informative_marker_count=int(tier_scores["informative_marker_count"]),
            core_marker_count=int(tier_scores["core_marker_count"]),
            shared_marker_count=int(tier_scores["shared_marker_count"]),
            core_fraction=float(tier_scores["core_fraction"]),
            shared_fraction=float(tier_scores["shared_fraction"]),
            baseline_n_refs=int(baseline_n_refs),
            group_level=group_def.level,
        )

        validation_mode = "fallback"
        ml_prediction_raw: float | None = None
        prediction = float(strategy2_norm)
        model = self.models.get(order_tax)
        meta = self.model_metadata.get(order_tax, {})
        if model is not None:
            validation_mode = str(meta.get("validation_type", "trained"))
            features = {
                "raw_score": float(strategy1_raw),
                "strategy1_score": float(strategy1_score),
                "strategy2_raw_score": float(tier_scores["strategy2_raw"]),
                "strategy2_score": float(strategy2_norm),
                "strategy2_baseline_mean": float(baseline_mean),
                "strategy2_baseline_n_refs": float(baseline_n_refs),
                "core_fraction": float(tier_scores["core_fraction"]),
                "shared_fraction": float(tier_scores["shared_fraction"]),
                "accessory_fraction": float(tier_scores["accessory_fraction"]),
                "n_expected_markers": float(self.order_expected_marker_counts.get(order_tax, 0)),
                "n_core_markers": float(tier_scores["core_marker_count"]),
                "n_shared_markers": float(tier_scores["shared_marker_count"]),
                "n_accessory_markers": float(tier_scores["accessory_marker_count"]),
                "group_is_family": 1.0 if group_def.level == "family" else 0.0,
            }
            X = pd.DataFrame([{name: features.get(name, 0.0) for name in self.feature_names}])
            ml_prediction_raw = float(model.predict(X)[0])
            prediction = ml_prediction_raw

        bounded = max(0.0, min(100.0, prediction))
        final_score = self.apply_ood_cap(
            prediction=bounded,
            strategy1_score=float(strategy1_score),
            strategy2_score=float(strategy2_norm),
            ood_flag=str(support_metrics["ood_flag"]),
            support_score=float(support_metrics["support_score"]),
            baseline_n_refs=int(baseline_n_refs),
        )
        novelty_score = round(final_score, 2)

        # Phase 3.3: gate the ML estimate on ``r2_holdout``. Below the
        # advisory floor we do NOT surface the ML prediction as the primary
        # estimate — the regressor is statistically no better than
        # strategy-2 tiers for this order. Surface the ML score as an
        # advisory column instead so users who trust it can still read it.
        r2_raw = meta.get("r2_holdout")
        try:
            if r2_raw in (None, ""):
                r2_holdout = None
            else:
                r2_holdout = float(r2_raw)
                # NaN (from "nan" / "NaN" literals in the TSV or a failed
                # training run) must be treated as missing, not as a
                # finite number that falls through to the high-quality
                # branch (Codex audit).
                import math as _math

                if _math.isnan(r2_holdout):
                    r2_holdout = None
        except (TypeError, ValueError):
            r2_holdout = None

        if model is None:
            # No ML regressor for this order at all.
            primary_score = float(strategy2_norm)
            primary_strategy = "strategy2_no_ml_profile"
            quality = "advisory_only"
            advisory_score = None
            ood_augmented = str(support_metrics["ood_flag"])
        elif r2_holdout is None:
            # Model present but no metadata — trust the ML score but flag
            # it as advisory-only since we cannot verify hold-out quality.
            primary_score = novelty_score
            primary_strategy = "novelty_aware_v1"
            quality = "advisory_only"
            advisory_score = (
                round(float(ml_prediction_raw), 2)
                if ml_prediction_raw is not None
                else None
            )
            ood_augmented = str(support_metrics["ood_flag"])
        elif r2_holdout < R2_ADVISORY_FLOOR:
            # Below the hard gate — ML is worse than strategy-2 in holdout.
            primary_score = float(strategy2_norm)
            primary_strategy = "strategy2_r2_below_gate"
            quality = "advisory_only"
            advisory_score = (
                round(float(ml_prediction_raw), 2)
                if ml_prediction_raw is not None
                else None
            )
            ood_augmented = self._augment_ood_flag(
                support_metrics["ood_flag"], "r2_below_gate"
            )
        elif r2_holdout < R2_HIGH_QUALITY:
            primary_score = novelty_score
            primary_strategy = "novelty_aware_v1"
            quality = "moderate"
            advisory_score = (
                round(float(ml_prediction_raw), 2)
                if ml_prediction_raw is not None
                else None
            )
            ood_augmented = str(support_metrics["ood_flag"])
        else:
            primary_score = novelty_score
            primary_strategy = "novelty_aware_v1"
            quality = "high"
            advisory_score = (
                round(float(ml_prediction_raw), 2)
                if ml_prediction_raw is not None
                else None
            )
            ood_augmented = str(support_metrics["ood_flag"])

        estimated_score = (
            strategy1_score if selected_mode == "legacy" else primary_score
        )
        estimated_strategy = (
            "order_baseline_ratio_v1" if selected_mode == "legacy" else primary_strategy
        )

        result: Dict[str, Any] = {
            "order_completeness_v2": novelty_score,
            "order_completeness_v2_strategy": "novelty_aware_v1",
            "order_completeness_v2_strategy2_raw": round(float(tier_scores["strategy2_raw"]), 2),
            "order_completeness_v2_strategy2_normalized": round(float(strategy2_norm), 2),
            "order_completeness_v2_support_score": round(float(support_metrics["support_score"]), 2),
            "order_completeness_v2_ood_flag": ood_augmented,
            "order_completeness_v2_reference_group": f"{group_def.level}:{group_def.name}",
            "order_completeness_v2_validation_mode": validation_mode,
            "order_completeness_v2_informative_fraction": round(float(support_metrics["informative_fraction"]), 2),
            "estimated_completeness": round(float(estimated_score), 2),
            "estimated_completeness_strategy": estimated_strategy,
            "estimated_completeness_quality": quality,
            "estimated_completeness_r2_holdout": (
                round(r2_holdout, 4) if r2_holdout is not None else ""
            ),
            "estimated_completeness_advisory": (
                advisory_score if advisory_score is not None else ""
            ),
        }
        return result

    @staticmethod
    def _augment_ood_flag(current: str, additional: str) -> str:
        """Merge ``additional`` into the existing ood_flag string without
        losing the previous flag. Flags are comma-separated for downstream
        consumers that parse them."""
        current = str(current or "").strip()
        if not current or current == "none":
            return additional
        if additional in current.split(","):
            return current
        return f"{current},{additional}"


def create_novelty_completeness_scorer(database_path: Path) -> NoveltyAwareCompletenessScorer:
    return NoveltyAwareCompletenessScorer(database_path)
