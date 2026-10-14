"""Build novelty-aware completeness resources for GVClass.

This is an offline resource-generation step. It derives strategy-2 tier sets and
baselines from the reference database and trains the novelty-aware strategy-3
models with leave-one-family-out validation where possible.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import shutil
import tempfile
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Any, Dict, List, Tuple

import pandas as pd
from joblib import dump
from sklearn.ensemble import (
    ExtraTreesRegressor,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.metrics import mean_absolute_error, r2_score

from src.core.novelty_completeness import GroupDefinition, Strategy2Config, TierSet
from src.utils.resource_store import ResourceStore

RANDOM_SEED = 42
SIMULATED_COMPLETENESS_LEVELS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
SIMULATIONS_PER_LEVEL_FINAL = 3
MAX_SYNTHETIC_SAMPLES_PER_FAMILY = 2500
MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY = 8000

CONFIG = Strategy2Config(
    name="cfg_241",
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

FEATURE_NAMES = [
    "raw_score",
    "strategy1_score",
    "strategy2_raw_score",
    "strategy2_score",
    "strategy2_baseline_mean",
    "strategy2_baseline_n_refs",
    "core_fraction",
    "shared_fraction",
    "accessory_fraction",
    "n_expected_markers",
    "n_core_markers",
    "n_shared_markers",
    "n_accessory_markers",
    "group_is_family",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resources-dir", default="resources")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="New directory for the generated resources; the source directory is never overwritten",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Use small model candidates for a fast builder check",
    )
    return parser.parse_args()


def load_order_table(
    resources_dir: Path,
) -> Tuple[Dict[str, List[str]], List[Dict[str, str]], List[str]]:
    order_markers: Dict[str, List[str]] = {}
    # v1.5.0 moved order_completeness.tab under markers/; fall back to the
    # pre-reorg top-level location for older resource layouts.
    path = resources_dir / "markers" / "order_completeness.tab"
    if not path.exists():
        path = resources_dir / "order_completeness.tab"
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
        for row in rows:
            order_name = row["Order"]
            order_markers[order_name] = [
                marker.strip()
                for marker in row["Orthogroups"].split(",")
                if marker.strip()
            ]
    return order_markers, rows, fieldnames


def load_labels(resources_dir: Path) -> Dict[str, Dict[str, str]]:
    labels: Dict[str, Dict[str, str]] = {}
    path = ResourceStore(resources_dir).label_path("labels.tsv")
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            genome_id = parts[0]
            tax = parts[1].split("|")
            while len(tax) < 7:
                tax.append("")
            labels[genome_id] = {
                "domain": tax[0],
                "phylum": tax[1],
                "class": tax[2],
                "order": tax[3],
                "family": tax[4],
                "genus": tax[5],
                "species": tax[6],
            }
    return labels


def parse_reference_marker_profiles(
    resources_dir: Path,
    labels: Dict[str, Dict[str, str]],
    order_markers: Dict[str, List[str]],
) -> Dict[str, Dict[str, Dict[str, int]]]:
    relevant_markers = sorted({marker for markers in order_markers.values() for marker in markers})
    order_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    store = ResourceStore(resources_dir)

    for marker in relevant_markers:
        try:
            marker_path = store.marker_faa_path(marker)
        except FileNotFoundError:
            continue
        genome_copy_counts: Dict[str, int] = Counter()
        with marker_path.open() as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                genome_id = line[1:].split("|", 1)[0].strip()
                genome_copy_counts[genome_id] += 1

        for genome_id, copy_count in genome_copy_counts.items():
            tax = labels.get(genome_id)
            if not tax:
                continue
            if "MIRUS" in order_markers and (
                tax["domain"] == "MIRUS" or tax["phylum"] == "MIRUS"
            ):
                order_name = "MIRUS"
            else:
                order_name = tax["order"]
            if order_name not in order_markers:
                continue
            order_profiles[order_name][genome_id][marker] = copy_count
    return order_profiles


def marker_stats_for_group(
    genome_profiles: Dict[str, Dict[str, int]], markers: List[str]
) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    genome_count = len(genome_profiles)
    if genome_count == 0:
        return pd.DataFrame()
    for marker in markers:
        copy_counts = [profile.get(marker, 0) for profile in genome_profiles.values()]
        present = sum(1 for count in copy_counts if count > 0)
        prevalence = present / genome_count
        avg_copy = sum(copy_counts) / genome_count
        rows.append(
            {
                "marker": marker,
                "prevalence": prevalence,
                "avg_copy": avg_copy,
                "n_present": present,
                "n_genomes": genome_count,
            }
        )
    return pd.DataFrame(rows)


def build_tier_set(group_stats: pd.DataFrame) -> TierSet:
    if group_stats.empty:
        return TierSet((), (), ())
    core = sorted(
        group_stats[
            (group_stats["prevalence"] >= CONFIG.core_prevalence)
            & (group_stats["avg_copy"] <= CONFIG.core_max_copy)
        ]["marker"].tolist()
    )
    shared = sorted(
        group_stats[
            (group_stats["prevalence"] >= CONFIG.shared_prevalence)
            & (group_stats["avg_copy"] <= CONFIG.shared_max_copy)
            & (group_stats["prevalence"] < CONFIG.core_prevalence)
        ]["marker"].tolist()
    )
    accessory = sorted(
        group_stats[
            (group_stats["prevalence"] >= CONFIG.accessory_prevalence)
            & (group_stats["prevalence"] < CONFIG.shared_prevalence)
        ]["marker"].tolist()
    )
    return TierSet(tuple(core), tuple(shared), tuple(accessory))


def informative_marker_count(tier_set: TierSet) -> int:
    return len(tier_set.core) + len(tier_set.shared)


def presence_fraction(marker_counts: Dict[str, int], markers: Tuple[str, ...]) -> float:
    if not markers:
        return 0.0
    present = sum(1 for marker in markers if marker_counts.get(marker, 0) > 0)
    return present / len(markers)


def score_tier_set(
    marker_counts: Dict[str, int], tier_set: TierSet
) -> Dict[str, float]:
    core_fraction = presence_fraction(marker_counts, tier_set.core)
    shared_fraction = presence_fraction(marker_counts, tier_set.shared)
    accessory_fraction = presence_fraction(marker_counts, tier_set.accessory)

    used_weights = []
    weighted_values = []
    if tier_set.core:
        used_weights.append(CONFIG.weight_core)
        weighted_values.append(CONFIG.weight_core * core_fraction)
    if tier_set.shared:
        used_weights.append(CONFIG.weight_shared)
        weighted_values.append(CONFIG.weight_shared * shared_fraction)
    if tier_set.accessory:
        used_weights.append(CONFIG.weight_accessory)
        weighted_values.append(CONFIG.weight_accessory * accessory_fraction)

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
    }


def build_group_resources(
    order_markers: Dict[str, List[str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_profiles: Dict[Tuple[str, str], Dict[str, Dict[str, int]]],
) -> Tuple[
    Dict[Tuple[str, str, str], TierSet],
    Dict[Tuple[str, str, str], Dict[str, float]],
]:
    tiers: Dict[Tuple[str, str, str], TierSet] = {}
    baselines: Dict[Tuple[str, str, str], Dict[str, float]] = {}

    for order_name, genome_profiles in order_profiles.items():
        order_stats = marker_stats_for_group(genome_profiles, order_markers[order_name])
        tier_set = build_tier_set(order_stats)
        tiers[("order", order_name, order_name)] = tier_set
        scores = [
            float(score_tier_set(profile, tier_set)["strategy2_raw"])
            for profile in genome_profiles.values()
        ]
        baselines[("order", order_name, order_name)] = {
            "baseline_mean": round(mean(scores), 4),
            "baseline_median": round(median(scores), 4),
            "n_refs": len(scores),
        }

    for (order_name, family_name), genome_profiles in family_profiles.items():
        if len(genome_profiles) < CONFIG.family_min_ref_genomes:
            continue
        family_stats = marker_stats_for_group(
            genome_profiles, order_markers[order_name]
        )
        tier_set = build_tier_set(family_stats)
        if informative_marker_count(tier_set) < CONFIG.family_min_informative_markers:
            continue
        tiers[("family", family_name, order_name)] = tier_set
        scores = [
            float(score_tier_set(profile, tier_set)["strategy2_raw"])
            for profile in genome_profiles.values()
        ]
        baselines[("family", family_name, order_name)] = {
            "baseline_mean": round(mean(scores), 4),
            "baseline_median": round(median(scores), 4),
            "n_refs": len(scores),
        }

    return tiers, baselines


def select_group(
    tiers: Dict[Tuple[str, str, str], TierSet],
    order_name: str,
    family_name: str,
) -> Tuple[GroupDefinition, TierSet]:
    family_key = ("family", family_name, order_name)
    if family_name and family_key in tiers:
        return GroupDefinition("family", family_name, order_name), tiers[family_key]
    order_key = ("order", order_name, order_name)
    return GroupDefinition("order", order_name, order_name), tiers[order_key]


def normalize_strategy2_score(
    raw_score: float,
    group_def: GroupDefinition,
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
) -> Tuple[float, float, int]:
    key = (group_def.level, group_def.name, group_def.order_name)
    baseline = baselines[key]
    baseline_mean = float(baseline["baseline_mean"])
    n_refs = int(baseline["n_refs"])
    if baseline_mean <= 0:
        return raw_score, baseline_mean, n_refs
    normalized = min(100.0, (raw_score / baseline_mean) * 100.0)
    return round(normalized, 2), baseline_mean, n_refs


def simulate_samples(
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    seed: int = RANDOM_SEED,
) -> List[Dict[str, object]]:
    rng = random.Random(seed)
    samples: List[Dict[str, object]] = []
    for order_name, genome_profiles in order_profiles.items():
        markers = order_markers[order_name]
        for genome_id, profile in genome_profiles.items():
            present_markers = [
                marker for marker in markers if profile.get(marker, 0) > 0
            ]
            if not present_markers:
                continue
            family_name = labels.get(genome_id, {}).get("family", "")
            for completeness_level in SIMULATED_COMPLETENESS_LEVELS:
                retain_n = max(1, round(len(present_markers) * completeness_level))
                repetitions = (
                    1
                    if retain_n == len(present_markers)
                    else SIMULATIONS_PER_LEVEL_FINAL
                )
                for rep in range(repetitions):
                    retained = set(rng.sample(present_markers, retain_n))
                    samples.append(
                        {
                            "source_genome_id": genome_id,
                            "order_name": order_name,
                            "family_name": family_name,
                            "true_completeness": 100.0
                            * retain_n
                            / len(present_markers),
                            "present_markers": retained,
                            "rep": rep,
                        }
                    )
    return samples


def downsample_training_order_df(order_df: pd.DataFrame) -> pd.DataFrame:
    family_names = [
        family for family in sorted(order_df["family_name"].unique()) if family
    ]
    if len(family_names) >= 2:
        parts = []
        for _, family_df in order_df.groupby("family_name", sort=False):
            if len(family_df) > MAX_SYNTHETIC_SAMPLES_PER_FAMILY:
                parts.append(
                    family_df.sample(
                        n=MAX_SYNTHETIC_SAMPLES_PER_FAMILY, random_state=RANDOM_SEED
                    )
                )
            else:
                parts.append(family_df)
        return pd.concat(parts, ignore_index=True)
    if len(order_df) > MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY:
        return order_df.sample(
            n=MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY, random_state=RANDOM_SEED
        ).reset_index(drop=True)
    return order_df.reset_index(drop=True)


def make_validation_plan(
    genome_profiles: Dict[str, Dict[str, int]],
    labels: Dict[str, Dict[str, str]],
) -> Dict[str, Any] | None:
    """Split source genomes before any reference-derived feature is calculated."""
    genome_ids = sorted(genome_profiles)
    families: Dict[str, List[str]] = defaultdict(list)
    for genome_id in genome_ids:
        family = labels.get(genome_id, {}).get("family", "")
        if family:
            families[family].append(genome_id)

    if len(families) >= 3:
        target_size = max(1, round(len(genome_ids) * 0.2))
        final_family = min(
            families,
            key=lambda family: (
                abs(len(families[family]) - target_size),
                len(families[family]),
                family,
            ),
        )
        final_ids = sorted(families[final_family])
        development_ids = [gid for gid in genome_ids if gid not in set(final_ids)]
        folds = []
        for family in sorted(
            {
                labels.get(genome_id, {}).get("family", "")
                for genome_id in development_ids
                if labels.get(genome_id, {}).get("family", "")
            }
        ):
            validation_ids = [
                genome_id
                for genome_id in development_ids
                if labels.get(genome_id, {}).get("family", "") == family
            ]
            training_ids = [
                genome_id for genome_id in development_ids if genome_id not in set(validation_ids)
            ]
            if len(training_ids) >= 2 and validation_ids:
                folds.append(
                    {
                        "held_out_group": family,
                        "training_ids": training_ids,
                        "validation_ids": validation_ids,
                    }
                )
        if len(folds) >= 2:
            return {
                "split_unit": "family",
                "final_held_out_group": final_family,
                "development_ids": development_ids,
                "final_test_ids": final_ids,
                "selection_folds": folds,
            }

    if len(genome_ids) < 4:
        return None
    final_n = min(max(1, round(len(genome_ids) * 0.2)), len(genome_ids) - 3)
    development_ids = genome_ids[:-final_n]
    final_ids = genome_ids[-final_n:]
    n_folds = min(5, len(development_ids))
    folds = []
    for fold_index in range(n_folds):
        validation_ids = development_ids[fold_index::n_folds]
        training_ids = [
            genome_id for genome_id in development_ids if genome_id not in set(validation_ids)
        ]
        if len(training_ids) >= 2 and validation_ids:
            folds.append(
                {
                    "held_out_group": ",".join(validation_ids),
                    "training_ids": training_ids,
                    "validation_ids": validation_ids,
                }
            )
    if len(folds) < 2:
        return None
    return {
        "split_unit": "genome",
        "final_held_out_group": ",".join(final_ids),
        "development_ids": development_ids,
        "final_test_ids": final_ids,
        "selection_folds": folds,
    }


def reference_resources_for_ids(
    order_name: str,
    genome_profiles: Dict[str, Dict[str, int]],
    genome_ids: List[str],
    labels: Dict[str, Dict[str, str]],
    order_markers: Dict[str, List[str]],
) -> Tuple[
    Dict[str, Dict[str, Dict[str, int]]],
    Dict[Tuple[str, str], Dict[str, Dict[str, int]]],
    Dict[Tuple[str, str, str], TierSet],
    Dict[Tuple[str, str, str], Dict[str, float]],
    float,
]:
    selected = {
        genome_id: genome_profiles[genome_id]
        for genome_id in genome_ids
        if genome_id in genome_profiles
    }
    selected_order_profiles = {order_name: selected}
    selected_family_profiles: Dict[
        Tuple[str, str], Dict[str, Dict[str, int]]
    ] = defaultdict(dict)
    for genome_id, profile in selected.items():
        family = labels.get(genome_id, {}).get("family", "")
        selected_family_profiles[(order_name, family)][genome_id] = profile
    tiers, baselines = build_group_resources(
        order_markers, selected_order_profiles, selected_family_profiles
    )
    strategy1_baseline, _ = order_recovery_stats(
        selected, order_markers[order_name]
    )
    return (
        selected_order_profiles,
        selected_family_profiles,
        tiers,
        baselines,
        strategy1_baseline,
    )


def order_recovery_stats(
    genome_profiles: Dict[str, Dict[str, int]], markers: List[str]
) -> Tuple[float, float]:
    marker_recovery = [
        100.0
        * sum(1 for marker in markers if profile.get(marker, 0) > 0)
        / len(markers)
        for profile in genome_profiles.values()
    ]
    if not marker_recovery:
        return 0.0, 0.0
    return mean(marker_recovery), pstdev(marker_recovery)


def build_training_df(
    labels: Dict[str, Dict[str, str]],
    sample_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    tiers: Dict[Tuple[str, str, str], TierSet],
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
    strategy1_baselines: Dict[str, float],
    seed: int = RANDOM_SEED,
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for sample in simulate_samples(labels, sample_profiles, order_markers, seed=seed):
        order_name = str(sample["order_name"])
        family_name = str(sample["family_name"])
        group_def, tier_set = select_group(tiers, order_name, family_name)
        raw_score = 100.0 * (len(sample["present_markers"]) / len(order_markers[order_name]))
        baseline = strategy1_baselines[order_name]
        strategy1_score = 0.0 if baseline <= 0 else min(100.0, (raw_score / baseline) * 100.0)
        marker_counts = {marker: 1 for marker in sample["present_markers"]}
        tier_scores = score_tier_set(marker_counts, tier_set)
        strategy2_raw_score = float(tier_scores["strategy2_raw"])
        strategy2_norm_score, strategy2_baseline_mean, strategy2_baseline_n_refs = (
            normalize_strategy2_score(
                strategy2_raw_score,
                group_def,
                baselines,
            )
        )
        rows.append(
            {
                "source_genome_id": sample["source_genome_id"],
                "order_name": order_name,
                "family_name": family_name,
                "raw_score": raw_score,
                "strategy1_score": strategy1_score,
                "strategy2_raw_score": strategy2_raw_score,
                "strategy2_score": strategy2_norm_score,
                "strategy2_baseline_mean": strategy2_baseline_mean,
                "strategy2_baseline_n_refs": strategy2_baseline_n_refs,
                "core_fraction": tier_scores["core_fraction"],
                "shared_fraction": tier_scores["shared_fraction"],
                "accessory_fraction": tier_scores["accessory_fraction"],
                "n_expected_markers": len(order_markers[order_name]),
                "n_core_markers": tier_scores["core_marker_count"],
                "n_shared_markers": tier_scores["shared_marker_count"],
                "n_accessory_markers": tier_scores["accessory_marker_count"],
                "group_is_family": 1 if group_def.level == "family" else 0,
                "true_completeness": float(sample["true_completeness"]),
                "target_definition": "marker_dropout_fraction_of_source_profile",
                "target_is_nucleotide_completeness": False,
            }
        )
    return pd.DataFrame(rows)


def _candidate_models(smoke: bool):
    trees = 20 if smoke else 400
    return [
        (
            "random_forest",
            lambda: RandomForestRegressor(
                n_estimators=trees,
                random_state=RANDOM_SEED,
                min_samples_leaf=2,
            ),
        ),
        (
            "extra_trees",
            lambda: ExtraTreesRegressor(
                n_estimators=trees if smoke else 500,
                random_state=RANDOM_SEED,
                min_samples_leaf=2,
            ),
        ),
        (
            "hist_gbm",
            lambda: HistGradientBoostingRegressor(
                random_state=RANDOM_SEED,
                max_depth=4 if smoke else 8,
                min_samples_leaf=2 if smoke else 10,
            ),
        ),
    ]


def _metrics(truth: List[float], predictions: List[float]) -> Tuple[float, float]:
    mae = float(mean_absolute_error(truth, predictions))
    r2 = float(r2_score(truth, predictions)) if len(set(truth)) > 1 else float("nan")
    return mae, r2


def train_models(
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    smoke: bool = False,
) -> Tuple[
    Dict[str, object],
    pd.DataFrame,
    Dict[Tuple[str, str, str], TierSet],
    Dict[Tuple[str, str, str], Dict[str, float]],
    Dict[str, Dict[str, Dict[str, int]]],
    List[Dict[str, str]],
]:
    models: Dict[str, object] = {}
    metadata_rows: List[Dict[str, object]] = []
    output_tiers: Dict[Tuple[str, str, str], TierSet] = {}
    output_baselines: Dict[Tuple[str, str, str], Dict[str, float]] = {}
    output_profiles: Dict[str, Dict[str, Dict[str, int]]] = {}
    split_rows: List[Dict[str, str]] = []

    for order_name, genome_profiles in sorted(order_profiles.items()):
        plan = make_validation_plan(genome_profiles, labels)
        common_metadata = {
            "order_name": order_name,
            "training_target": "marker_dropout_fraction_of_source_profile",
            "target_is_nucleotide_completeness": False,
            "reference_ascertainment": "labeled_references_with_detected_order_markers",
            "marker_panel_ascertainment": "fixed_input_panel_predates_validation_holdouts",
            "model_reliability": "advisory_only",
            "test_r2": "",
            "test_mae": "",
        }
        if plan is None:
            fallback_profiles, _, fallback_tiers, fallback_baselines, _ = (
                reference_resources_for_ids(
                    order_name,
                    genome_profiles,
                    list(genome_profiles),
                    labels,
                    order_markers,
                )
            )
            output_tiers.update(fallback_tiers)
            output_baselines.update(fallback_baselines)
            output_profiles.update(fallback_profiles)
            metadata_rows.append(
                {
                    **common_metadata,
                    "training_status": "unsupported",
                    "unsupported_reason": "fewer_than_four_source_genomes_or_two_selection_folds",
                    "model_name": "",
                    "validation_type": "unsupported_grouped_validation",
                    "n_training_samples": 0,
                    "n_training_samples_used": 0,
                }
            )
            continue

        fold_frames = []
        for fold_index, fold in enumerate(plan["selection_folds"]):
            (
                training_profiles,
                _,
                fold_tiers,
                fold_baselines,
                strategy1_baseline,
            ) = reference_resources_for_ids(
                order_name,
                genome_profiles,
                fold["training_ids"],
                labels,
                order_markers,
            )
            validation_profiles = {
                order_name: {
                    genome_id: genome_profiles[genome_id]
                    for genome_id in fold["validation_ids"]
                }
            }
            training_df = build_training_df(
                labels,
                training_profiles,
                order_markers,
                fold_tiers,
                fold_baselines,
                {order_name: strategy1_baseline},
                seed=RANDOM_SEED + fold_index,
            )
            training_df = downsample_training_order_df(training_df)
            validation_df = build_training_df(
                labels,
                validation_profiles,
                order_markers,
                fold_tiers,
                fold_baselines,
                {order_name: strategy1_baseline},
                seed=RANDOM_SEED + 100 + fold_index,
            )
            if len(training_df) >= 20 and not validation_df.empty:
                fold_frames.append((training_df, validation_df))
                for role, genome_ids in (
                    ("selection_train", fold["training_ids"]),
                    ("selection_validation", fold["validation_ids"]),
                ):
                    for genome_id in genome_ids:
                        split_rows.append(
                            {
                                "order_name": order_name,
                                "stage": "candidate_selection",
                                "fold": str(fold_index),
                                "role": role,
                                "split_unit": str(plan["split_unit"]),
                                "genome_id": genome_id,
                                "family_name": labels.get(genome_id, {}).get(
                                    "family", ""
                                ),
                            }
                        )

        if len(fold_frames) < 2:
            fallback_profiles, _, fallback_tiers, fallback_baselines, _ = (
                reference_resources_for_ids(
                    order_name,
                    genome_profiles,
                    plan["development_ids"],
                    labels,
                    order_markers,
                )
            )
            output_tiers.update(fallback_tiers)
            output_baselines.update(fallback_baselines)
            output_profiles.update(fallback_profiles)
            metadata_rows.append(
                {
                    **common_metadata,
                    "training_status": "unsupported",
                    "unsupported_reason": "fewer_than_two_usable_grouped_selection_folds",
                    "model_name": "",
                    "validation_type": "unsupported_grouped_validation",
                    "n_training_samples": 0,
                    "n_training_samples_used": 0,
                }
            )
            continue

        best_name = ""
        best_factory = None
        best_mae = None
        best_r2 = float("nan")
        for model_name, model_factory in _candidate_models(smoke):
            predictions: List[float] = []
            truth: List[float] = []
            for training_df, validation_df in fold_frames:
                model = model_factory()
                model.fit(training_df[FEATURE_NAMES], training_df["true_completeness"])
                predictions.extend(model.predict(validation_df[FEATURE_NAMES]).tolist())
                truth.extend(validation_df["true_completeness"].tolist())
            mae, r2 = _metrics(truth, predictions)
            if best_mae is None or mae < best_mae:
                best_name = model_name
                best_factory = model_factory
                best_mae = mae
                best_r2 = r2

        development_ids = plan["development_ids"]
        (
            development_profiles,
            _,
            development_tiers,
            development_baselines,
            strategy1_baseline,
        ) = reference_resources_for_ids(
            order_name,
            genome_profiles,
            development_ids,
            labels,
            order_markers,
        )
        development_df = downsample_training_order_df(
            build_training_df(
                labels,
                development_profiles,
                order_markers,
                development_tiers,
                development_baselines,
                {order_name: strategy1_baseline},
            )
        )
        final_profiles = {
            order_name: {
                genome_id: genome_profiles[genome_id]
                for genome_id in plan["final_test_ids"]
            }
        }
        final_df = build_training_df(
            labels,
            final_profiles,
            order_markers,
            development_tiers,
            development_baselines,
            {order_name: strategy1_baseline},
            seed=RANDOM_SEED + 1000,
        )
        if best_factory is None or development_df.empty or final_df.empty:
            output_tiers.update(development_tiers)
            output_baselines.update(development_baselines)
            output_profiles.update(development_profiles)
            metadata_rows.append(
                {
                    **common_metadata,
                    "training_status": "unsupported",
                    "unsupported_reason": "empty_development_or_final_test_features",
                    "model_name": "",
                    "validation_type": "unsupported_grouped_validation",
                    "n_training_samples": len(development_df),
                    "n_training_samples_used": len(development_df),
                }
            )
            continue

        best_model = best_factory()
        best_model.fit(development_df[FEATURE_NAMES], development_df["true_completeness"])
        final_predictions = best_model.predict(final_df[FEATURE_NAMES]).tolist()
        final_truth = final_df["true_completeness"].tolist()
        final_mae, final_r2 = _metrics(final_truth, final_predictions)
        baseline_mae, baseline_r2 = _metrics(
            final_truth, final_df["strategy1_score"].astype(float).tolist()
        )
        models[order_name] = best_model
        output_tiers.update(development_tiers)
        output_baselines.update(development_baselines)
        output_profiles.update(development_profiles)
        for role, genome_ids in (
            ("development", development_ids),
            ("final_test", plan["final_test_ids"]),
        ):
            for genome_id in genome_ids:
                split_rows.append(
                    {
                        "order_name": order_name,
                        "stage": "final_evaluation",
                        "fold": "final",
                        "role": role,
                        "split_unit": str(plan["split_unit"]),
                        "genome_id": genome_id,
                        "family_name": labels.get(genome_id, {}).get("family", ""),
                    }
                )
        metadata_rows.append(
            {
                **common_metadata,
                "training_status": "trained_marker_dropout_only",
                "unsupported_reason": "",
                "model_name": best_name,
                "validation_type": f"independent_{plan['split_unit']}_holdout_marker_dropout",
                "selection_validation_type": f"grouped_{plan['split_unit']}_folds",
                "n_selection_folds": len(fold_frames),
                "n_development_genomes": len(development_ids),
                "n_final_test_genomes": len(plan["final_test_ids"]),
                "n_training_samples": len(development_df),
                "n_training_samples_used": len(development_df),
                "selection_marker_dropout_mae": round(float(best_mae), 4),
                "selection_marker_dropout_r2": round(float(best_r2), 4),
                "marker_dropout_final_test_mae": round(final_mae, 4),
                "marker_dropout_final_test_r2": round(final_r2, 4),
                "test_r2": round(final_r2, 4),
                "test_mae": round(final_mae, 4),
                "test_metric_scope": "marker_dropout_only_not_nucleotide_completeness",
                "marker_recovery_baseline_final_test_mae": round(baseline_mae, 4),
                "marker_recovery_baseline_final_test_r2": round(baseline_r2, 4),
            }
        )
    return (
        models,
        pd.DataFrame(metadata_rows),
        output_tiers,
        output_baselines,
        output_profiles,
        split_rows,
    )


def compute_marker_conservation(
    order_markers: Dict[str, List[str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
) -> List[Dict[str, object]]:
    """Per-(order, marker) conservation used by the weighted-completeness scorer.

    Reuses the same prevalence definition (present / genome_count) that the
    strategy-2 tier sets are derived from, so ``weighted_order_completeness``
    and the tier model agree on marker conservation rather than the scorer
    silently falling back to uniform 1.0 weights (the file the scorer used to
    read, ``markers/stats.tsv``, is an HMM size table with the wrong schema).
    """
    rows: List[Dict[str, object]] = []
    for order_name, markers in sorted(order_markers.items()):
        genome_profiles = order_profiles.get(order_name, {})
        stats = marker_stats_for_group(genome_profiles, markers)
        if stats.empty:
            continue
        for _, row in stats.iterrows():
            rows.append(
                {
                    "marker_name": row["marker"],
                    "order_name": order_name,
                    "percent_genomes_with_marker": round(
                        float(row["prevalence"]) * 100.0, 4
                    ),
                    "avg_duplication_factor": round(float(row["avg_copy"]), 4),
                    "n_genomes": int(row["n_genomes"]),
                }
            )
    return rows


def write_marker_conservation(
    resources_dir: Path, rows: List[Dict[str, object]]
) -> None:
    path = resources_dir / "markers" / "marker_conservation.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "marker_name",
        "order_name",
        "percent_genomes_with_marker",
        "avg_duplication_factor",
        "n_genomes",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_order_table(
    resources_dir: Path,
    rows: List[Dict[str, str]],
    fieldnames: List[str],
    order_markers: Dict[str, List[str]],
    development_profiles: Dict[str, Dict[str, Dict[str, int]]],
    model_metadata: pd.DataFrame,
) -> None:
    status_by_order = {
        str(row["order_name"]): str(row["training_status"])
        for _, row in model_metadata.iterrows()
    }
    output_fields = list(fieldnames)
    for name in ("Average_Percent", "Std_Percent", "Builder_Status", "Baseline_Source"):
        if name not in output_fields:
            output_fields.append(name)
    path = resources_dir / "markers" / "order_completeness.tab"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        for source_row in rows:
            row = dict(source_row)
            order_name = row["Order"]
            profiles = development_profiles.get(order_name, {})
            status = status_by_order.get(order_name, "unsupported")
            row["Builder_Status"] = status
            if profiles:
                baseline_mean, baseline_std = order_recovery_stats(
                    profiles, order_markers[order_name]
                )
                row["Average_Percent"] = str(round(baseline_mean, 8))
                row["Std_Percent"] = str(round(baseline_std, 8))
                if "n_genomes" in row:
                    row["n_genomes"] = str(len(profiles))
                row["Baseline_Source"] = (
                    "development_only_reference_profiles"
                    if status == "trained_marker_dropout_only"
                    else "unvalidated_reference_profiles_no_ml"
                )
            else:
                row["Baseline_Source"] = "source_table_unvalidated_unsupported_group"
            writer.writerow(row)


def write_resources(
    resources_dir: Path,
    order_markers: Dict[str, List[str]],
    tiers: Dict[Tuple[str, str, str], TierSet],
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
    models: Dict[str, object],
    model_metadata: pd.DataFrame,
) -> None:
    (resources_dir / "completeness").mkdir(parents=True, exist_ok=True)
    config_path = resources_dir / "completeness/config.json"
    config_payload = {
        "resource_version": "v2_marker_dropout",
        "strategy2_config": asdict(CONFIG),
        "order_expected_marker_counts": {order: len(markers) for order, markers in order_markers.items()},
        "training_target": "marker_dropout_fraction_of_source_profile",
        "target_is_nucleotide_completeness": False,
        "reference_ascertainment": "labeled_references_with_detected_order_markers",
        "marker_panel_ascertainment": "fixed_input_panel_predates_validation_holdouts",
    }
    config_path.write_text(json.dumps(config_payload, indent=2, sort_keys=True) + "\n")

    tiers_path = resources_dir / "completeness/tiers.tsv"
    with tiers_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "group_level",
                "group_name",
                "order_name",
                "tier",
                "marker_name",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for (group_level, group_name, order_name), tier_set in sorted(tiers.items()):
            for tier_name, markers in (
                ("core", tier_set.core),
                ("shared", tier_set.shared),
                ("accessory", tier_set.accessory),
            ):
                for marker_name in markers:
                    writer.writerow(
                        {
                            "group_level": group_level,
                            "group_name": group_name,
                            "order_name": order_name,
                            "tier": tier_name,
                            "marker_name": marker_name,
                        }
                    )

    baselines_path = resources_dir / "completeness/baselines.tsv"
    with baselines_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "group_level",
                "group_name",
                "order_name",
                "baseline_mean",
                "baseline_median",
                "n_refs",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for (group_level, group_name, order_name), baseline in sorted(
            baselines.items()
        ):
            writer.writerow(
                {
                    "group_level": group_level,
                    "group_name": group_name,
                    "order_name": order_name,
                    **baseline,
                }
            )

    model_metadata.to_csv(
        resources_dir / "completeness/model_metadata.tsv", sep="\t", index=False
    )
    dump(
        {
            "models": models,
            "feature_names": FEATURE_NAMES,
            "training_target": "marker_dropout_fraction_of_source_profile",
            "target_is_nucleotide_completeness": False,
        },
        resources_dir / "completeness/model.joblib",
        compress=3,
    )


def write_validation_splits(resources_dir: Path, rows: List[Dict[str, str]]) -> None:
    path = resources_dir / "completeness" / "validation_splits.tsv"
    with path.open("w", newline="") as handle:
        fieldnames = [
            "order_name",
            "stage",
            "fold",
            "role",
            "split_unit",
            "genome_id",
            "family_name",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_resource_bundle(
    resources_dir: Path, output_dir: Path, smoke: bool = False
) -> Tuple[int, int]:
    resources_dir = resources_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir == resources_dir:
        raise ValueError("--output-dir must differ from --resources-dir")
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")

    old_cache = os.environ.get("GVCLASS_RESOURCE_CACHE")
    with tempfile.TemporaryDirectory(prefix="gvclass_resource_builder_cache_") as cache:
        os.environ["GVCLASS_RESOURCE_CACHE"] = cache
        try:
            order_markers, order_rows, order_fieldnames = load_order_table(
                resources_dir
            )
            labels = load_labels(resources_dir)
            order_profiles = parse_reference_marker_profiles(
                resources_dir, labels, order_markers
            )
        finally:
            if old_cache is None:
                os.environ.pop("GVCLASS_RESOURCE_CACHE", None)
            else:
                os.environ["GVCLASS_RESOURCE_CACHE"] = old_cache

    profile_count = sum(len(profiles) for profiles in order_profiles.values())
    if not order_markers:
        raise RuntimeError("No order marker definitions were loaded")
    if not labels:
        raise RuntimeError("No reference labels were loaded")
    if profile_count == 0:
        raise RuntimeError("No reference marker profiles were loaded")

    (
        models,
        model_metadata,
        tiers,
        baselines,
        development_profiles,
        split_rows,
    ) = train_models(labels, order_profiles, order_markers, smoke=smoke)
    if not models or not tiers or not baselines:
        raise RuntimeError(
            "No model group has enough independent source genomes for grouped "
            "selection and a final test"
        )

    marker_conservation = compute_marker_conservation(
        order_markers, development_profiles
    )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        write_marker_conservation(staging_dir, marker_conservation)
        write_order_table(
            staging_dir,
            order_rows,
            order_fieldnames,
            order_markers,
            development_profiles,
            model_metadata,
        )
        write_resources(
            staging_dir,
            order_markers,
            tiers,
            baselines,
            models,
            model_metadata,
        )
        write_validation_splits(staging_dir, split_rows)
        if output_dir.exists():
            raise FileExistsError(f"Output directory already exists: {output_dir}")
        os.replace(staging_dir, output_dir)
    except Exception:
        shutil.rmtree(staging_dir, ignore_errors=True)
        raise
    return len(tiers), len(models)


def main() -> None:
    args = parse_args()
    resources_dir = Path(args.resources_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    group_count, model_count = build_resource_bundle(
        resources_dir, output_dir, smoke=args.smoke
    )
    print(f"output_dir\t{output_dir}")
    print(f"strategy2_groups\t{group_count}")
    print(f"strategy3_models\t{model_count}")
    model_digest = hashlib.sha256(
        (output_dir / "completeness" / "model.joblib").read_bytes()
    ).hexdigest()
    print(f"model_sha256_for_manual_review\t{model_digest}")


if __name__ == "__main__":
    main()
