"""Build novelty-aware completeness resources for GVClass.

This is an offline resource-generation step. It derives strategy-2 tier sets and
baselines from the reference database and trains the novelty-aware strategy-3
models with leave-one-family-out validation where possible.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from joblib import dump
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split

from src.core.novelty_completeness import GroupDefinition, Strategy2Config, TierSet

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
    return parser.parse_args()


def load_order_table(resources_dir: Path) -> Tuple[Dict[str, List[str]], Dict[str, float]]:
    order_markers: Dict[str, List[str]] = {}
    order_baselines: Dict[str, float] = {}
    path = resources_dir / "order_completeness.tab"
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            order_name = row["Order"]
            order_markers[order_name] = [
                marker.strip() for marker in row["Orthogroups"].split(",") if marker.strip()
            ]
            order_baselines[order_name] = float(row["Average_Percent"])
    return order_markers, order_baselines


def load_labels(resources_dir: Path) -> Dict[str, Dict[str, str]]:
    labels: Dict[str, Dict[str, str]] = {}
    path = resources_dir / "gvclassFeb26_labels.tsv"
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
) -> Tuple[
    Dict[str, Dict[str, Dict[str, int]]],
    Dict[Tuple[str, str], Dict[str, Dict[str, int]]],
]:
    relevant_markers = sorted({marker for markers in order_markers.values() for marker in markers})
    order_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    family_profiles: Dict[Tuple[str, str], Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    faa_dir = resources_dir / "database" / "faa"

    for marker in relevant_markers:
        marker_path = faa_dir / f"{marker}.faa"
        if not marker_path.exists():
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
            order_name = tax["order"]
            family_name = tax["family"]
            if order_name not in order_markers:
                continue
            order_profiles[order_name][genome_id][marker] = copy_count
            family_profiles[(order_name, family_name)][genome_id][marker] = copy_count

    return order_profiles, family_profiles


def marker_stats_for_group(genome_profiles: Dict[str, Dict[str, int]], markers: List[str]) -> pd.DataFrame:
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


def score_tier_set(marker_counts: Dict[str, int], tier_set: TierSet) -> Dict[str, float]:
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
        family_stats = marker_stats_for_group(genome_profiles, order_markers[order_name])
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
) -> List[Dict[str, object]]:
    rng = random.Random(RANDOM_SEED)
    samples: List[Dict[str, object]] = []
    for order_name, genome_profiles in order_profiles.items():
        markers = order_markers[order_name]
        for genome_id, profile in genome_profiles.items():
            present_markers = [marker for marker in markers if profile.get(marker, 0) > 0]
            if not present_markers:
                continue
            family_name = labels.get(genome_id, {}).get("family", "")
            for completeness_level in SIMULATED_COMPLETENESS_LEVELS:
                retain_n = max(1, round(len(present_markers) * completeness_level))
                for rep in range(SIMULATIONS_PER_LEVEL_FINAL):
                    retained = set(rng.sample(present_markers, retain_n))
                    samples.append(
                        {
                            "order_name": order_name,
                            "family_name": family_name,
                            "true_completeness": completeness_level * 100.0,
                            "present_markers": retained,
                            "rep": rep,
                        }
                    )
    return samples


def downsample_training_order_df(order_df: pd.DataFrame) -> pd.DataFrame:
    family_names = [family for family in sorted(order_df["family_name"].unique()) if family]
    if len(family_names) >= 2:
        parts = []
        for _, family_df in order_df.groupby("family_name", sort=False):
            if len(family_df) > MAX_SYNTHETIC_SAMPLES_PER_FAMILY:
                parts.append(family_df.sample(n=MAX_SYNTHETIC_SAMPLES_PER_FAMILY, random_state=RANDOM_SEED))
            else:
                parts.append(family_df)
        return pd.concat(parts, ignore_index=True)
    if len(order_df) > MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY:
        return order_df.sample(n=MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY, random_state=RANDOM_SEED).reset_index(drop=True)
    return order_df.reset_index(drop=True)


def build_training_df(
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    order_baselines: Dict[str, float],
    tiers: Dict[Tuple[str, str, str], TierSet],
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for sample in simulate_samples(labels, order_profiles, order_markers):
        order_name = str(sample["order_name"])
        family_name = str(sample["family_name"])
        group_def, tier_set = select_group(tiers, order_name, family_name)
        raw_score = 100.0 * (len(sample["present_markers"]) / len(order_markers[order_name]))
        baseline = order_baselines[order_name]
        strategy1_score = 0.0 if baseline <= 0 else min(100.0, (raw_score / baseline) * 100.0)
        marker_counts = {marker: 1 for marker in sample["present_markers"]}
        tier_scores = score_tier_set(marker_counts, tier_set)
        strategy2_raw_score = float(tier_scores["strategy2_raw"])
        strategy2_norm_score, strategy2_baseline_mean, strategy2_baseline_n_refs = normalize_strategy2_score(
            strategy2_raw_score,
            group_def,
            baselines,
        )
        rows.append(
            {
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
            }
        )
    return pd.DataFrame(rows)


def train_models(training_df: pd.DataFrame) -> Tuple[Dict[str, object], pd.DataFrame]:
    models: Dict[str, object] = {}
    metadata_rows: List[Dict[str, object]] = []
    for order_name, order_df in training_df.groupby("order_name"):
        model_order_df = downsample_training_order_df(order_df)
        if len(model_order_df) < 20:
            continue
        X = model_order_df[FEATURE_NAMES]
        y = model_order_df["true_completeness"]
        family_names = [family for family in sorted(model_order_df["family_name"].unique()) if family]
        use_lofo = len(family_names) >= 2

        candidate_models = [
            (
                "random_forest",
                lambda: RandomForestRegressor(n_estimators=400, random_state=RANDOM_SEED, min_samples_leaf=2),
            ),
            (
                "extra_trees",
                lambda: ExtraTreesRegressor(n_estimators=500, random_state=RANDOM_SEED, min_samples_leaf=2),
            ),
            (
                "hist_gbm",
                lambda: HistGradientBoostingRegressor(random_state=RANDOM_SEED, max_depth=8, min_samples_leaf=10),
            ),
        ]

        best_model = None
        best_model_name = ""
        best_mae = None
        best_r2 = None
        best_validation_type = ""
        best_validation_families = 0

        for model_name, model_factory in candidate_models:
            preds: List[float] = []
            truth: List[float] = []
            validation_type = "random_split"
            validation_families = 0
            if use_lofo:
                validation_type = "leave_one_family_out"
                for family_name in family_names:
                    train_mask = model_order_df["family_name"] != family_name
                    test_mask = model_order_df["family_name"] == family_name
                    if not train_mask.any() or not test_mask.any():
                        continue
                    model = model_factory()
                    model.fit(model_order_df.loc[train_mask, FEATURE_NAMES], model_order_df.loc[train_mask, "true_completeness"])
                    pred = model.predict(model_order_df.loc[test_mask, FEATURE_NAMES])
                    preds.extend(pred.tolist())
                    truth.extend(model_order_df.loc[test_mask, "true_completeness"].tolist())
                    validation_families += 1
            if not preds:
                X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=RANDOM_SEED)
                model = model_factory()
                model.fit(X_train, y_train)
                pred = model.predict(X_test)
                preds = pred.tolist()
                truth = y_test.tolist()
                validation_type = "random_split"

            mae = float(mean_absolute_error(truth, preds))
            r2 = float(r2_score(truth, preds))
            if best_mae is None or mae < best_mae or (mae == best_mae and r2 > best_r2):
                best_model = model_factory()
                best_model.fit(X, y)
                best_model_name = model_name
                best_mae = mae
                best_r2 = r2
                best_validation_type = validation_type
                best_validation_families = validation_families

        models[order_name] = best_model
        metadata_rows.append(
            {
                "order_name": order_name,
                "n_training_samples": len(order_df),
                "n_training_samples_used": len(model_order_df),
                "model_name": best_model_name,
                "validation_type": best_validation_type,
                "n_validation_families": best_validation_families,
                "test_r2": round(float(best_r2), 4),
                "test_mae": round(float(best_mae), 4),
            }
        )
    return models, pd.DataFrame(metadata_rows)


def write_resources(
    resources_dir: Path,
    order_markers: Dict[str, List[str]],
    tiers: Dict[Tuple[str, str, str], TierSet],
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
    models: Dict[str, object],
    model_metadata: pd.DataFrame,
) -> None:
    config_path = resources_dir / "novelty_completeness_config.json"
    config_payload = {
        "resource_version": "v1",
        "strategy2_config": asdict(CONFIG),
        "order_expected_marker_counts": {order: len(markers) for order, markers in order_markers.items()},
    }
    config_path.write_text(json.dumps(config_payload, indent=2, sort_keys=True) + "\n")

    tiers_path = resources_dir / "novelty_strategy2_tiers.tsv"
    with tiers_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["group_level", "group_name", "order_name", "tier", "marker_name"],
            delimiter="\t",
        )
        writer.writeheader()
        for (group_level, group_name, order_name), tier_set in sorted(tiers.items()):
            for tier_name, markers in (("core", tier_set.core), ("shared", tier_set.shared), ("accessory", tier_set.accessory)):
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

    baselines_path = resources_dir / "novelty_strategy2_baselines.tsv"
    with baselines_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["group_level", "group_name", "order_name", "baseline_mean", "baseline_median", "n_refs"],
            delimiter="\t",
        )
        writer.writeheader()
        for (group_level, group_name, order_name), baseline in sorted(baselines.items()):
            writer.writerow(
                {
                    "group_level": group_level,
                    "group_name": group_name,
                    "order_name": order_name,
                    **baseline,
                }
            )

    model_metadata.to_csv(resources_dir / "novelty_strategy3_model_metadata.tsv", sep="\t", index=False)
    dump(
        {"models": models, "feature_names": FEATURE_NAMES},
        resources_dir / "novelty_strategy3_models.joblib",
        compress=3,
    )


def main() -> None:
    args = parse_args()
    resources_dir = Path(args.resources_dir).resolve()
    order_markers, order_baselines = load_order_table(resources_dir)
    labels = load_labels(resources_dir)
    order_profiles, family_profiles = parse_reference_marker_profiles(resources_dir, labels, order_markers)
    tiers, baselines = build_group_resources(order_markers, order_profiles, family_profiles)
    training_df = build_training_df(labels, order_profiles, order_markers, order_baselines, tiers, baselines)
    models, model_metadata = train_models(training_df)
    write_resources(resources_dir, order_markers, tiers, baselines, models, model_metadata)
    print(f"strategy2_groups\t{len(tiers)}")
    print(f"strategy3_models\t{len(models)}")


if __name__ == "__main__":
    main()
