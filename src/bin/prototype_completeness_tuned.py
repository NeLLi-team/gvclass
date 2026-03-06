"""Tune prototype completeness strategies 2 and 3 on an existing GVClass run.

This script does three things:
1. Grid-searches strategy 2 configurations against the curated isolate set and
   synthetic degraded reference genomes.
2. Selects the best practical strategy 2 configuration.
3. Retrains strategy 3 with the improved strategy 2 features and writes
   side-by-side comparison tables.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from sklearn.ensemble import (
    ExtraTreesRegressor,
    HistGradientBoostingRegressor,
    RandomForestRegressor,
)
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split


RANDOM_SEED = 42
SIMULATED_COMPLETENESS_LEVELS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
SIMULATIONS_PER_LEVEL_TUNING = 2
SIMULATIONS_PER_LEVEL_FINAL = 3
MAX_SYNTHETIC_SAMPLES_PER_FAMILY = 2500
MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY = 8000


@dataclass(frozen=True)
class GroupDefinition:
    level: str
    name: str
    order_name: str


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
class TierSet:
    core: Tuple[str, ...]
    shared: Tuple[str, ...]
    accessory: Tuple[str, ...]


OOD_STRICT_FLAGS = {
    "unassigned",
    "no_informative_markers",
    "zero_informative_hits",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10")
    parser.add_argument("--resources-dir", default="resources")
    parser.add_argument(
        "--best-config-json",
        default="",
        help="Optional path to a Strategy 2 config JSON. If set, skip grid tuning.",
    )
    parser.add_argument(
        "--output-suffix",
        default="tuned_iter2",
        help="Suffix for output files from this run",
    )
    return parser.parse_args()


def load_strategy1_summary(run_dir: Path) -> List[Dict[str, str]]:
    path = run_dir / "gvclass_extended_results" / "gvclass_summary_strategy1.tsv"
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_metadata(run_dir: Path) -> Dict[str, Dict[str, str]]:
    path = run_dir / "metadata_table.tsv"
    metadata: Dict[str, Dict[str, str]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            taxonomy_parts = [part.strip() for part in row["taxonomy"].split(";")]
            metadata[row["accession_number"]] = {
                **row,
                "ncbi_family": taxonomy_parts[6] if len(taxonomy_parts) > 6 else "NA",
                "ncbi_genus": row.get("ncbi_genus", "NA") or "NA",
            }
    return metadata


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


def extract_primary_token(value: str) -> str:
    if not value:
        return ""
    token = value.split(",")[0].split(":")[0]
    return token.split("__", 1)[1] if "__" in token else token


def build_query_counts(run_dir: Path, accessions: Iterable[str]) -> Dict[str, Dict[str, int]]:
    results_dir = run_dir / "gvclass_extended_results"
    counts_by_accession: Dict[str, Dict[str, int]] = {}
    for accession in accessions:
        tar_path = results_dir / f"{accession}.tar.gz"
        counts: Dict[str, int] = {}
        with tarfile.open(tar_path, "r:gz") as tar_handle:
            member = tar_handle.extractfile(f"{accession}/hmmout/models.counts")
            if member is None:
                counts_by_accession[accession] = counts
                continue
            for line in member.read().decode().splitlines():
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                try:
                    counts[parts[0]] = int(parts[1])
                except ValueError:
                    continue
        counts_by_accession[accession] = counts
    return counts_by_accession


def collect_relevant_orders(summary_rows: List[Dict[str, str]]) -> List[str]:
    return sorted(
        {
            extract_primary_token(row["order"])
            for row in summary_rows
            if extract_primary_token(row["order"])
        }
    )


def parse_reference_marker_profiles(
    resources_dir: Path,
    labels: Dict[str, Dict[str, str]],
    order_markers: Dict[str, List[str]],
    relevant_orders: Iterable[str],
) -> Tuple[
    Dict[str, Dict[str, Dict[str, int]]],
    Dict[str, Dict[str, Dict[str, int]]],
]:
    relevant_order_set = set(relevant_orders)
    relevant_markers = sorted(
        {
            marker
            for order_name in relevant_order_set
            for marker in order_markers.get(order_name, [])
        }
    )

    order_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    family_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(lambda: defaultdict(dict))
    faa_dir = resources_dir / "database" / "faa"

    for marker in relevant_markers:
        marker_path = faa_dir / f"{marker}.faa"
        genome_copy_counts: Dict[str, int] = Counter()
        if not marker_path.exists():
            continue
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
            if order_name not in relevant_order_set:
                continue
            order_profiles[order_name][genome_id][marker] = copy_count
            family_profiles[family_name][genome_id][marker] = copy_count

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


def build_stats_cache(
    relevant_orders: List[str],
    order_markers: Dict[str, List[str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
) -> Tuple[Dict[str, pd.DataFrame], Dict[Tuple[str, str], pd.DataFrame]]:
    order_stats_cache: Dict[str, pd.DataFrame] = {}
    family_stats_cache: Dict[Tuple[str, str], pd.DataFrame] = {}
    for order_name in relevant_orders:
        order_stats_cache[order_name] = marker_stats_for_group(
            order_profiles.get(order_name, {}),
            order_markers.get(order_name, []),
        )
    for family_name, genome_profiles in family_profiles.items():
        orders = set()
        for genome_id in genome_profiles:
            orders.add(next((o for o, prof in order_profiles.items() if genome_id in prof), None))
        orders.discard(None)
        if len(orders) != 1:
            continue
        order_name = next(iter(orders))
        family_stats_cache[(order_name, family_name)] = marker_stats_for_group(
            genome_profiles,
            order_markers.get(order_name, []),
        )
    return order_stats_cache, family_stats_cache


def build_tier_set(group_stats: pd.DataFrame, config: Strategy2Config) -> TierSet:
    if group_stats.empty:
        return TierSet((), (), ())

    core = sorted(
        group_stats[
            (group_stats["prevalence"] >= config.core_prevalence)
            & (group_stats["avg_copy"] <= config.core_max_copy)
        ]["marker"].tolist()
    )
    shared = sorted(
        group_stats[
            (group_stats["prevalence"] >= config.shared_prevalence)
            & (group_stats["avg_copy"] <= config.shared_max_copy)
            & (group_stats["prevalence"] < config.core_prevalence)
        ]["marker"].tolist()
    )
    accessory = sorted(
        group_stats[
            (group_stats["prevalence"] >= config.accessory_prevalence)
            & (group_stats["prevalence"] < config.shared_prevalence)
        ]["marker"].tolist()
    )
    return TierSet(tuple(core), tuple(shared), tuple(accessory))


def informative_marker_count(tier_set: TierSet) -> int:
    return len(tier_set.core) + len(tier_set.shared)


def select_group_and_tiers(
    order_name: str,
    family_name: str,
    config: Strategy2Config,
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_stats_cache: Dict[str, pd.DataFrame],
    family_stats_cache: Dict[Tuple[str, str], pd.DataFrame],
) -> Tuple[GroupDefinition, TierSet]:
    family_profile = family_profiles.get(family_name, {})
    family_stats = family_stats_cache.get((order_name, family_name), pd.DataFrame())
    family_tiers = build_tier_set(family_stats, config)
    if (
        len(family_profile) >= config.family_min_ref_genomes
        and informative_marker_count(family_tiers) >= config.family_min_informative_markers
    ):
        return GroupDefinition("family", family_name, order_name), family_tiers

    order_stats = order_stats_cache.get(order_name, pd.DataFrame())
    order_tiers = build_tier_set(order_stats, config)
    return GroupDefinition("order", order_name, order_name), order_tiers


def presence_fraction(present_markers: set[str], markers: Tuple[str, ...]) -> float:
    if not markers:
        return 0.0
    present = sum(1 for marker in markers if marker in present_markers)
    return present / len(markers)


def score_tier_set(present_markers: set[str], tier_set: TierSet, config: Strategy2Config) -> Dict[str, float]:
    core_fraction = presence_fraction(present_markers, tier_set.core)
    shared_fraction = presence_fraction(present_markers, tier_set.shared)
    accessory_fraction = presence_fraction(present_markers, tier_set.accessory)

    used_weights = []
    weighted_values = []
    if tier_set.core:
        used_weights.append(config.weight_core)
        weighted_values.append(config.weight_core * core_fraction)
    if tier_set.shared:
        used_weights.append(config.weight_shared)
        weighted_values.append(config.weight_shared * shared_fraction)
    if tier_set.accessory:
        used_weights.append(config.weight_accessory)
        weighted_values.append(config.weight_accessory * accessory_fraction)

    if not used_weights:
        strategy2_score = 0.0
    else:
        strategy2_score = 100.0 * (sum(weighted_values) / sum(used_weights))

    return {
        "strategy2": round(strategy2_score, 2),
        "core_fraction": round(core_fraction * 100.0, 2),
        "shared_fraction": round(shared_fraction * 100.0, 2),
        "accessory_fraction": round(accessory_fraction * 100.0, 2),
        "core_marker_count": len(tier_set.core),
        "shared_marker_count": len(tier_set.shared),
        "accessory_marker_count": len(tier_set.accessory),
        "informative_marker_count": informative_marker_count(tier_set),
    }


def simulate_samples(
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    reps: int,
) -> List[Dict[str, object]]:
    rng = random.Random(RANDOM_SEED)
    samples: List[Dict[str, object]] = []
    for order_name, genome_profiles in order_profiles.items():
        markers = order_markers.get(order_name, [])
        if not markers:
            continue
        for genome_id, profile in genome_profiles.items():
            present_markers = [marker for marker in markers if profile.get(marker, 0) > 0]
            if not present_markers:
                continue
            family_name = labels.get(genome_id, {}).get("family", "")
            for completeness_level in SIMULATED_COMPLETENESS_LEVELS:
                retain_n = max(1, round(len(present_markers) * completeness_level))
                for rep in range(reps):
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


def build_config_grid() -> List[Strategy2Config]:
    configs: List[Strategy2Config] = []
    idx = 1
    for core_prev in (0.85, 0.9):
        for shared_prev in (0.65, 0.75):
            for accessory_prev in (0.35, 0.5):
                if not (accessory_prev < shared_prev < core_prev):
                    continue
                for core_max_copy in (1.2, 1.5):
                    for shared_max_copy in (1.5, 2.0):
                        if shared_max_copy < core_max_copy:
                            continue
                        for family_min_ref in (3, 5):
                            for family_min_info in (2, 3):
                                for weights in ((0.7, 0.2, 0.1), (0.6, 0.3, 0.1), (0.5, 0.3, 0.2)):
                                    configs.append(
                                        Strategy2Config(
                                            name=f"cfg_{idx:03d}",
                                            core_prevalence=core_prev,
                                            shared_prevalence=shared_prev,
                                            accessory_prevalence=accessory_prev,
                                            core_max_copy=core_max_copy,
                                            shared_max_copy=shared_max_copy,
                                            family_min_ref_genomes=family_min_ref,
                                            family_min_informative_markers=family_min_info,
                                            weight_core=weights[0],
                                            weight_shared=weights[1],
                                            weight_accessory=weights[2],
                                        )
                                    )
                                    idx += 1
    return configs


def load_config_from_json(path: Path) -> Strategy2Config:
    data = json.loads(path.read_text())
    return Strategy2Config(
        name=str(data["name"]),
        core_prevalence=float(data["core_prevalence"]),
        shared_prevalence=float(data["shared_prevalence"]),
        accessory_prevalence=float(data["accessory_prevalence"]),
        core_max_copy=float(data["core_max_copy"]),
        shared_max_copy=float(data["shared_max_copy"]),
        family_min_ref_genomes=int(data["family_min_ref_genomes"]),
        family_min_informative_markers=int(data["family_min_informative_markers"]),
        weight_core=float(data["weight_core"]),
        weight_shared=float(data["weight_shared"]),
        weight_accessory=float(data["weight_accessory"]),
    )


def build_strategy2_baselines(
    config: Strategy2Config,
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_stats_cache: Dict[str, pd.DataFrame],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_stats_cache: Dict[Tuple[str, str], pd.DataFrame],
) -> Dict[Tuple[str, str, str], Dict[str, float]]:
    tier_cache: Dict[Tuple[str, str], Tuple[GroupDefinition, TierSet]] = {}
    group_scores: Dict[Tuple[str, str, str], List[float]] = defaultdict(list)

    def get_group(order_name: str, family_name: str) -> Tuple[GroupDefinition, TierSet]:
        key = (order_name, family_name)
        if key not in tier_cache:
            tier_cache[key] = select_group_and_tiers(
                order_name,
                family_name,
                config,
                order_profiles,
                family_profiles,
                order_stats_cache,
                family_stats_cache,
            )
        return tier_cache[key]

    for order_name, genome_profiles in order_profiles.items():
        for genome_id, profile in genome_profiles.items():
            family_name = labels.get(genome_id, {}).get("family", "")
            group_def, tier_set = get_group(order_name, family_name)
            present_markers = {marker for marker, count in profile.items() if count > 0}
            raw_score = float(score_tier_set(present_markers, tier_set, config)["strategy2"])
            group_key = (group_def.level, group_def.name, order_name)
            group_scores[group_key].append(raw_score)

    baselines: Dict[Tuple[str, str, str], Dict[str, float]] = {}
    for group_key, scores in group_scores.items():
        baselines[group_key] = {
            "baseline_mean": round(mean(scores), 4),
            "baseline_median": round(median(scores), 4),
            "n_refs": len(scores),
        }
    return baselines


def normalize_strategy2_score(
    raw_score: float,
    group_def: GroupDefinition,
    baselines: Dict[Tuple[str, str, str], Dict[str, float]],
) -> Tuple[float, float, int]:
    group_key = (group_def.level, group_def.name, group_def.order_name)
    baseline = baselines.get(group_key, {})
    baseline_mean = float(baseline.get("baseline_mean", 0.0))
    n_refs = int(baseline.get("n_refs", 0))
    if baseline_mean <= 0:
        return raw_score, baseline_mean, n_refs
    normalized = min(100.0, (raw_score / baseline_mean) * 100.0)
    return round(normalized, 2), baseline_mean, n_refs


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


def evaluate_config(
    config: Strategy2Config,
    summary_rows: List[Dict[str, str]],
    metadata: Dict[str, Dict[str, str]],
    query_counts: Dict[str, Dict[str, int]],
    order_markers: Dict[str, List[str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_stats_cache: Dict[str, pd.DataFrame],
    family_stats_cache: Dict[Tuple[str, str], pd.DataFrame],
    synthetic_samples: List[Dict[str, object]],
) -> Dict[str, object]:
    isolate_scores: List[float] = []
    isolate_lt50 = 0
    isolate_ge90 = 0
    informative_zero = 0
    group_levels = Counter()

    tier_cache: Dict[Tuple[str, str], Tuple[GroupDefinition, TierSet]] = {}

    def get_group(order_name: str, family_name: str) -> Tuple[GroupDefinition, TierSet]:
        key = (order_name, family_name)
        if key not in tier_cache:
            tier_cache[key] = select_group_and_tiers(
                order_name,
                family_name,
                config,
                order_profiles,
                family_profiles,
                order_stats_cache,
                family_stats_cache,
            )
        return tier_cache[key]

    for row in summary_rows:
        order_token = extract_primary_token(row["order"])
        family_token = extract_primary_token(row["family"])
        accession = row["query"]
        if order_token not in order_markers:
            score = float(row["order_completeness"])
            group_levels["unavailable"] += 1
        else:
            group_def, tier_set = get_group(order_token, family_token)
            group_levels[group_def.level] += 1
            if informative_marker_count(tier_set) == 0:
                informative_zero += 1
            present_markers = {m for m, c in query_counts[accession].items() if c > 0}
            score = score_tier_set(present_markers, tier_set, config)["strategy2"]
        isolate_scores.append(score)
        if score < 50.0:
            isolate_lt50 += 1
        if score >= 90.0:
            isolate_ge90 += 1

    synthetic_scores: List[float] = []
    synthetic_truth: List[float] = []
    for sample in synthetic_samples:
        order_name = str(sample["order_name"])
        family_name = str(sample["family_name"])
        if order_name not in order_markers:
            continue
        _, tier_set = get_group(order_name, family_name)
        score = score_tier_set(sample["present_markers"], tier_set, config)["strategy2"]
        synthetic_scores.append(score)
        synthetic_truth.append(float(sample["true_completeness"]))

    synthetic_mae = mean_absolute_error(synthetic_truth, synthetic_scores) if synthetic_scores else 100.0
    isolate_mae_100 = mean(abs(100.0 - score) for score in isolate_scores)
    composite = synthetic_mae + (0.5 * isolate_mae_100) + (0.2 * informative_zero)

    return {
        **asdict(config),
        "n_genomes": len(isolate_scores),
        "isolate_mean": round(mean(isolate_scores), 2),
        "isolate_median": round(median(isolate_scores), 2),
        "isolate_ge90": isolate_ge90,
        "isolate_lt50": isolate_lt50,
        "isolate_mae_to_100": round(isolate_mae_100, 2),
        "synthetic_mae": round(float(synthetic_mae), 2),
        "synthetic_r2": round(float(r2_score(synthetic_truth, synthetic_scores)), 4),
        "informative_zero_groups": informative_zero,
        "group_level_family": group_levels.get("family", 0),
        "group_level_order": group_levels.get("order", 0),
        "group_level_unavailable": group_levels.get("unavailable", 0),
        "composite_rank_score": round(float(composite), 2),
    }


def build_strategy3_training_table(
    config: Strategy2Config,
    labels: Dict[str, Dict[str, str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    order_baselines: Dict[str, float],
    order_stats_cache: Dict[str, pd.DataFrame],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_stats_cache: Dict[Tuple[str, str], pd.DataFrame],
    strategy2_baselines: Dict[Tuple[str, str, str], Dict[str, float]],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    synthetic_samples = simulate_samples(labels, order_profiles, order_markers, SIMULATIONS_PER_LEVEL_FINAL)
    tier_cache: Dict[Tuple[str, str], Tuple[GroupDefinition, TierSet]] = {}

    def get_group(order_name: str, family_name: str) -> Tuple[GroupDefinition, TierSet]:
        key = (order_name, family_name)
        if key not in tier_cache:
            tier_cache[key] = select_group_and_tiers(
                order_name,
                family_name,
                config,
                order_profiles,
                family_profiles,
                order_stats_cache,
                family_stats_cache,
            )
        return tier_cache[key]

    training_rows: List[Dict[str, float | str | int]] = []
    for sample in synthetic_samples:
        order_name = str(sample["order_name"])
        family_name = str(sample["family_name"])
        if order_name not in order_markers:
            continue
        group_def, tier_set = get_group(order_name, family_name)
        raw_score = 100.0 * (len(sample["present_markers"]) / len(order_markers[order_name]))
        baseline = order_baselines.get(order_name, 0.0)
        strategy1_score = 0.0 if baseline <= 0 else min(100.0, (raw_score / baseline) * 100.0)
        tier_scores = score_tier_set(sample["present_markers"], tier_set, config)
        strategy2_raw_score = float(tier_scores["strategy2"])
        strategy2_norm_score, strategy2_baseline_mean, strategy2_baseline_n_refs = normalize_strategy2_score(
            strategy2_raw_score,
            group_def,
            strategy2_baselines,
        )
        training_rows.append(
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

    training_df = pd.DataFrame(training_rows)
    model_rows: List[Dict[str, object]] = []
    for order_name, order_df in training_df.groupby("order_name"):
        model_order_df = downsample_training_order_df(order_df)
        if len(model_order_df) < 20:
            continue
        feature_cols = [
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
        X = model_order_df[feature_cols]
        y = model_order_df["true_completeness"]
        candidate_models = [
            (
                "random_forest",
                lambda: RandomForestRegressor(
                    n_estimators=400,
                    random_state=RANDOM_SEED,
                    min_samples_leaf=2,
                ),
            ),
            (
                "extra_trees",
                lambda: ExtraTreesRegressor(
                    n_estimators=500,
                    random_state=RANDOM_SEED,
                    min_samples_leaf=2,
                ),
            ),
            (
                "hist_gbm",
                lambda: HistGradientBoostingRegressor(
                    random_state=RANDOM_SEED,
                    max_depth=8,
                    min_samples_leaf=10,
                ),
            ),
        ]

        best_model = None
        best_model_name = ""
        best_mae = None
        best_r2 = None
        best_validation_type = ""
        best_validation_families = 0

        family_names = [family for family in sorted(model_order_df["family_name"].unique()) if family]
        use_lofo = len(family_names) >= 2

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
                    model.fit(
                        model_order_df.loc[train_mask, feature_cols],
                        model_order_df.loc[train_mask, "true_completeness"],
                    )
                    pred = model.predict(model_order_df.loc[test_mask, feature_cols])
                    preds.extend(pred.tolist())
                    truth.extend(model_order_df.loc[test_mask, "true_completeness"].tolist())
                    validation_families += 1

            if not preds:
                X_train, X_test, y_train, y_test = train_test_split(
                    X,
                    y,
                    test_size=0.2,
                    random_state=RANDOM_SEED,
                )
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

        model_rows.append(
            {
                "order_name": order_name,
                "n_training_samples": len(order_df),
                "n_training_samples_used": len(model_order_df),
                "model_name": best_model_name,
                "validation_type": best_validation_type,
                "n_validation_families": best_validation_families,
                "test_r2": round(float(best_r2), 4),
                "test_mae": round(float(best_mae), 4),
                "model": best_model,
            }
        )
    return training_df, pd.DataFrame(model_rows)


def summarize_by_level(rows: List[Dict[str, object]], level_key: str, output_path: Path) -> None:
    groups: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row[level_key])].append(row)

    summary_rows: List[Dict[str, object]] = []
    for group_name, items in groups.items():
        summary_rows.append(
            {
                level_key: group_name,
                "n_genomes": len(items),
                "mean_strategy1": round(mean(float(item["strategy1"]) for item in items), 2),
                "median_strategy1": round(median(float(item["strategy1"]) for item in items), 2),
                "mean_strategy2": round(mean(float(item["strategy2"]) for item in items), 2),
                "median_strategy2": round(median(float(item["strategy2"]) for item in items), 2),
                "mean_strategy3": round(mean(float(item["strategy3"]) for item in items), 2),
                "median_strategy3": round(median(float(item["strategy3"]) for item in items), 2),
            }
        )

    summary_rows.sort(key=lambda row: (-float(row["mean_strategy3"]), row[level_key]))
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                level_key,
                "n_genomes",
                "mean_strategy1",
                "median_strategy1",
                "mean_strategy2",
                "median_strategy2",
                "mean_strategy3",
                "median_strategy3",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)


def describe(values: List[float]) -> Tuple[float, float, int]:
    return round(mean(values), 2), round(median(values), 2), sum(v >= 90.0 for v in values)


def downsample_training_order_df(order_df: pd.DataFrame) -> pd.DataFrame:
    family_names = [family for family in sorted(order_df["family_name"].unique()) if family]
    if len(family_names) >= 2:
        parts = []
        for family_name, family_df in order_df.groupby("family_name", sort=False):
            if len(family_df) > MAX_SYNTHETIC_SAMPLES_PER_FAMILY:
                parts.append(
                    family_df.sample(
                        n=MAX_SYNTHETIC_SAMPLES_PER_FAMILY,
                        random_state=RANDOM_SEED,
                    )
                )
            else:
                parts.append(family_df)
        return pd.concat(parts, ignore_index=True)

    if len(order_df) > MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY:
        return order_df.sample(
            n=MAX_SYNTHETIC_SAMPLES_SINGLE_FAMILY,
            random_state=RANDOM_SEED,
        ).reset_index(drop=True)
    return order_df.reset_index(drop=True)


def main() -> None:
    args = parse_args()
    repo_dir = Path.cwd()
    run_dir = (repo_dir / args.run_dir).resolve()
    resources_dir = (repo_dir / args.resources_dir).resolve()

    summary_rows = load_strategy1_summary(run_dir)
    metadata = load_metadata(run_dir)
    order_markers, order_baselines = load_order_table(resources_dir)
    labels = load_labels(resources_dir)
    accessions = [row["query"] for row in summary_rows]
    query_counts = build_query_counts(run_dir, accessions)
    relevant_orders = collect_relevant_orders(summary_rows)
    order_profiles, family_profiles = parse_reference_marker_profiles(
        resources_dir,
        labels,
        order_markers,
        relevant_orders,
    )
    order_stats_cache, family_stats_cache = build_stats_cache(
        relevant_orders,
        order_markers,
        order_profiles,
        family_profiles,
    )

    suffix = args.output_suffix
    best_config_json = args.best_config_json.strip()
    if best_config_json:
        best_config = load_config_from_json((repo_dir / best_config_json).resolve())
    else:
        synthetic_tuning_samples = simulate_samples(
            labels,
            order_profiles,
            order_markers,
            SIMULATIONS_PER_LEVEL_TUNING,
        )

        tuning_results: List[Dict[str, object]] = []
        for config in build_config_grid():
            tuning_results.append(
                evaluate_config(
                    config,
                    summary_rows,
                    metadata,
                    query_counts,
                    order_markers,
                    order_profiles,
                    family_profiles,
                    order_stats_cache,
                    family_stats_cache,
                    synthetic_tuning_samples,
                )
            )

        tuning_results.sort(
            key=lambda row: (
                float(row["composite_rank_score"]),
                float(row["synthetic_mae"]),
                -float(row["isolate_mean"]),
                -int(row["isolate_ge90"]),
            )
        )

        tuning_path = run_dir / f"strategy2_tuning_results_{suffix}.tsv"
        with tuning_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(tuning_results[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(tuning_results)

        best_row = tuning_results[0]
        best_config = Strategy2Config(
            name=str(best_row["name"]),
            core_prevalence=float(best_row["core_prevalence"]),
            shared_prevalence=float(best_row["shared_prevalence"]),
            accessory_prevalence=float(best_row["accessory_prevalence"]),
            core_max_copy=float(best_row["core_max_copy"]),
            shared_max_copy=float(best_row["shared_max_copy"]),
            family_min_ref_genomes=int(best_row["family_min_ref_genomes"]),
            family_min_informative_markers=int(best_row["family_min_informative_markers"]),
            weight_core=float(best_row["weight_core"]),
            weight_shared=float(best_row["weight_shared"]),
            weight_accessory=float(best_row["weight_accessory"]),
        )

    best_config_path = run_dir / f"strategy2_best_config_{suffix}.json"
    best_config_path.write_text(json.dumps(asdict(best_config), indent=2, sort_keys=True) + "\n")
    strategy2_baselines = build_strategy2_baselines(
        best_config,
        labels,
        order_profiles,
        order_stats_cache,
        family_profiles,
        family_stats_cache,
    )

    tier_cache: Dict[Tuple[str, str], Tuple[GroupDefinition, TierSet]] = {}

    def get_group(order_name: str, family_name: str) -> Tuple[GroupDefinition, TierSet]:
        key = (order_name, family_name)
        if key not in tier_cache:
            tier_cache[key] = select_group_and_tiers(
                order_name,
                family_name,
                best_config,
                order_profiles,
                family_profiles,
                order_stats_cache,
                family_stats_cache,
            )
        return tier_cache[key]

    comparison_rows: List[Dict[str, object]] = []
    for row in summary_rows:
        accession = row["query"]
        order_token = extract_primary_token(row["order"])
        family_token = extract_primary_token(row["family"])
        fallback_score = round(float(row["order_completeness"]), 2)
        if order_token not in order_markers:
            support_metrics = compute_support_metrics(
                informative_marker_count=0,
                core_marker_count=0,
                shared_marker_count=0,
                core_fraction=0.0,
                shared_fraction=0.0,
                baseline_n_refs=0,
                group_level="unavailable",
            )
            comparison_rows.append(
                {
                    "accession_number": accession,
                    "genome_name": metadata[accession]["genome_name"],
                    "ncbi_family": metadata[accession]["ncbi_family"],
                    "ncbi_genus": metadata[accession]["ncbi_genus"],
                    "gvclass_order_reference": order_token or "NA",
                    "gvclass_family_reference": family_token or "NA",
                    "strategy2_group_level": "unavailable",
                    "strategy2_group_name": "fallback_strategy1",
                    "strategy2_core_marker_count": 0,
                    "strategy2_shared_marker_count": 0,
                    "strategy2_accessory_marker_count": 0,
                    "strategy1": fallback_score,
                    "strategy1_raw": float(row["order_completeness_raw"]),
                    "strategy2_raw": fallback_score,
                    "strategy2": fallback_score,
                    "strategy2_baseline_mean": 100.0,
                    "strategy2_baseline_n_refs": 0,
                    "strategy2_core_fraction": 0.0,
                    "strategy2_shared_fraction": 0.0,
                    "strategy2_accessory_fraction": 0.0,
                    "strategy2_informative_marker_count": support_metrics["informative_marker_count"],
                    "strategy2_informative_present_count": support_metrics["informative_present_count"],
                    "strategy2_informative_fraction": support_metrics["informative_fraction"],
                    "strategy2_support_score": support_metrics["support_score"],
                    "strategy3_ood_flag": support_metrics["ood_flag"],
                    "strategy3": fallback_score,
                }
            )
            continue

        group_def, tier_set = get_group(order_token, family_token)
        present_markers = {m for m, c in query_counts[accession].items() if c > 0}
        tier_scores = score_tier_set(present_markers, tier_set, best_config)
        strategy2_raw = float(tier_scores["strategy2"])
        strategy2_norm, strategy2_baseline_mean, strategy2_baseline_n_refs = normalize_strategy2_score(
            strategy2_raw,
            group_def,
            strategy2_baselines,
        )
        support_metrics = compute_support_metrics(
            informative_marker_count=int(tier_scores["informative_marker_count"]),
            core_marker_count=int(tier_scores["core_marker_count"]),
            shared_marker_count=int(tier_scores["shared_marker_count"]),
            core_fraction=float(tier_scores["core_fraction"]),
            shared_fraction=float(tier_scores["shared_fraction"]),
            baseline_n_refs=int(strategy2_baseline_n_refs),
            group_level=group_def.level,
        )
        comparison_rows.append(
            {
                "accession_number": accession,
                "genome_name": metadata[accession]["genome_name"],
                "ncbi_family": metadata[accession]["ncbi_family"],
                "ncbi_genus": metadata[accession]["ncbi_genus"],
                "gvclass_order_reference": order_token,
                "gvclass_family_reference": family_token,
                "strategy2_group_level": group_def.level,
                "strategy2_group_name": group_def.name,
                "strategy2_core_marker_count": tier_scores["core_marker_count"],
                "strategy2_shared_marker_count": tier_scores["shared_marker_count"],
                "strategy2_accessory_marker_count": tier_scores["accessory_marker_count"],
                "strategy1": float(row["order_completeness"]),
                "strategy1_raw": float(row["order_completeness_raw"]),
                "strategy2_raw": strategy2_raw,
                "strategy2": strategy2_norm,
                "strategy2_baseline_mean": round(strategy2_baseline_mean, 2),
                "strategy2_baseline_n_refs": strategy2_baseline_n_refs,
                "strategy2_core_fraction": tier_scores["core_fraction"],
                "strategy2_shared_fraction": tier_scores["shared_fraction"],
                "strategy2_accessory_fraction": tier_scores["accessory_fraction"],
                "strategy2_informative_marker_count": support_metrics["informative_marker_count"],
                "strategy2_informative_present_count": support_metrics["informative_present_count"],
                "strategy2_informative_fraction": support_metrics["informative_fraction"],
                "strategy2_support_score": support_metrics["support_score"],
                "strategy3_ood_flag": support_metrics["ood_flag"],
                "strategy3": 0.0,
            }
        )

    training_df, model_df = build_strategy3_training_table(
        best_config,
        labels,
        order_profiles,
        order_markers,
        order_baselines,
        order_stats_cache,
        family_profiles,
        family_stats_cache,
        strategy2_baselines,
    )
    models = {row["order_name"]: row["model"] for _, row in model_df.iterrows()}

    for row in comparison_rows:
        order_name = str(row["gvclass_order_reference"])
        model = models.get(order_name)
        if model is None:
            row["strategy3"] = round(float(row["strategy1"]), 2)
            continue
        group_is_family = 1 if row["strategy2_group_level"] == "family" else 0
        X = pd.DataFrame(
            [
                {
                    "raw_score": float(row["strategy1_raw"]),
                    "strategy1_score": float(row["strategy1"]),
                    "strategy2_raw_score": float(row["strategy2_raw"]),
                    "strategy2_score": float(row["strategy2"]),
                    "strategy2_baseline_mean": float(row["strategy2_baseline_mean"]),
                    "strategy2_baseline_n_refs": float(row["strategy2_baseline_n_refs"]),
                    "core_fraction": float(row["strategy2_core_fraction"]),
                    "shared_fraction": float(row["strategy2_shared_fraction"]),
                    "accessory_fraction": float(row["strategy2_accessory_fraction"]),
                    "n_expected_markers": len(order_markers.get(order_name, [])),
                    "n_core_markers": int(row["strategy2_core_marker_count"]),
                    "n_shared_markers": int(row["strategy2_shared_marker_count"]),
                    "n_accessory_markers": int(row["strategy2_accessory_marker_count"]),
                    "group_is_family": group_is_family,
                }
            ]
        )
        prediction = float(model.predict(X)[0])
        bounded = max(0.0, min(100.0, prediction))
        capped = apply_ood_cap(
            prediction=bounded,
            strategy1_score=float(row["strategy1"]),
            strategy2_score=float(row["strategy2"]),
            ood_flag=str(row["strategy3_ood_flag"]),
            support_score=float(row["strategy2_support_score"]),
            baseline_n_refs=int(row["strategy2_baseline_n_refs"]),
        )
        row["strategy3"] = round(capped, 2)

    comparison_rows.sort(key=lambda row: str(row["accession_number"]))

    comparison_path = run_dir / f"completeness_strategy_comparison_{suffix}.tsv"
    with comparison_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(comparison_rows)

    summarize_by_level(
        comparison_rows,
        "ncbi_family",
        run_dir / f"completeness_by_family_all_strategies_{suffix}.tsv",
    )
    summarize_by_level(
        comparison_rows,
        "ncbi_genus",
        run_dir / f"completeness_by_genus_all_strategies_{suffix}.tsv",
    )

    metrics_path = run_dir / f"completeness_strategy3_model_metrics_{suffix}.tsv"
    if not model_df.empty:
        model_df.drop(columns=["model"]).to_csv(metrics_path, sep="\t", index=False)

    s1 = describe([float(row["strategy1"]) for row in comparison_rows])
    s2 = describe([float(row["strategy2"]) for row in comparison_rows])
    s3 = describe([float(row["strategy3"]) for row in comparison_rows])

    print(f"best_config\t{best_config.name}")
    print(f"strategy1_mean\t{s1[0]}")
    print(f"strategy1_median\t{s1[1]}")
    print(f"strategy1_ge90\t{s1[2]}")
    print(f"strategy2_mean\t{s2[0]}")
    print(f"strategy2_median\t{s2[1]}")
    print(f"strategy2_ge90\t{s2[2]}")
    print(f"strategy3_mean\t{s3[0]}")
    print(f"strategy3_median\t{s3[1]}")
    print(f"strategy3_ge90\t{s3[2]}")


if __name__ == "__main__":
    main()
