"""Prototype completeness scoring strategies 2 and 3 for GVClass outputs.

This script compares three completeness strategies on an existing GVClass run:

1. Strategy 1: order-baseline normalization (already produced by strategy1 summary)
2. Strategy 2: hierarchical core-marker completeness with family/order fallback
3. Strategy 3: order-specific calibration model trained on synthetically degraded
   complete-reference marker profiles
"""

from __future__ import annotations

import argparse
import csv
import random
import tarfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split


CORE_PREVALENCE_THRESHOLD = 0.90
ACCESSORY_PREVALENCE_THRESHOLD = 0.50
FAMILY_MIN_REF_GENOMES = 5
FAMILY_MIN_CORE_MARKERS = 3
SIMULATED_COMPLETENESS_LEVELS = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]
SIMULATIONS_PER_LEVEL = 3
RANDOM_SEED = 42


@dataclass(frozen=True)
class GroupDefinition:
    level: str
    name: str
    order_name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-dir",
        default="benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10",
        help="Directory containing metadata_table.tsv and gvclass_extended_results/",
    )
    parser.add_argument(
        "--resources-dir",
        default="resources",
        help="GVClass resources directory",
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
    path = resources_dir / "labels.tsv"
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
    orders = sorted(
        {
            extract_primary_token(row["order"])
            for row in summary_rows
            if extract_primary_token(row["order"])
        }
    )
    return orders


def parse_reference_marker_profiles(
    resources_dir: Path,
    labels: Dict[str, Dict[str, str]],
    order_markers: Dict[str, List[str]],
    relevant_orders: Iterable[str],
) -> Tuple[
    Dict[str, Dict[str, Dict[str, int]]],
    Dict[str, Dict[str, Dict[str, int]]],
]:
    """Return nested copy-count matrices for order and family reference groups."""
    relevant_order_set = set(relevant_orders)
    relevant_markers = sorted(
        {
            marker
            for order_name in relevant_order_set
            for marker in order_markers.get(order_name, [])
        }
    )

    order_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    family_profiles: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(
        lambda: defaultdict(dict)
    )
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


def marker_stats_for_group(
    genome_profiles: Dict[str, Dict[str, int]],
    markers: List[str],
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


def select_group_definition(
    order_name: str,
    family_name: str,
    order_markers: Dict[str, List[str]],
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    family_profiles: Dict[str, Dict[str, Dict[str, int]]],
) -> Tuple[GroupDefinition, pd.DataFrame]:
    family_profile = family_profiles.get(family_name, {})
    family_stats = marker_stats_for_group(family_profile, order_markers[order_name])
    family_core_count = 0
    if not family_stats.empty:
        family_core_count = int(
            (
                (family_stats["prevalence"] >= CORE_PREVALENCE_THRESHOLD)
                & (family_stats["avg_copy"] <= 1.2)
            ).sum()
        )
    if len(family_profile) >= FAMILY_MIN_REF_GENOMES and family_core_count >= FAMILY_MIN_CORE_MARKERS:
        return (
            GroupDefinition(level="family", name=family_name, order_name=order_name),
            family_stats,
        )

    order_profile = order_profiles.get(order_name, {})
    order_stats = marker_stats_for_group(order_profile, order_markers[order_name])
    return (
        GroupDefinition(level="order", name=order_name, order_name=order_name),
        order_stats,
    )


def build_strategy2_sets(group_stats: pd.DataFrame) -> Tuple[List[str], List[str]]:
    if group_stats.empty:
        return [], []
    core = sorted(
        group_stats[
            (group_stats["prevalence"] >= CORE_PREVALENCE_THRESHOLD)
            & (group_stats["avg_copy"] <= 1.2)
        ]["marker"].tolist()
    )
    accessory = sorted(
        group_stats[
            (group_stats["prevalence"] >= ACCESSORY_PREVALENCE_THRESHOLD)
            & (group_stats["prevalence"] < CORE_PREVALENCE_THRESHOLD)
        ]["marker"].tolist()
    )
    return core, accessory


def compute_marker_fraction(marker_counts: Dict[str, int], markers: List[str]) -> float:
    if not markers:
        return 0.0
    present = sum(1 for marker in markers if marker_counts.get(marker, 0) > 0)
    return (present / len(markers)) * 100.0


def build_strategy3_training_table(
    order_profiles: Dict[str, Dict[str, Dict[str, int]]],
    order_markers: Dict[str, List[str]],
    strategy3_sets: Dict[str, Tuple[List[str], List[str]]],
    order_baselines: Dict[str, float],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rng = random.Random(RANDOM_SEED)
    training_rows: List[Dict[str, float | str]] = []
    model_rows: List[Dict[str, float | str]] = []

    for order_name, genome_profiles in order_profiles.items():
        markers = order_markers.get(order_name, [])
        if not markers or len(genome_profiles) < 5:
            continue
        core_markers, accessory_markers = strategy3_sets.get(order_name, ([], []))
        for genome_id, profile in genome_profiles.items():
            present_markers = [marker for marker in markers if profile.get(marker, 0) > 0]
            if not present_markers:
                continue
            for completeness_level in SIMULATED_COMPLETENESS_LEVELS:
                retain_n = max(1, round(len(present_markers) * completeness_level))
                for _ in range(SIMULATIONS_PER_LEVEL):
                    retained = set(rng.sample(present_markers, retain_n))
                    synthetic_counts = {marker: 1 for marker in retained}
                    raw_score = compute_marker_fraction(synthetic_counts, markers)
                    baseline = order_baselines.get(order_name, 0.0)
                    strategy1_score = 0.0
                    if baseline > 0:
                        strategy1_score = min(100.0, (raw_score / baseline) * 100.0)
                    strategy2_score = compute_marker_fraction(synthetic_counts, core_markers)
                    accessory_score = compute_marker_fraction(
                        synthetic_counts, accessory_markers
                    )
                    training_rows.append(
                        {
                            "order_name": order_name,
                            "raw_score": raw_score,
                            "strategy1_score": strategy1_score,
                            "strategy2_score": strategy2_score,
                            "accessory_score": accessory_score,
                            "n_expected_markers": len(markers),
                            "n_core_markers": len(core_markers),
                            "true_completeness": completeness_level * 100.0,
                        }
                    )

    if not training_rows:
        return pd.DataFrame(), pd.DataFrame()

    training_df = pd.DataFrame(training_rows)
    for order_name, order_df in training_df.groupby("order_name"):
        if len(order_df) < 20:
            continue
        X = order_df[
            [
                "raw_score",
                "strategy1_score",
                "strategy2_score",
                "accessory_score",
                "n_expected_markers",
                "n_core_markers",
            ]
        ]
        y = order_df["true_completeness"]
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=RANDOM_SEED
        )
        model = RandomForestRegressor(
            n_estimators=300,
            random_state=RANDOM_SEED,
            min_samples_leaf=2,
        )
        model.fit(X_train, y_train)
        pred = model.predict(X_test)
        model_rows.append(
            {
                "order_name": order_name,
                "n_training_samples": len(order_df),
                "test_r2": round(r2_score(y_test, pred), 4),
                "test_mae": round(mean_absolute_error(y_test, pred), 4),
                "model": model,
            }
        )
    model_df = pd.DataFrame(model_rows)
    return training_df, model_df


def summarize_by_level(
    rows: List[Dict[str, str | float]],
    level_key: str,
    output_path: Path,
) -> None:
    groups: Dict[str, List[Dict[str, str | float]]] = defaultdict(list)
    for row in rows:
        groups[str(row[level_key])].append(row)

    summary_rows: List[Dict[str, str | float]] = []
    for group_name, items in groups.items():
        summary_rows.append(
            {
                level_key: group_name,
                "n_genomes": len(items),
                "mean_strategy1": round(mean(float(item["strategy1"]) for item in items), 2),
                "median_strategy1": round(
                    median(float(item["strategy1"]) for item in items), 2
                ),
                "mean_strategy2": round(mean(float(item["strategy2"]) for item in items), 2),
                "median_strategy2": round(
                    median(float(item["strategy2"]) for item in items), 2
                ),
                "mean_strategy3": round(mean(float(item["strategy3"]) for item in items), 2),
                "median_strategy3": round(
                    median(float(item["strategy3"]) for item in items), 2
                ),
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
        resources_dir, labels, order_markers, relevant_orders
    )

    strategy3_sets: Dict[str, Tuple[List[str], List[str]]] = {}
    for order_name in relevant_orders:
        order_stats = marker_stats_for_group(
            order_profiles.get(order_name, {}),
            order_markers.get(order_name, []),
        )
        strategy3_sets[order_name] = build_strategy2_sets(order_stats)

    comparison_rows: List[Dict[str, str | float]] = []

    for row in summary_rows:
        accession = row["query"]
        order_token = extract_primary_token(row["order"])
        family_token = extract_primary_token(row["family"])
        if order_token not in order_markers:
            fallback_score = round(float(row["order_completeness"]), 2)
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
                    "strategy2_accessory_marker_count": 0,
                    "strategy1": fallback_score,
                    "strategy1_raw": float(row["order_completeness_raw"]),
                    "strategy2": fallback_score,
                    "strategy2_accessory_support": 0.0,
                    "strategy3": fallback_score,
                }
            )
            continue

        group_def, group_stats = select_group_definition(
            order_token,
            family_token,
            order_markers,
            order_profiles,
            family_profiles,
        )
        core_markers, accessory_markers = build_strategy2_sets(group_stats)

        marker_counts = query_counts[accession]
        strategy2_score = compute_marker_fraction(marker_counts, core_markers)
        accessory_score = compute_marker_fraction(marker_counts, accessory_markers)

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
                "strategy2_core_marker_count": len(core_markers),
                "strategy2_accessory_marker_count": len(accessory_markers),
                "strategy1": float(row["order_completeness"]),
                "strategy1_raw": float(row["order_completeness_raw"]),
                "strategy2": round(strategy2_score, 2),
                "strategy2_accessory_support": round(accessory_score, 2),
                "strategy3": 0.0,
            }
        )

    training_df, model_df = build_strategy3_training_table(
        order_profiles, order_markers, strategy3_sets, order_baselines
    )
    models = {row["order_name"]: row["model"] for _, row in model_df.iterrows()}

    for row in comparison_rows:
        model = models.get(str(row["gvclass_order_reference"]))
        if model is None:
            row["strategy3"] = row["strategy1"]
            continue
        X = pd.DataFrame(
            [
                {
                    "raw_score": row["strategy1_raw"],
                    "strategy1_score": row["strategy1"],
                    "strategy2_score": row["strategy2"],
                    "accessory_score": row["strategy2_accessory_support"],
                    "n_expected_markers": len(
                        order_markers[str(row["gvclass_order_reference"])]
                    ),
                    "n_core_markers": row["strategy2_core_marker_count"],
                }
            ]
        )
        prediction = float(model.predict(X)[0])
        row["strategy3"] = round(max(0.0, min(100.0, prediction)), 2)

    comparison_rows.sort(key=lambda row: str(row["accession_number"]))
    comparison_path = run_dir / "completeness_strategy_comparison.tsv"
    with comparison_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(comparison_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(comparison_rows)

    summarize_by_level(
        comparison_rows,
        "ncbi_family",
        run_dir / "completeness_by_family_all_strategies.tsv",
    )
    summarize_by_level(
        comparison_rows,
        "ncbi_genus",
        run_dir / "completeness_by_genus_all_strategies.tsv",
    )

    training_metrics_path = run_dir / "completeness_strategy3_model_metrics.tsv"
    if not model_df.empty:
        output_model_df = model_df.drop(columns=["model"]).copy()
        output_model_df.to_csv(training_metrics_path, sep="\t", index=False)

    def describe(strategy_key: str) -> Tuple[float, float, int]:
        values = [float(row[strategy_key]) for row in comparison_rows]
        return (
            round(mean(values), 2),
            round(median(values), 2),
            sum(value >= 90.0 for value in values),
        )

    s1 = describe("strategy1")
    s2 = describe("strategy2")
    s3 = describe("strategy3")
    print(f"comparison_rows\t{len(comparison_rows)}")
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
