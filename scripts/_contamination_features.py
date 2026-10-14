"""Shared feature recomputation for the contamination benchmark/retraining scripts.

The contamination model's training features (``CONTAMINATION_MODEL_FEATURES``) are
no longer all present in ``gvclass_summary.tsv``: the final-summary overhaul
dropped the v1 rule-based signals, and the marker-panel standardization dropped the
cellular/phage marker counts. Scraping the summary table therefore silently zeros
those features, which would degrade any retrained model.

This helper instead recomputes the full feature dict directly from a query's
per-query output directory (``hmmout/``, ``stats/``, ``blastp_out/``, ``query_*``)
using the runtime :class:`~src.core.summarize_full.FullSummarizer` — the single
source of truth that the live pipeline uses to feed the model.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

from src.core.contamination_scoring import CONTAMINATION_MODEL_FEATURES
from src.core.summarize_full import FullSummarizer

COMPLETENESS_MODE = "novelty-aware"
CONTAMINATION_FEATURE_PROFILE = "feature_v2_unvalidated"
BENCHMARK_RUN_PROFILE = {
    "completeness_mode": COMPLETENESS_MODE,
    "sensitive_mode": True,
    "mode_fast": True,
    "contamination_feature_profile": CONTAMINATION_FEATURE_PROFILE,
    "contamination_model_feature_names": CONTAMINATION_MODEL_FEATURES,
}


def load_tree_nn_results(
    tree_nn_file: Path,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Reconstruct ``tree_nn_results`` (marker -> query_protein -> {neighbor: dist})
    from a ``stats/<query>.tree_nn`` TSV written by the runtime."""
    results: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(dict)
    if not tree_nn_file.exists():
        return results
    with tree_nn_file.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "marker":
                continue
            marker, query_protein, neighbor, distance = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
            )
            try:
                dist = float(distance)
            except ValueError:
                continue
            results[marker].setdefault(query_protein, {})[neighbor] = dist
    return results


def compute_query_features(
    query_id: str,
    query_output_dir: Path,
    database_path: Path,
    summarizer: Optional[FullSummarizer] = None,
    sensitive_mode: bool = True,
    completeness_mode: str = COMPLETENESS_MODE,
    mode_fast: bool = True,
) -> Dict[str, object]:
    """Recompute the full feature dict for one query from its per-query output dir.

    Returns the complete ``summary_data`` dict, a superset of
    ``CONTAMINATION_MODEL_FEATURES`` (v1 rule-based signals, ``order_dup`` /
    ``gvog8_dup``, the cellular/phage marker counts, and the per-contig
    ``suspicious_*`` / ``cellular_coherent_*`` features). Pass a shared
    ``summarizer`` to avoid reloading the database for every query.
    """
    if summarizer is None:
        summarizer = FullSummarizer(
            database_path,
            completeness_mode=completeness_mode,
            sensitive_mode=sensitive_mode,
        )
    elif (
        summarizer.completeness_mode != completeness_mode
        or summarizer.sensitive_mode != sensitive_mode
    ):
        raise ValueError(
            "Benchmark summarizer settings do not match the requested runtime profile"
        )
    tree_nn = load_tree_nn_results(query_output_dir / "stats" / f"{query_id}.tree_nn")
    result = summarizer.summarize_query_full(
        query_id, query_output_dir, tree_nn, mode_fast=mode_fast
    )
    missing = [name for name in CONTAMINATION_MODEL_FEATURES if name not in result]
    if missing:
        raise ValueError(
            "Runtime summary is missing contamination model features: "
            + ", ".join(missing)
        )
    runtime_profile = result.get("contamination_feature_version")
    if runtime_profile != CONTAMINATION_FEATURE_PROFILE:
        raise ValueError(
            "Runtime contamination feature profile does not match the benchmark: "
            f"{runtime_profile!r}"
        )
    result["contamination_feature_profile"] = runtime_profile
    result["benchmark_completeness_mode"] = completeness_mode
    result["benchmark_sensitive_mode"] = sensitive_mode
    result["benchmark_mode_fast"] = mode_fast
    result["contamination_model_feature_names_csv"] = ",".join(
        CONTAMINATION_MODEL_FEATURES
    )
    return result


def contamination_feature_vector(summary: Dict[str, object]) -> Dict[str, object]:
    """Return the runtime model vector in its declared order."""
    missing = [name for name in CONTAMINATION_MODEL_FEATURES if name not in summary]
    if missing:
        raise ValueError(
            "Runtime summary is missing contamination model features: "
            + ", ".join(missing)
        )
    return {name: summary[name] for name in CONTAMINATION_MODEL_FEATURES}


def validate_source_split_contract(feature_df, source_column: str = "base_accession"):
    """Return train and final-test rows after a disjoint source-level check."""
    required = {"split", source_column}
    missing = sorted(required.difference(feature_df.columns))
    if missing:
        raise ValueError("Missing source split columns: " + ", ".join(missing))
    split_names = set(feature_df["split"])
    if split_names != {"train", "test"}:
        raise ValueError("Source split must contain only train and test rows")
    train_df = feature_df[feature_df["split"] == "train"].copy()
    test_df = feature_df[feature_df["split"] == "test"].copy()
    if train_df.empty or test_df.empty:
        raise ValueError(
            "Both train and independent final-test source splits are required"
        )
    overlap = set(train_df[source_column]).intersection(test_df[source_column])
    if overlap:
        raise ValueError(
            "Source genomes occur in both train and final-test splits: "
            + ", ".join(sorted(str(value) for value in overlap))
        )
    return train_df, test_df
