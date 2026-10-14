"""Fit the frozen GVClass experimental quality candidates on development data."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import platform
import warnings
from collections.abc import Mapping, Sequence
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.config.marker_sets import GROUP_TO_MODELS, MODEL_TO_GROUP
from src.core.quality_features import (
    COMPLETENESS_MODEL_FEATURES,
    COMPLETENESS_RAW_INPUTS,
    COMPOSITION_FEATURES,
    CONTAMINATION_FEATURES,
    build_completeness_features_from_counts,
)
from src.core.quality_models import (
    COMPLETENESS_BUNDLE_SCHEMA,
    CONTAMINATION_BUNDLE_SCHEMA,
    EXPERIMENTAL_COMPLETENESS_RELIABILITY,
    EXPERIMENTAL_COMPLETENESS_STRATEGY,
    EXPERIMENTAL_CONTAMINATION_RELIABILITY,
    EXPERIMENTAL_CONTAMINATION_STRATEGY,
    FAMILY_CATALOG_SHA256,
)
from src.utils.common import sha256_file

PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_DIR = PROJECT_ROOT / "tasks/quality-improvement-20260910"
COLLECTED_PATH = EXPERIMENT_DIR / "panel_collected_development.tsv"
COLLECTED_METADATA_PATH = EXPERIMENT_DIR / "panel_collected_development.json"
FAMILY_COUNTS_PATH = EXPERIMENT_DIR / "panel_family_counts.tsv"
FAMILY_METADATA_PATH = EXPERIMENT_DIR / "panel_family_counts.json"
COMPLETENESS_SUMMARY_PATH = EXPERIMENT_DIR / "panel_completeness_families_summary.json"
COMPOSITION_FEATURE_PATH = EXPERIMENT_DIR / "panel_composition_window_features.parquet"
COMPOSITION_QUERY_PATH = EXPERIMENT_DIR / "panel_composition_query_windows.tsv"
COMPOSITION_PARENT_PATH = EXPERIMENT_DIR / "panel_composition_parent_windows.tsv"
COMPOSITION_CACHE_PATH = EXPERIMENT_DIR / "panel_composition_cache.json"
COMPOSITION_OOF_PATH = EXPERIMENT_DIR / "panel_composition_oof_predictions.tsv"
COMPOSITION_RUN_PATH = EXPERIMENT_DIR / "panel_composition_run.json"
MARKER_GROUPS_PATH = PROJECT_ROOT / "src/config/marker_groups.tsv"
COMBINED_HMM_PATH = PROJECT_ROOT / "resources/hmm/combined.hmm"
CURRENT_REFERENCE_PATHS = (
    PROJECT_ROOT / "resources/parquet/labels/labels.parquet",
    PROJECT_ROOT / "resources/parquet/faa.parquet",
)
EXPECTED_INPUT_SHA256 = {
    "panel_collected_development.tsv": (
        "b5475e611610e712067e173a5d8ec50da0b20611011782c0d22c9dadf4b3594f"
    ),
    "panel_collected_development.json": (
        "991fa1fd614cb5e7c750f74ab76b60c3fa569512fe8f980969f1c551da400515"
    ),
    "panel_family_counts.tsv": (
        "7f2d7483e394f56b92168523d4c6ce591999b4e2458d1b8f666fb8f041c1bfcb"
    ),
    "panel_family_counts.json": (
        "cfd287f178d99b37fa6b5b1825f9faa589222a8ea732f4de1bc28ea40684306d"
    ),
    "panel_completeness_families_summary.json": (
        "65a4ef45a5c171cdcf2a23e31275d09bde61b948036785e13475ec7ca8bf588c"
    ),
    "panel_composition_window_features.parquet": (
        "b6988b96b08051ec88e54c3f03bb6a8822f0f55ff6230d43693f6d80f38c2986"
    ),
    "panel_composition_query_windows.tsv": (
        "16db5a662dd648bbb4b268fcc1cc833b65b26d4363fbd8e89b01bac6577f007c"
    ),
    "panel_composition_parent_windows.tsv": (
        "3396fd54d908f33aba9eca33284c4cce023d569c19b6ac129ef7b6a82b50b6bd"
    ),
    "panel_composition_cache.json": (
        "caf18c0df00d957e499380bc949d0e4ea9daa869364378becd3f5a73bf99adf3"
    ),
    "panel_composition_oof_predictions.tsv": (
        "edfc6e2ad3346eb7658692d768de8bc7cddb03149eaf58c266ee12b138950223"
    ),
    "panel_composition_run.json": (
        "2695f59681724cdc5efffe6946dce4437e5ce5934e740a4a61510839dce7deab"
    ),
    "marker_groups.tsv": (
        "60f893bb5951c4e3911df91f892849ffaa81649a385143b48cf5641af8b3b1b5"
    ),
    "combined.hmm": (
        "78d57a14c66edd0f7cb39c5a33ca88158329aafcae6ad98c5976c5fd39cac2b1"
    ),
}
EXPECTED_ROWS = 867
EXPECTED_FAMILY_VOCABULARY = 494
SEED = 20260911
LOGGER = logging.getLogger(__name__)


def build_models(output_dir: Path, evidence_dir: Path) -> dict[str, object]:
    """Fit both candidates and write fresh bundles plus review evidence."""
    paths = _output_paths(output_dir, evidence_dir)
    _refuse_overwrite(paths.values())
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    input_hashes = _verify_inputs()
    current_reference = _current_reference_snapshot()
    collected = _load_collected()

    completeness_bundle, completeness_record = _fit_completeness(
        collected,
        input_hashes,
        current_reference,
    )
    contamination_bundle, contamination_record, bridge = _fit_contamination(
        collected,
        input_hashes,
        current_reference,
    )
    _dump_joblib(completeness_bundle, paths["completeness_model"])
    _dump_joblib(contamination_bundle, paths["contamination_model"])
    _write_tsv(bridge, paths["bridge"])
    model_records = {
        "completeness": _finish_model_record(
            completeness_record,
            paths["completeness_model"],
        ),
        "contamination": _finish_model_record(
            contamination_record,
            paths["contamination_model"],
        ),
    }
    _write_json(model_records["completeness"], paths["completeness_card"])
    _write_json(model_records["contamination"], paths["contamination_card"])
    manifest = _build_manifest(
        paths,
        input_hashes,
        current_reference,
        model_records,
    )
    _write_json(manifest, paths["manifest"])
    return manifest


def _output_paths(output_dir: Path, evidence_dir: Path) -> dict[str, Path]:
    return {
        "completeness_model": output_dir / "completeness.joblib",
        "contamination_model": output_dir / "contamination.joblib",
        "completeness_card": evidence_dir / "completeness-model-card.json",
        "contamination_card": evidence_dir / "contamination-model-card.json",
        "bridge": evidence_dir / "composition-fullfit-bridge.tsv",
        "manifest": evidence_dir / "build-manifest.json",
    }


def _refuse_overwrite(paths: Sequence[Path]) -> None:
    existing = [path for path in paths if path.exists()]
    if existing:
        names = ", ".join(map(str, existing))
        raise FileExistsError(f"Refusing to overwrite model outputs: {names}")


def _verify_inputs() -> dict[str, str]:
    paths = (
        COLLECTED_PATH,
        COLLECTED_METADATA_PATH,
        FAMILY_COUNTS_PATH,
        FAMILY_METADATA_PATH,
        COMPLETENESS_SUMMARY_PATH,
        COMPOSITION_FEATURE_PATH,
        COMPOSITION_QUERY_PATH,
        COMPOSITION_PARENT_PATH,
        COMPOSITION_CACHE_PATH,
        COMPOSITION_OOF_PATH,
        COMPOSITION_RUN_PATH,
        MARKER_GROUPS_PATH,
        COMBINED_HMM_PATH,
    )
    observed = {path.name: sha256_file(path) for path in paths}
    if observed != EXPECTED_INPUT_SHA256:
        changed = sorted(
            name
            for name, expected in EXPECTED_INPUT_SHA256.items()
            if observed.get(name) != expected
        )
        raise ValueError(f"Frozen model inputs changed: {', '.join(changed)}")
    return observed


def _current_reference_snapshot() -> dict[str, object]:
    members = {
        str(path.relative_to(PROJECT_ROOT)): {
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        }
        for path in CURRENT_REFERENCE_PATHS
    }
    return {
        "interpretation": (
            "current runtime reference bytes; development profiles did not hash "
            "these members, so this does not prove historical byte identity"
        ),
        "members": members,
    }


def _load_collected() -> pd.DataFrame:
    metadata = json.loads(COLLECTED_METADATA_PATH.read_text(encoding="utf-8"))
    if (
        metadata.get("status") != "complete_development_manifest"
        or metadata.get("sample_count") != EXPECTED_ROWS
        or metadata.get("table_sha256") != EXPECTED_INPUT_SHA256[COLLECTED_PATH.name]
    ):
        raise ValueError(
            "Collected development metadata differs from the frozen contract"
        )
    frame = pd.read_csv(COLLECTED_PATH, sep="\t", low_memory=False)
    if len(frame) != EXPECTED_ROWS or frame["sample_id"].duplicated().any():
        raise ValueError("Collected development table does not contain 867 unique rows")
    if set(frame["split"].astype(str)) != {"development"}:
        raise ValueError("Collected model input contains a non-development row")
    return frame


def _fit_completeness(
    collected: pd.DataFrame,
    input_hashes: Mapping[str, str],
    current_reference: Mapping[str, object],
) -> tuple[dict[str, object], dict[str, object]]:
    family_table, family_metadata = _load_family_table(collected)
    observations = family_table["family_copy_counts_json"].map(_parse_counts)
    vocabulary = tuple(sorted(set().union(*(set(item) for item in observations))))
    if len(vocabulary) != EXPECTED_FAMILY_VOCABULARY:
        raise ValueError("Full-development family vocabulary does not contain 494 keys")
    summaries = [_runtime_summary(row) for row in collected.itertuples(index=False)]
    feature_rows = [
        build_completeness_features_from_counts(
            summary,
            counts,
            family_vocabulary=vocabulary,
        )
        for summary, counts in zip(summaries, observations, strict=True)
    ]
    features = pd.concat(feature_rows, ignore_index=True)
    expected_features = (
        *COMPLETENESS_MODEL_FEATURES,
        *(f"family_present::{name}" for name in vocabulary),
        *(f"family_count::{name}" for name in vocabulary),
    )
    if tuple(features.columns) != expected_features or features.shape != (867, 1_011):
        raise ValueError(
            "Completeness training feature matrix differs from the contract"
        )
    truth = pd.to_numeric(collected["true_source_retention_pct"], errors="raise")
    weights = _balanced_level_weights(collected, truth)
    pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("regressor", HistGradientBoostingRegressor(random_state=SEED)),
        ]
    )
    pipeline.fit(features, truth, regressor__sample_weight=weights)
    model_to_family, review_models = _family_catalog(family_metadata)
    metadata = _completeness_metadata(
        collected,
        weights,
        input_hashes,
        current_reference,
    )
    bundle: dict[str, object] = {
        "bundle_schema": COMPLETENESS_BUNDLE_SCHEMA,
        "family_catalog_sha256": FAMILY_CATALOG_SHA256,
        "family_vocabulary": vocabulary,
        "feature_names": expected_features,
        "libraries": _library_versions(),
        "metadata": metadata,
        "model_to_family": model_to_family,
        "pipeline": pipeline,
        "raw_input_names": COMPLETENESS_RAW_INPUTS,
        "review_models": review_models,
    }
    record = {
        "bundle_schema": COMPLETENESS_BUNDLE_SCHEMA,
        "metadata": metadata,
        "feature_count": features.shape[1],
        "family_vocabulary_size": len(vocabulary),
        "family_vocabulary_sha256": _line_digest(vocabulary),
        "training_prediction_summary": _prediction_summary(
            np.clip(pipeline.predict(features), 0.0, 100.0)
        ),
    }
    return bundle, record


def _load_family_table(
    collected: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    metadata = json.loads(FAMILY_METADATA_PATH.read_text(encoding="utf-8"))
    if (
        metadata.get("status") != "complete_full_development"
        or metadata.get("sample_count") != EXPECTED_ROWS
        or metadata.get("catalog_sha256") != FAMILY_CATALOG_SHA256
    ):
        raise ValueError("Family-count metadata differs from the frozen contract")
    table = pd.read_csv(FAMILY_COUNTS_PATH, sep="\t", low_memory=False)
    merged = collected[["sample_id", "archive_sha256"]].merge(
        table,
        on="sample_id",
        how="left",
        validate="one_to_one",
        suffixes=("_collected", "_family"),
    )
    if (
        len(merged) != EXPECTED_ROWS
        or merged["hmm_artifact_status"].ne("available").any()
    ):
        raise ValueError("Family counts omit a development query")
    if not merged["archive_sha256_collected"].equals(merged["archive_sha256_family"]):
        raise ValueError("Family counts do not match the collected query archives")
    return merged, metadata


def _parse_counts(value: object) -> dict[str, int]:
    parsed = json.loads(str(value))
    if not isinstance(parsed, dict):
        raise TypeError("Family counts are not a JSON mapping")
    counts = {str(name): int(count) for name, count in parsed.items()}
    if any(count <= 0 for count in counts.values()):
        raise ValueError("Family counts contain a nonpositive value")
    return counts


def _runtime_summary(row: object) -> dict[str, object]:
    return {
        name.removeprefix("runtime_"): getattr(row, name)
        for name in COMPLETENESS_RAW_INPUTS
    }


def _family_catalog(
    family_metadata: Mapping[str, object],
) -> tuple[dict[str, str], tuple[str, ...]]:
    review_models = tuple(map(str, family_metadata["review_excluded_models"]))
    hmm_models = _hmm_models()
    model_to_family = {
        model: (
            f"group::{MODEL_TO_GROUP[model]}"
            if model in MODEL_TO_GROUP
            else f"model::{model}"
        )
        for model in hmm_models
        if model not in review_models
    }
    singletons = sorted(set(hmm_models) - set(MODEL_TO_GROUP) - set(review_models))
    catalog_body = {
        "groups": {
            key: sorted(value) for key, value in sorted(GROUP_TO_MODELS.items())
        },
        "review_models": sorted(review_models),
        "singletons": singletons,
    }
    digest = hashlib.sha256(_canonical_json(catalog_body).encode()).hexdigest()
    if digest != FAMILY_CATALOG_SHA256:
        raise ValueError(
            "Current production family catalogue differs from the frozen fit"
        )
    return model_to_family, review_models


def _hmm_models() -> tuple[str, ...]:
    with COMBINED_HMM_PATH.open(encoding="utf-8") as handle:
        models = tuple(
            line.split(maxsplit=1)[1].strip()
            for line in handle
            if line.startswith("NAME  ")
        )
    if len(models) != len(set(models)):
        raise ValueError("Combined HMM resource contains duplicate model names")
    return models


def _completeness_metadata(
    collected: pd.DataFrame,
    weights: np.ndarray,
    input_hashes: Mapping[str, str],
    current_reference: Mapping[str, object],
) -> dict[str, object]:
    support_summary = _completeness_support_summary(collected)
    return {
        "applicability": "nucleotide_extended_sensitive_only",
        "avgdist_meaning": "mean accepted nearest-neighbor tree distance",
        "estimator": "median_imputer_then_HistGradientBoostingRegressor_defaults",
        "frozen_input_sha256": dict(input_hashes),
        "current_runtime_reference_snapshot": dict(current_reference),
        "model_name": EXPERIMENTAL_COMPLETENESS_STRATEGY,
        "random_seed": SEED,
        "query_support_policy": (
            "accepted_order_and_at_least_one_filtered_marker_protein"
        ),
        "development_support_summary": support_summary,
        "refit_scope": "all_867_frozen_development_rows",
        "reliability": EXPERIMENTAL_COMPLETENESS_RELIABILITY,
        "source_count": int(collected["source_genome_id"].nunique()),
        "target": "true_source_retention_pct",
        "target_group_count": int(collected["target_group_id"].nunique()),
        "training_row_count": len(collected),
        "weight_sha256": hashlib.sha256(weights.astype("<f8").tobytes()).hexdigest(),
        "weight_mean": float(weights.mean()),
    }


def _completeness_support_summary(
    collected: pd.DataFrame,
) -> dict[str, int]:
    marker_totals = collected["hmm_marker_counts_json"].map(
        lambda value: sum(json.loads(value).values())
    )
    has_markers = marker_totals.gt(0)
    has_order = collected["accepted_order"].fillna("").astype(str).str.strip().ne("")
    installed_numeric = pd.to_numeric(
        collected["estimated_completeness"],
        errors="coerce",
    ).notna()
    runtime_supported = has_markers & has_order
    return {
        "frozen_numeric_prediction_count": len(collected),
        "runtime_supported_count": int(runtime_supported.sum()),
        "runtime_insufficient_marker_count": int((~has_markers).sum()),
        "runtime_unassigned_lineage_count": int((has_markers & ~has_order).sum()),
        "saved_installed_paired_count": int(installed_numeric.sum()),
        "saved_installed_paired_runtime_supported_count": int(
            (installed_numeric & runtime_supported).sum()
        ),
        "saved_installed_paired_removed_count": int(
            (installed_numeric & ~runtime_supported).sum()
        ),
    }


def _fit_contamination(
    collected: pd.DataFrame,
    input_hashes: Mapping[str, str],
    current_reference: Mapping[str, object],
) -> tuple[dict[str, object], dict[str, object], pd.DataFrame]:
    features, queries, parents = _load_composition_cache()
    window_scaler, classifier, window_record = _fit_window_classifier(features, parents)
    full_scores = _score_all_queries(features, queries, window_scaler, classifier)
    oof_scores = _load_oof_scores(collected)
    bridge = _composition_bridge(collected, oof_scores, full_scores)
    support_summary = _contamination_support_summary(
        collected,
        oof_scores,
    )
    regression_features = collected[list(CONTAMINATION_FEATURES[:-1])].apply(
        pd.to_numeric,
        errors="coerce",
    )
    regression_features["composition_cellular_score_mean"] = oof_scores[
        "composition_cellular_score_mean"
    ].to_numpy()
    if tuple(regression_features.columns) != CONTAMINATION_FEATURES:
        raise ValueError("Contamination regression feature order differs")
    truth = pd.to_numeric(collected["true_contamination_pct"], errors="raise")
    weights = _balanced_contamination_weights(collected, truth)
    regressor_pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
            (
                "regressor",
                ExtraTreesRegressor(
                    n_estimators=300,
                    random_state=SEED,
                    n_jobs=1,
                ),
            ),
        ]
    )
    regressor_pipeline.fit(
        regression_features,
        truth,
        regressor__sample_weight=weights,
    )
    bridge_summary = _bridge_summary(bridge)
    fit_evidence = {
        "development_support_summary": support_summary,
        "stack_bridge": bridge_summary,
        "window_fit": dict(window_record),
    }
    metadata = _contamination_metadata(
        collected,
        weights,
        input_hashes,
        current_reference,
        fit_evidence,
    )
    bundle: dict[str, object] = {
        "bundle_schema": CONTAMINATION_BUNDLE_SCHEMA,
        "composition_feature_names": COMPOSITION_FEATURES,
        "feature_names": CONTAMINATION_FEATURES,
        "libraries": _library_versions(),
        "metadata": metadata,
        "regressor_pipeline": regressor_pipeline,
        "window_classifier": classifier,
        "window_scaler": window_scaler,
    }
    record = {
        "bundle_schema": CONTAMINATION_BUNDLE_SCHEMA,
        "metadata": metadata,
        "composition_bridge": bridge_summary,
        "training_prediction_summary": _prediction_summary(
            np.clip(regressor_pipeline.predict(regression_features), 0.0, 100.0)
        ),
    }
    return bundle, record, bridge


def _load_composition_cache() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metadata = json.loads(COMPOSITION_CACHE_PATH.read_text(encoding="utf-8"))
    if (
        metadata.get("status") != "complete_full_development_composition_cache"
        or tuple(metadata.get("feature_columns", ())) != COMPOSITION_FEATURES
    ):
        raise ValueError("Composition cache metadata differs from the frozen contract")
    features = pd.read_parquet(COMPOSITION_FEATURE_PATH).set_index("feature_id")
    queries = pd.read_csv(COMPOSITION_QUERY_PATH, sep="\t", low_memory=False)
    parents = pd.read_csv(COMPOSITION_PARENT_PATH, sep="\t", low_memory=False)
    expected_rows = metadata["artifacts"]
    for path, frame in (
        (COMPOSITION_FEATURE_PATH, features),
        (COMPOSITION_QUERY_PATH, queries),
        (COMPOSITION_PARENT_PATH, parents),
    ):
        if len(frame) != int(expected_rows[path.name]["row_count"]):
            raise ValueError(f"Composition cache row count differs for {path.name}")
    if tuple(features.columns) != COMPOSITION_FEATURES:
        raise ValueError("Composition cache numerical feature order differs")
    return features, queries, parents


def _fit_window_classifier(
    features: pd.DataFrame,
    parents: pd.DataFrame,
) -> tuple[StandardScaler, LogisticRegression, dict[str, object]]:
    eligible = parents.loc[_booleans(parents["feature_eligible"])].copy()
    labels = ~eligible["parent_class"].eq("viral")
    if set(labels) != {False, True}:
        raise ValueError("Composition parent cache does not contain both classes")
    weights, weight_record = _balanced_window_weights(eligible, labels)
    matrix = features.loc[
        eligible["feature_id"].astype(str), list(COMPOSITION_FEATURES)
    ]
    matrix.index = eligible.index
    scaler = StandardScaler().fit(matrix, sample_weight=weights)
    classifier = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1_000)
    with warnings.catch_warnings():
        warnings.simplefilter("error", ConvergenceWarning)
        classifier.fit(
            scaler.transform(matrix),
            labels.astype(int),
            sample_weight=weights,
        )
    iterations = int(classifier.n_iter_[0])
    if iterations >= classifier.max_iter:
        raise RuntimeError("Full-development composition classifier did not converge")
    record = {
        "cellular_parent_group_count": int(
            eligible.loc[labels, "parent_group_id"].nunique()
        ),
        "cellular_window_count": int(labels.sum()),
        "convergence_iterations": iterations,
        "exact_hash_training_rows_removed": 0,
        "fit_rule": "all eligible parent occurrences; no query-dependent hash removal",
        "training_parent_group_count": int(eligible["parent_group_id"].nunique()),
        "training_window_count": len(eligible),
        "viral_parent_group_count": int(
            eligible.loc[~labels, "parent_group_id"].nunique()
        ),
        "viral_window_count": int((~labels).sum()),
        "weight": weight_record,
    }
    return scaler, classifier, record


def _balanced_window_weights(
    rows: pd.DataFrame,
    cellular: pd.Series,
) -> tuple[np.ndarray, dict[str, object]]:
    classes = pd.Series(np.where(cellular, "cellular", "viral"), index=rows.index)
    groups = rows["parent_group_id"].astype(str)
    weights = pd.Series(0.0, index=rows.index)
    for class_name in ("cellular", "viral"):
        class_rows = classes.eq(class_name)
        class_groups = sorted(set(groups[class_rows]))
        for group in class_groups:
            selected = class_rows & groups.eq(group)
            weights.loc[selected] = 0.5 / len(class_groups) / int(selected.sum())
    weights *= len(weights) / float(weights.sum())
    values = weights.to_numpy()
    if not np.isclose(values.mean(), 1.0):
        raise AssertionError("Composition window weights do not have mean one")
    return values, {
        "class_totals": {
            name: float(weights[classes.eq(name)].sum())
            for name in ("cellular", "viral")
        },
        "mean": float(values.mean()),
        "sha256": hashlib.sha256(values.astype("<f8").tobytes()).hexdigest(),
    }


def _score_all_queries(
    features: pd.DataFrame,
    queries: pd.DataFrame,
    scaler: StandardScaler,
    classifier: LogisticRegression,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for sample_id, sample in queries.groupby("sample_id", sort=True):
        support = _query_support(sample)
        eligible = sample.loc[_booleans(sample["feature_eligible"])]
        status = "available"
        score = float("nan")
        if eligible.empty:
            status = "no_supported_windows"
        elif support["composition_bp_coverage_pct"] < 80.0:
            status = "coverage_below_80_pct"
        else:
            matrix = features.loc[
                eligible["feature_id"].astype(str),
                list(COMPOSITION_FEATURES),
            ]
            probabilities = classifier.predict_proba(scaler.transform(matrix))[:, 1]
            score = float(np.average(probabilities, weights=eligible["segment_bp"]))
        rows.append(
            {
                "sample_id": str(sample_id),
                "fullfit_composition_score": score,
                "fullfit_composition_status": status,
                **support,
            }
        )
    output = pd.DataFrame(rows)
    if len(output) != EXPECTED_ROWS or output["sample_id"].duplicated().any():
        raise ValueError("Full-fit composition scoring omitted development queries")
    return output


def _query_support(sample: pd.DataFrame) -> dict[str, float | int]:
    fields = (
        "query_total_bp",
        "query_scored_bp",
        "query_tail_bp",
        "query_unscored_full_window_bp",
        "query_scored_window_count",
        "query_full_window_count",
    )
    values: dict[str, float | int] = {}
    for field in fields:
        numeric = pd.to_numeric(sample[field], errors="raise")
        if not numeric.eq(numeric.iloc[0]).all():
            raise ValueError(f"Query composition cache has inconsistent {field}")
        values[field] = int(numeric.iloc[0])
    total_bp = int(values["query_total_bp"])
    values["composition_bp_coverage_pct"] = (
        100.0 * int(values["query_scored_bp"]) / total_bp if total_bp else 0.0
    )
    return values


def _load_oof_scores(collected: pd.DataFrame) -> pd.DataFrame:
    scores = pd.read_csv(COMPOSITION_OOF_PATH, sep="\t", low_memory=False)
    required = {
        "sample_id",
        "composition_cellular_score_mean",
        "composition_status",
        "installed_contamination",
        "panel_contamination_ratio_regional_composition_full",
    }
    if not required <= set(scores):
        raise ValueError("Saved held-fold composition scores lack required fields")
    ordered = collected[["sample_id"]].merge(
        scores[list(required)],
        on="sample_id",
        how="left",
        validate="one_to_one",
    )
    if len(ordered) != EXPECTED_ROWS or ordered["composition_status"].isna().any():
        raise ValueError("Saved held-fold composition scores omit development queries")
    return ordered


def _composition_bridge(
    collected: pd.DataFrame,
    oof_scores: pd.DataFrame,
    full_scores: pd.DataFrame,
) -> pd.DataFrame:
    context = collected[
        [
            "sample_id",
            "source_genome_id",
            "development_cv_fold",
            "donor_development_cv_fold",
        ]
    ].copy()
    output = context.merge(oof_scores, on="sample_id", validate="one_to_one")
    output = output.merge(full_scores, on="sample_id", validate="one_to_one")
    output = output.rename(
        columns={
            "composition_cellular_score_mean": "held_fold_composition_score",
            "composition_status": "held_fold_composition_status",
        }
    )
    output["fullfit_minus_held_fold_score"] = (
        output["fullfit_composition_score"] - output["held_fold_composition_score"]
    )
    return output


def _bridge_summary(bridge: pd.DataFrame) -> dict[str, object]:
    available = bridge[
        bridge["held_fold_composition_status"].eq("available")
        & bridge["fullfit_composition_status"].eq("available")
    ]
    differences = available["fullfit_minus_held_fold_score"].to_numpy(float)
    return {
        "comparison_is_label_free": True,
        "held_fold_status_counts": bridge["held_fold_composition_status"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "fullfit_status_counts": bridge["fullfit_composition_status"]
        .value_counts()
        .sort_index()
        .to_dict(),
        "paired_numeric_count": len(available),
        "score_difference": _distribution(differences),
        "pearson_correlation": float(
            available["fullfit_composition_score"].corr(
                available["held_fold_composition_score"]
            )
        ),
        "spearman_correlation": float(
            available["fullfit_composition_score"].corr(
                available["held_fold_composition_score"],
                method="spearman",
            )
        ),
        "by_development_fold": _grouped_bridge(bridge, "development_cv_fold"),
        "by_donor_fold": _grouped_bridge(bridge, "donor_development_cv_fold"),
        "by_source": _grouped_bridge(bridge, "source_genome_id"),
    }


def _grouped_bridge(bridge: pd.DataFrame, column: str) -> dict[str, object]:
    output: dict[str, object] = {}
    grouped = bridge.assign(**{column: bridge[column].fillna("")}).groupby(
        column,
        dropna=False,
        sort=True,
    )
    for name, group in grouped:
        values = group["fullfit_minus_held_fold_score"].dropna().to_numpy(float)
        output[str(name) or "none"] = {
            "row_count": len(group),
            "score_difference": _distribution(values),
        }
    return output


def _distribution(values: np.ndarray) -> dict[str, float | int | None]:
    if values.size == 0:
        return {"count": 0, "mean": None, "median": None, "q05": None, "q95": None}
    return {
        "count": int(values.size),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "q05": float(np.quantile(values, 0.05)),
        "q95": float(np.quantile(values, 0.95)),
    }


def _contamination_metadata(
    collected: pd.DataFrame,
    weights: np.ndarray,
    input_hashes: Mapping[str, str],
    current_reference: Mapping[str, object],
    fit_evidence: Mapping[str, object],
) -> dict[str, object]:
    return {
        "applicability": "nucleotide_extended_sensitive_only",
        "composition_support": "2kb_windows_90pct_valid_4mers_80pct_scored_bp",
        "frozen_input_sha256": dict(input_hashes),
        "current_runtime_reference_snapshot": dict(current_reference),
        "model_name": EXPERIMENTAL_CONTAMINATION_STRATEGY,
        "random_seed": SEED,
        "query_support_policy": (
            "usable_accepted_family_and_genus_with_positive_tree_vote_coverage_"
            "and_supported_composition"
        ),
        "development_support_summary": fit_evidence["development_support_summary"],
        "refit_scope": (
            "all_867_frozen_development_rows_with_saved_held_fold_stack_scores"
        ),
        "reliability": EXPERIMENTAL_CONTAMINATION_RELIABILITY,
        "stack_bridge": fit_evidence["stack_bridge"],
        "stack_construction": (
            "full-data window classifier for inference; ExtraTrees fit to saved "
            "held-outer-"
            "fold composition scores"
        ),
        "target": "true_contamination_pct",
        "target_group_count": int(collected["target_group_id"].nunique()),
        "training_row_count": len(collected),
        "weight_sha256": hashlib.sha256(weights.astype("<f8").tobytes()).hexdigest(),
        "weight_mean": float(weights.mean()),
        "window_fit": fit_evidence["window_fit"],
    }


def _contamination_support_summary(
    collected: pd.DataFrame,
    oof_scores: pd.DataFrame,
) -> dict[str, object]:
    family_supported = _usable_lineage_rank(collected["accepted_family"])
    genus_supported = _usable_lineage_rank(collected["accepted_genus"])
    family_coverage = pd.to_numeric(
        collected["spatial_viral_family_vote_coverage_bp_fraction"],
        errors="coerce",
    ).gt(0)
    genus_coverage = pd.to_numeric(
        collected["spatial_viral_genus_vote_coverage_bp_fraction"],
        errors="coerce",
    ).gt(0)
    composition_supported = oof_scores["composition_status"].eq("available")
    remaining = pd.Series(True, index=collected.index)
    summary = {"frozen_numeric_prediction_count": len(collected)}
    for name, supported in (
        ("runtime_accepted_family_missing_count", family_supported),
        ("runtime_accepted_genus_missing_count", genus_supported),
        ("runtime_family_vote_coverage_missing_count", family_coverage),
        ("runtime_genus_vote_coverage_missing_count", genus_coverage),
        ("runtime_composition_unavailable_count", composition_supported),
    ):
        summary[name] = int((remaining & ~supported).sum())
        remaining &= supported
    summary["runtime_supported_count"] = int(remaining.sum())
    candidate = pd.to_numeric(
        oof_scores["panel_contamination_ratio_regional_composition_full"],
        errors="coerce",
    )
    installed = pd.to_numeric(
        oof_scores["installed_contamination"],
        errors="coerce",
    )
    truth = pd.to_numeric(collected["true_contamination_pct"], errors="raise")
    paired = candidate.notna() & installed.notna()
    fixed_mask = paired & remaining
    summary["saved_installed_paired_count"] = int(paired.sum())
    summary["saved_installed_paired_runtime_supported_count"] = int(fixed_mask.sum())
    summary["saved_installed_paired_removed_count"] = int((paired & ~remaining).sum())
    summary["saved_paired_removed_status_counts"] = {
        name: int((paired & status).sum())
        for name, status in _contamination_status_masks(
            family_supported,
            genus_supported,
            family_coverage,
            genus_coverage,
            composition_supported,
        ).items()
        if name != "supported"
    }
    summary["posthoc_fixed_mask_paired_mae"] = _posthoc_paired_mae(
        collected,
        truth,
        candidate,
        installed,
        fixed_mask,
    )
    return summary


def _contamination_status_masks(
    family_supported: pd.Series,
    genus_supported: pd.Series,
    family_coverage: pd.Series,
    genus_coverage: pd.Series,
    composition_supported: pd.Series,
) -> dict[str, pd.Series]:
    remaining = pd.Series(True, index=family_supported.index)
    masks: dict[str, pd.Series] = {}
    for name, supported in (
        ("accepted_family_missing", family_supported),
        ("accepted_genus_missing", genus_supported),
        ("family_vote_coverage_missing", family_coverage),
        ("genus_vote_coverage_missing", genus_coverage),
        ("composition_unavailable", composition_supported),
    ):
        masks[name] = remaining & ~supported
        remaining &= supported
    masks["supported"] = remaining
    return masks


def _posthoc_paired_mae(
    collected: pd.DataFrame,
    truth: pd.Series,
    candidate: pd.Series,
    installed: pd.Series,
    fixed_mask: pd.Series,
) -> dict[str, object]:
    fixed_collected = collected.loc[fixed_mask].reset_index(drop=True)
    fixed_truth = truth.loc[fixed_mask].reset_index(drop=True)
    weights = _balanced_contamination_weights(
        fixed_collected,
        fixed_truth,
    )
    candidate_errors = np.abs(
        candidate.loc[fixed_mask].to_numpy() - fixed_truth.to_numpy()
    )
    installed_errors = np.abs(
        installed.loc[fixed_mask].to_numpy() - fixed_truth.to_numpy()
    )
    candidate_mae = float(np.average(candidate_errors, weights=weights))
    installed_mae = float(np.average(installed_errors, weights=weights))
    return {
        "interpretation": (
            "post-hoc description on the fixed runtime-support mask; not model "
            "selection or validation of the full-fit bytes"
        ),
        "paired_count": int(fixed_mask.sum()),
        "candidate_balanced_mae": candidate_mae,
        "installed_balanced_mae": installed_mae,
        "candidate_minus_installed_balanced_mae": candidate_mae - installed_mae,
    }


def _usable_lineage_rank(values: pd.Series) -> pd.Series:
    normalized = values.fillna("").astype(str).str.strip().str.lower()
    unavailable = normalized.isin(("", "nan", "none"))
    contextual = normalized.str.contains("unclassified|singleton|other", regex=True)
    return ~(unavailable | contextual)


def _balanced_contamination_weights(
    frame: pd.DataFrame,
    truth: pd.Series,
) -> np.ndarray:
    levels = np.array([0.0, 10.0, 25.0, 50.0])
    distances = np.abs(truth.to_numpy()[:, None] - levels[None, :])
    return _balanced_level_weights(
        frame, pd.Series(levels[np.argmin(distances, axis=1)])
    )


def _balanced_level_weights(frame: pd.DataFrame, levels: pd.Series) -> np.ndarray:
    work = frame[["target_group_id", "source_genome_id"]].copy()
    work["level"] = np.rint(pd.to_numeric(levels, errors="raise")).astype(int)
    group_count = work["target_group_id"].nunique()
    sources_per_group = work.groupby("target_group_id")["source_genome_id"].nunique()
    levels_per_source = work.groupby("source_genome_id")["level"].nunique()
    rows_per_level = work.groupby(["source_genome_id", "level"])["level"].size()
    values = np.asarray(
        [
            1.0
            / group_count
            / int(sources_per_group[group])
            / int(levels_per_source[source])
            / int(rows_per_level[(source, level)])
            for group, source, level in work.itertuples(index=False, name=None)
        ],
        dtype=float,
    )
    return values * len(values) / values.sum()


def _booleans(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.astype(bool)
    normalized = values.astype(str).str.strip().str.lower()
    if not set(normalized) <= {"true", "false"}:
        raise ValueError("Composition eligibility contains a non-boolean value")
    return normalized.eq("true")


def _finish_model_record(
    record: dict[str, object],
    path: Path,
) -> dict[str, object]:
    loaded = joblib.load(path)
    if not isinstance(loaded, dict):
        raise TypeError(f"Serialized candidate is not a mapping: {path}")
    return {
        **record,
        "artifact": {
            "path": str(path),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        },
        "serialization": (
            "joblib dictionary of plain sklearn estimators and standard data"
        ),
    }


def _build_manifest(
    paths: Mapping[str, Path],
    input_hashes: Mapping[str, str],
    current_reference: Mapping[str, object],
    model_records: Mapping[str, object],
) -> dict[str, object]:
    return {
        "status": "complete_full_development_candidate_refit",
        "input_sha256": dict(input_hashes),
        "current_runtime_reference_snapshot": dict(current_reference),
        "libraries": _library_versions(),
        "models": dict(model_records),
        "builder": {
            "path": str(Path(__file__)),
            "sha256": sha256_file(Path(__file__)),
        },
        "bridge_artifact": {
            "path": str(paths["bridge"]),
            "sha256": sha256_file(paths["bridge"]),
            "row_count": EXPECTED_ROWS,
        },
        "limitations": [
            "Both candidates are restricted to nucleotide, extended, sensitive runs.",
            "The serialized refits have no untouched-test error estimate.",
            "Nested development error does not directly validate the full-fit stack.",
            "The contamination candidate failed the development false-positive gate.",
            "Production queries may share exact windows with development parents.",
        ],
    }


def _prediction_summary(values: np.ndarray) -> dict[str, float | int]:
    return {
        "count": int(values.size),
        "min": float(np.min(values)),
        "median": float(np.median(values)),
        "max": float(np.max(values)),
    }


def _library_versions() -> dict[str, str]:
    return {
        "joblib": version("joblib"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "python": platform.python_version(),
        "scikit-learn": version("scikit-learn"),
    }


def _line_digest(values: Sequence[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _dump_joblib(value: object, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    joblib.dump(value, temporary, compress=3)
    temporary.replace(path)


def _write_tsv(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    temporary.replace(path)


def _write_json(value: object, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evidence-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Fit both candidates once and report the manifest path."""
    arguments = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    manifest = build_models(arguments.output_dir, arguments.evidence_dir)
    LOGGER.info("Built candidates with status %s", manifest["status"])
    print(arguments.evidence_dir / "build-manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
