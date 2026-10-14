"""Load and score the hash-bound experimental quality-model candidates."""

from __future__ import annotations

import hashlib
import io
import json
import platform
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import joblib
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.core.contamination_scoring import ContaminationScorer
from src.core.quality_features import (
    COMPLETENESS_MODEL_FEATURES,
    COMPLETENESS_RAW_INPUTS,
    COMPOSITION_FEATURES,
    CONTAMINATION_FEATURES,
    build_completeness_features,
    build_contamination_features,
)

__all__ = [
    "EXPERIMENTAL_COMPLETENESS_RELIABILITY",
    "EXPERIMENTAL_COMPLETENESS_STRATEGY",
    "EXPERIMENTAL_CONTAMINATION_RELIABILITY",
    "EXPERIMENTAL_CONTAMINATION_STRATEGY",
    "CandidateQualityPredictions",
    "ExperimentalQualityScorer",
    "QualityPrediction",
    "create_experimental_quality_scorer",
    "experimental_model_profile",
]

COMPLETENESS_BUNDLE_SCHEMA = "gvclass.experimental.completeness.family_union.v1"
CONTAMINATION_BUNDLE_SCHEMA = "gvclass.experimental.contamination.composition.v1"
EXPERIMENTAL_COMPLETENESS_STRATEGY = "experimental_family_union_hgb_v1"
EXPERIMENTAL_CONTAMINATION_STRATEGY = "experimental_composition_regional_et_v1"
EXPERIMENTAL_COMPLETENESS_RELIABILITY = "development_only_no_untouched_refit_test"
EXPERIMENTAL_CONTAMINATION_RELIABILITY = "development_only_failed_false_positive_gate"
COMPLETENESS_MODEL_SHA256 = (
    "a9ff2bae5e47042a7c93d17a1bf16d68cdcd72768e0541a492821a0a8e78a1e6"
)
CONTAMINATION_MODEL_SHA256 = (
    "6b5edd6921bb6115e1f8bfe18d43d408b843af88823f9644f41486c92b0218b0"
)
FAMILY_CATALOG_SHA256 = (
    "69aa643215edb66ef4e8c53fda3b91199a7f0a5c63ecac01cdb217fc59cf5eb8"
)
EXPECTED_FAMILY_VOCABULARY_SIZE = 494
EXPECTED_FAMILY_CATALOG_SIZE = 624
EXPECTED_COMPLETENESS_FEATURE_COUNT = 1_011
EXPECTED_COMPOSITION_FEATURE_COUNT = 137
EXPECTED_CONTAMINATION_FEATURE_COUNT = 14
EXPERIMENTAL_PROFILE_SCHEMA = "gvclass.experimental-quality-profile.v2"
EXPERIMENTAL_SUPPORTED_DOMAIN = "NCLDV"


@dataclass(frozen=True, slots=True)
class QualityPrediction:
    """One candidate estimate and its explicit support state."""

    value: float | None
    status: str
    detail: str


@dataclass(frozen=True, slots=True)
class CandidateQualityPredictions:
    """The completeness and contamination candidate results for one query."""

    completeness: QualityPrediction
    contamination: QualityPrediction


@dataclass(frozen=True, slots=True)
class ExperimentalQualityScorer:
    """Two validated full-development refits with one query scoring entry point."""

    completeness_bundle: dict[str, object]
    contamination_bundle: dict[str, object]

    def score_query(  # noqa: PLR0913  # FullSummarizer evidence boundary
        self,
        *,
        summary: Mapping[str, object],
        marker_hits: Mapping[str, set[str]] | None,
        query_fna: Path,
        query_output_dir: Path,
        tree_nn_results: Mapping[str, Mapping[str, Mapping[str, float]]],
        accepted_lineage: Mapping[str, str],
        mode_fast: bool,
        sensitive_mode: bool,
        contamination_scorer: ContaminationScorer,
    ) -> CandidateQualityPredictions:
        """Score one nucleotide query within the validated runtime profile."""
        applicability = _applicability(
            query_fna,
            mode_fast=mode_fast,
            sensitive_mode=sensitive_mode,
        )
        if applicability is not None:
            return CandidateQualityPredictions(applicability, applicability)
        if accepted_lineage.get("domain", "") != EXPERIMENTAL_SUPPORTED_DOMAIN:
            unsupported_domain = QualityPrediction(
                None,
                "unsupported_target_domain",
                "Experimental quality candidates require accepted domain NCLDV.",
            )
            return CandidateQualityPredictions(
                unsupported_domain,
                unsupported_domain,
            )
        completeness = self._score_completeness(
            summary,
            marker_hits,
            accepted_lineage,
        )
        contamination = self._score_contamination(
            summary,
            query_fna=query_fna,
            query_output_dir=query_output_dir,
            tree_nn_results=tree_nn_results,
            accepted_lineage=accepted_lineage,
            contamination_scorer=contamination_scorer,
        )
        return CandidateQualityPredictions(completeness, contamination)

    def _score_completeness(
        self,
        summary: Mapping[str, object],
        marker_hits: Mapping[str, set[str]] | None,
        accepted_lineage: Mapping[str, str],
    ) -> QualityPrediction:
        if marker_hits is None:
            return QualityPrediction(
                None,
                "missing_filtered_hmm_hits",
                "The exact family-union feature requires models.out.filtered.",
            )
        if not any(marker_hits.values()):
            return QualityPrediction(
                None,
                "insufficient_evidence",
                "No filtered marker protein supports a completeness estimate.",
            )
        if not str(accepted_lineage.get("order", "")).strip():
            return QualityPrediction(
                None,
                "unassigned",
                "An accepted order lineage is required for completeness support.",
            )
        features = build_completeness_features(
            summary,
            marker_hits,
            family_vocabulary=_string_tuple(
                self.completeness_bundle,
                "family_vocabulary",
            ),
            model_to_family=_string_mapping(
                self.completeness_bundle,
                "model_to_family",
            ),
            review_models=_string_tuple(self.completeness_bundle, "review_models"),
        )
        expected = _string_tuple(self.completeness_bundle, "feature_names")
        _require_columns(features, expected, "completeness")
        pipeline = _pipeline(self.completeness_bundle, "pipeline")
        value = _bounded_prediction(pipeline, features)
        return QualityPrediction(
            value,
            "available",
            (
                "Full-development refit; development cross-validation does not "
                "test these bytes."
            ),
        )

    def _score_contamination(
        self,
        summary: Mapping[str, object],
        *,
        query_fna: Path,
        query_output_dir: Path,
        tree_nn_results: Mapping[str, Mapping[str, Mapping[str, float]]],
        accepted_lineage: Mapping[str, str],
        contamination_scorer: ContaminationScorer,
    ) -> QualityPrediction:
        try:
            observations = build_contamination_features(
                summary,
                query_fna=query_fna,
                query_output_dir=query_output_dir,
                tree_nn_results=tree_nn_results,
                accepted_lineage=accepted_lineage,
                contamination_scorer=contamination_scorer,
            )
        except FileNotFoundError as error:
            return QualityPrediction(None, "missing_feature_artifact", str(error))
        if observations.lineage_status != "supported":
            return QualityPrediction(
                None,
                observations.lineage_status,
                "Accepted family/genus and viral tree-vote coverage are required.",
            )
        if observations.composition.status != "available":
            return QualityPrediction(
                None,
                observations.composition.status,
                (
                    "Composition requires at least 80% scored query bases after "
                    "the 90% valid-4-mer window filter."
                ),
            )
        scaler = _scaler(self.contamination_bundle, "window_scaler")
        classifier = _classifier(self.contamination_bundle, "window_classifier")
        transformed = scaler.transform(observations.composition.features)
        scores = classifier.predict_proba(transformed)[:, 1]
        composition_score = float(
            np.average(scores, weights=observations.composition.weights_bp)
        )
        features = observations.base_features.copy()
        features["composition_cellular_score_mean"] = composition_score
        expected = _string_tuple(self.contamination_bundle, "feature_names")
        _require_columns(features, expected, "contamination")
        pipeline = _pipeline(self.contamination_bundle, "regressor_pipeline")
        value = _bounded_prediction(pipeline, features)
        return QualityPrediction(
            value,
            "available",
            (
                "Full-development refit; nested development error does not directly "
                "validate the serialized full-fit stack."
            ),
        )


def create_experimental_quality_scorer(
    database_path: Path,
) -> ExperimentalQualityScorer:
    """Load both experimental candidates after exact byte-hash checks."""
    model_dir = database_path / "quality_models"
    completeness = _load_bundle(
        model_dir / "completeness.joblib",
        COMPLETENESS_MODEL_SHA256,
        "completeness",
    )
    contamination = _load_bundle(
        model_dir / "contamination.joblib",
        CONTAMINATION_MODEL_SHA256,
        "contamination",
    )
    _validate_completeness_bundle(completeness)
    _validate_contamination_bundle(contamination)
    return ExperimentalQualityScorer(completeness, contamination)


def experimental_model_profile(*, mode_fast: bool, sensitive_mode: bool) -> str:
    """Bind both candidate byte hashes and runtime settings into one stable ID."""
    profile = {
        "completeness_sha256": COMPLETENESS_MODEL_SHA256,
        "contamination_sha256": CONTAMINATION_MODEL_SHA256,
        "mode_fast": mode_fast,
        "schema": EXPERIMENTAL_PROFILE_SCHEMA,
        "sensitive_mode": sensitive_mode,
        "supported_domain": EXPERIMENTAL_SUPPORTED_DOMAIN,
    }
    body = json.dumps(profile, sort_keys=True, separators=(",", ":"))
    return f"experimental:{hashlib.sha256(body.encode()).hexdigest()}"


def _load_bundle(path: Path, expected_sha256: str, name: str) -> dict[str, object]:
    if len(expected_sha256) != 64:
        raise RuntimeError(f"The {name} candidate SHA-256 pin is not finalized")
    if not path.is_file():
        raise FileNotFoundError(f"Experimental {name} model is missing: {path}")
    model_bytes = path.read_bytes()
    observed = hashlib.sha256(model_bytes).hexdigest()
    if observed != expected_sha256:
        raise ValueError(
            f"Experimental {name} model SHA-256 mismatch: expected "
            f"{expected_sha256}, observed {observed}"
        )
    loaded = joblib.load(io.BytesIO(model_bytes))
    if not isinstance(loaded, dict) or any(not isinstance(key, str) for key in loaded):
        raise TypeError(
            f"Experimental {name} model bundle is not a string-keyed mapping"
        )
    return loaded


def _validate_completeness_bundle(bundle: dict[str, object]) -> None:
    required = {
        "bundle_schema",
        "family_catalog_sha256",
        "family_vocabulary",
        "feature_names",
        "libraries",
        "metadata",
        "model_to_family",
        "pipeline",
        "raw_input_names",
        "review_models",
    }
    _require_keys(bundle, required, "completeness")
    if bundle["bundle_schema"] != COMPLETENESS_BUNDLE_SCHEMA:
        raise ValueError("Experimental completeness bundle schema is unsupported")
    if bundle["family_catalog_sha256"] != FAMILY_CATALOG_SHA256:
        raise ValueError("Experimental completeness family catalogue hash differs")
    raw_inputs = _string_tuple(bundle, "raw_input_names")
    vocabulary = _string_tuple(bundle, "family_vocabulary")
    feature_names = _string_tuple(bundle, "feature_names")
    if raw_inputs != COMPLETENESS_RAW_INPUTS:
        raise ValueError("Experimental completeness raw input order differs")
    if len(vocabulary) != EXPECTED_FAMILY_VOCABULARY_SIZE:
        raise ValueError("Experimental completeness family vocabulary has wrong size")
    expected_features = (
        *COMPLETENESS_MODEL_FEATURES,
        *(f"family_present::{name}" for name in vocabulary),
        *(f"family_count::{name}" for name in vocabulary),
    )
    if (
        feature_names != expected_features
        or len(feature_names) != EXPECTED_COMPLETENESS_FEATURE_COUNT
    ):
        raise ValueError("Experimental completeness feature schema differs")
    model_to_family = _string_mapping(bundle, "model_to_family")
    _validate_metadata(
        bundle,
        expected={
            "applicability": "nucleotide_extended_sensitive_only",
            "avgdist_meaning": "mean accepted nearest-neighbor tree distance",
            "model_name": EXPERIMENTAL_COMPLETENESS_STRATEGY,
            "random_seed": 20260911,
            "query_support_policy": (
                "accepted_order_and_at_least_one_filtered_marker_protein"
            ),
            "refit_scope": "all_867_frozen_development_rows",
            "reliability": EXPERIMENTAL_COMPLETENESS_RELIABILITY,
            "target": "true_source_retention_pct",
            "training_row_count": 867,
        },
    )
    _validate_completeness_estimators(bundle, vocabulary, model_to_family)
    _validate_libraries(bundle)


def _validate_completeness_estimators(
    bundle: Mapping[str, object],
    vocabulary: tuple[str, ...],
    model_to_family: Mapping[str, str],
) -> None:
    if not set(vocabulary) <= set(model_to_family.values()):
        raise ValueError(
            "Experimental completeness vocabulary is absent from its catalogue"
        )
    if len(set(model_to_family.values())) != EXPECTED_FAMILY_CATALOG_SIZE:
        raise ValueError("Experimental completeness family catalogue has wrong size")
    pipeline = _pipeline(bundle, "pipeline")
    if tuple(pipeline.named_steps) != ("imputer", "regressor"):
        raise TypeError("Experimental completeness pipeline steps differ")
    imputer = pipeline.named_steps["imputer"]
    _validate_median_imputer(imputer, "completeness")
    regressor = pipeline.named_steps["regressor"]
    if type(regressor) is not HistGradientBoostingRegressor:
        raise TypeError("Experimental completeness regressor type differs")
    expected = HistGradientBoostingRegressor(random_state=20260911).get_params()
    if regressor.get_params() != expected:
        raise ValueError("Experimental completeness regressor parameters differ")
    if regressor.n_features_in_ != EXPECTED_COMPLETENESS_FEATURE_COUNT:
        raise ValueError("Experimental completeness fitted feature count differs")


def _validate_contamination_bundle(bundle: dict[str, object]) -> None:
    required = {
        "bundle_schema",
        "composition_feature_names",
        "feature_names",
        "libraries",
        "metadata",
        "regressor_pipeline",
        "window_classifier",
        "window_scaler",
    }
    _require_keys(bundle, required, "contamination")
    if bundle["bundle_schema"] != CONTAMINATION_BUNDLE_SCHEMA:
        raise ValueError("Experimental contamination bundle schema is unsupported")
    feature_names = _string_tuple(bundle, "feature_names")
    composition_features = _string_tuple(bundle, "composition_feature_names")
    if (
        feature_names != CONTAMINATION_FEATURES
        or len(feature_names) != EXPECTED_CONTAMINATION_FEATURE_COUNT
    ):
        raise ValueError("Experimental contamination feature schema differs")
    if (
        composition_features != COMPOSITION_FEATURES
        or len(composition_features) != EXPECTED_COMPOSITION_FEATURE_COUNT
    ):
        raise ValueError("Experimental composition feature schema differs")
    _validate_metadata(
        bundle,
        expected={
            "applicability": "nucleotide_extended_sensitive_only",
            "model_name": EXPERIMENTAL_CONTAMINATION_STRATEGY,
            "random_seed": 20260911,
            "refit_scope": (
                "all_867_frozen_development_rows_with_saved_held_fold_stack_scores"
            ),
            "reliability": EXPERIMENTAL_CONTAMINATION_RELIABILITY,
            "target": "true_contamination_pct",
            "training_row_count": 867,
        },
    )
    _validate_contamination_estimators(bundle)
    _validate_libraries(bundle)


def _validate_contamination_estimators(bundle: Mapping[str, object]) -> None:
    _validate_window_estimators(bundle)
    _validate_contamination_regressor(bundle)


def _validate_window_estimators(bundle: Mapping[str, object]) -> None:
    scaler = _scaler(bundle, "window_scaler")
    classifier = _classifier(bundle, "window_classifier")
    if scaler.get_params() != StandardScaler().get_params():
        raise ValueError("Experimental composition scaler parameters differ")
    if tuple(scaler.feature_names_in_) != COMPOSITION_FEATURES:
        raise ValueError("Experimental composition scaler feature order differs")
    if tuple(classifier.classes_) != (0, 1):
        raise ValueError("Experimental composition classifier classes differ")
    expected_classifier = LogisticRegression(
        C=1.0,
        solver="lbfgs",
        max_iter=1_000,
    ).get_params()
    if classifier.get_params() != expected_classifier:
        raise ValueError("Experimental composition classifier parameters differ")


def _validate_contamination_regressor(bundle: Mapping[str, object]) -> None:
    pipeline = _pipeline(bundle, "regressor_pipeline")
    if tuple(pipeline.named_steps) != ("imputer", "regressor"):
        raise TypeError("Experimental contamination pipeline steps differ")
    imputer = pipeline.named_steps["imputer"]
    _validate_median_imputer(imputer, "contamination")
    regressor = pipeline.named_steps["regressor"]
    if type(regressor) is not ExtraTreesRegressor:
        raise TypeError("Experimental contamination regressor type differs")
    expected_regressor = ExtraTreesRegressor(
        n_estimators=300,
        random_state=20260911,
        n_jobs=1,
    ).get_params()
    if regressor.get_params() != expected_regressor:
        raise ValueError("Experimental contamination regressor parameters differ")
    if regressor.n_features_in_ != EXPECTED_CONTAMINATION_FEATURE_COUNT:
        raise ValueError("Experimental contamination fitted feature count differs")


def _validate_median_imputer(imputer: object, name: str) -> None:
    if type(imputer) is not SimpleImputer:
        raise TypeError(f"Experimental {name} imputer type differs")
    observed = imputer.get_params()
    expected = SimpleImputer(
        strategy="median",
        keep_empty_features=True,
    ).get_params()
    observed_missing = observed.pop("missing_values")
    expected_missing = expected.pop("missing_values")
    missing_values_match = bool(
        np.isnan(observed_missing) and np.isnan(expected_missing)
    )
    if observed != expected or not missing_values_match:
        raise ValueError(f"Experimental {name} imputer parameters differ")


def _validate_libraries(bundle: Mapping[str, object]) -> None:
    libraries = bundle["libraries"]
    if not isinstance(libraries, dict):
        raise TypeError("Experimental model library metadata is not a mapping")
    expected = _library_versions()
    if libraries != expected:
        raise ValueError(
            f"Experimental model library versions differ: expected {libraries}, "
            f"running {expected}"
        )


def _validate_metadata(
    bundle: Mapping[str, object],
    *,
    expected: Mapping[str, object],
) -> None:
    metadata = bundle["metadata"]
    if not isinstance(metadata, dict) or not _is_standard_data(metadata):
        raise TypeError("Experimental model metadata is not standard data")
    mismatches = [
        name
        for name, expected_value in expected.items()
        if metadata.get(name) != expected_value
    ]
    if mismatches:
        raise ValueError(
            f"Experimental model metadata differs: {', '.join(sorted(mismatches))}"
        )


def _is_standard_data(value: object) -> bool:
    if value is None or isinstance(value, (str, int, float, bool)):
        return True
    if isinstance(value, (list, tuple)):
        return all(_is_standard_data(item) for item in value)
    if isinstance(value, dict):
        return all(
            isinstance(key, str) and _is_standard_data(item)
            for key, item in value.items()
        )
    return False


def _library_versions() -> dict[str, str]:
    return {
        "joblib": version("joblib"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "python": platform.python_version(),
        "scikit-learn": version("scikit-learn"),
    }


def _require_keys(bundle: Mapping[str, object], expected: set[str], name: str) -> None:
    if set(bundle) != expected:
        raise ValueError(f"Experimental {name} model bundle keys differ")


def _string_tuple(bundle: Mapping[str, object], key: str) -> tuple[str, ...]:
    value = bundle[key]
    if not isinstance(value, (list, tuple)) or any(
        not isinstance(item, str) for item in value
    ):
        raise TypeError(f"Experimental model field {key} is not a string sequence")
    return tuple(value)


def _string_mapping(bundle: Mapping[str, object], key: str) -> dict[str, str]:
    value = bundle[key]
    if not isinstance(value, dict) or any(
        not isinstance(name, str) or not isinstance(item, str)
        for name, item in value.items()
    ):
        raise TypeError(f"Experimental model field {key} is not a string mapping")
    return value


def _pipeline(bundle: Mapping[str, object], key: str) -> Pipeline:
    value = bundle[key]
    if type(value) is not Pipeline:
        raise TypeError(
            f"Experimental model field {key} is not a plain sklearn Pipeline"
        )
    return value


def _scaler(bundle: Mapping[str, object], key: str) -> StandardScaler:
    value = bundle[key]
    if type(value) is not StandardScaler:
        raise TypeError(f"Experimental model field {key} is not a StandardScaler")
    return value


def _classifier(bundle: Mapping[str, object], key: str) -> LogisticRegression:
    value = bundle[key]
    if type(value) is not LogisticRegression:
        raise TypeError(f"Experimental model field {key} is not LogisticRegression")
    return value


def _require_columns(frame: object, expected: tuple[str, ...], name: str) -> None:
    columns = getattr(frame, "columns", None)
    if columns is None or tuple(columns) != expected:
        raise ValueError(f"Experimental {name} runtime feature order differs")


def _bounded_prediction(pipeline: Pipeline, features: object) -> float:
    prediction = float(pipeline.predict(features)[0])
    if not np.isfinite(prediction):
        raise ValueError("Experimental quality model returned a nonfinite prediction")
    return float(np.clip(prediction, 0.0, 100.0))


def _applicability(
    query_fna: Path,
    *,
    mode_fast: bool,
    sensitive_mode: bool,
) -> QualityPrediction | None:
    if mode_fast:
        return QualityPrediction(
            None,
            "unsupported_fast_mode",
            "The development feature distribution used extended marker trees.",
        )
    if not sensitive_mode:
        return QualityPrediction(
            None,
            "unsupported_nonsensitive_search",
            "The development feature distribution used sensitive HMM search.",
        )
    if not query_fna.is_file():
        return QualityPrediction(
            None,
            "unsupported_protein_input",
            "The development rows and composition model require nucleotide input.",
        )
    return None
