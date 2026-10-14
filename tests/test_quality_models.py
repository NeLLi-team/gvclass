"""Contracts for the experimental quality-model features and strict loaders."""

from __future__ import annotations

import hashlib
import io
from pathlib import Path
from unittest.mock import Mock

import joblib
import numpy as np
import pytest
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from src.core import quality_models
from src.core.quality_features import (
    COMPOSITION_FEATURES,
    build_completeness_features,
    build_composition_windows,
    build_family_copy_counts,
)


def test_family_copy_counts_unions_proteins_across_member_models() -> None:
    marker_hits = {
        "model_a": {"protein_1", "protein_2"},
        "model_b": {"protein_2", "protein_3"},
        "review": {"protein_4"},
    }

    counts = build_family_copy_counts(
        marker_hits,
        {"model_a": "group::one", "model_b": "group::one"},
        ("review",),
    )

    assert counts == {"group::one": 3}


def test_family_copy_counts_rejects_model_outside_frozen_catalogue() -> None:
    with pytest.raises(ValueError, match="absent from the frozen catalogue"):
        build_family_copy_counts(
            {"unknown": {"protein_1"}},
            {"known": "model::known"},
            (),
        )


def test_completeness_features_reject_missing_summary_field() -> None:
    with pytest.raises(KeyError, match="LENbp"):
        build_completeness_features(
            {},
            {},
            family_vocabulary=(),
            model_to_family={},
            review_models=(),
        )


def test_completeness_scorer_requires_a_filtered_marker_protein() -> None:
    scorer = quality_models.ExperimentalQualityScorer({}, {})

    prediction = scorer._score_completeness(
        {},
        {},
        {"order": "Imitervirales"},
    )

    assert prediction.value is None
    assert prediction.status == "insufficient_evidence"


def test_completeness_scorer_requires_an_accepted_order() -> None:
    scorer = quality_models.ExperimentalQualityScorer({}, {})

    prediction = scorer._score_completeness(
        {},
        {"GVOGm0001": {"protein_1"}},
        {"order": ""},
    )

    assert prediction.value is None
    assert prediction.status == "unassigned"


def test_composition_windows_preserve_reverse_complement_features(
    tmp_path: Path,
) -> None:
    query = tmp_path / "query.fna"
    query.write_text(">forward\n" + "A" * 2_000 + "\n>reverse\n" + "T" * 2_000 + "\n")

    windows = build_composition_windows(query)

    assert windows.status == "available"
    assert windows.coverage_pct == 100.0
    assert tuple(windows.features.columns) == COMPOSITION_FEATURES
    np.testing.assert_array_equal(windows.features.iloc[0], windows.features.iloc[1])


def test_composition_windows_preserve_coverage_failure(tmp_path: Path) -> None:
    query = tmp_path / "query.fna"
    query.write_text(">query\n" + "A" * 2_000 + "N" * 2_000 + "\n")

    windows = build_composition_windows(query)

    assert windows.status == "coverage_below_80_pct"
    assert windows.scored_window_count == 1
    assert windows.unscored_full_window_bp == 2_000
    assert windows.coverage_pct == 50.0


@pytest.mark.parametrize(
    ("mode_fast", "sensitive_mode", "has_fna", "expected"),
    [
        (True, True, True, "unsupported_fast_mode"),
        (False, False, True, "unsupported_nonsensitive_search"),
        (False, True, False, "unsupported_protein_input"),
    ],
)
def test_candidate_pair_preserves_profile_applicability(
    tmp_path: Path,
    mode_fast: bool,
    sensitive_mode: bool,
    has_fna: bool,
    expected: str,
) -> None:
    query_fna = tmp_path / "query.fna"
    if has_fna:
        query_fna.write_text(">query\nACGT\n")
    scorer = quality_models.ExperimentalQualityScorer({}, {})

    predictions = scorer.score_query(
        summary={},
        marker_hits={},
        query_fna=query_fna,
        query_output_dir=tmp_path,
        tree_nn_results={},
        accepted_lineage={},
        mode_fast=mode_fast,
        sensitive_mode=sensitive_mode,
        contamination_scorer=Mock(),
    )

    assert predictions.completeness.status == expected
    assert predictions.contamination.status == expected
    assert predictions.completeness.value is None
    assert predictions.contamination.value is None


@pytest.mark.parametrize("domain", ["EUK", "BAC", "PHAGE", "MIRUS", ""])
def test_candidate_pair_rejects_non_ncldv_domain_before_features(
    tmp_path: Path,
    domain: str,
) -> None:
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(">query\nACGT\n")
    scorer = quality_models.ExperimentalQualityScorer({}, {})

    predictions = scorer.score_query(
        summary={},
        marker_hits={"GVOGm0001": {"protein_1"}},
        query_fna=query_fna,
        query_output_dir=tmp_path,
        tree_nn_results={},
        accepted_lineage={"domain": domain},
        mode_fast=False,
        sensitive_mode=True,
        contamination_scorer=Mock(),
    )

    assert predictions.completeness.status == "unsupported_target_domain"
    assert predictions.contamination.status == "unsupported_target_domain"
    assert predictions.completeness.value is None
    assert predictions.contamination.value is None


def test_model_profile_is_stable_and_binds_settings() -> None:
    first = quality_models.experimental_model_profile(
        mode_fast=False,
        sensitive_mode=True,
    )
    repeated = quality_models.experimental_model_profile(
        mode_fast=False,
        sensitive_mode=True,
    )
    changed = quality_models.experimental_model_profile(
        mode_fast=True,
        sensitive_mode=True,
    )

    assert first == repeated
    assert first.startswith("experimental:")
    assert len(first.removeprefix("experimental:")) == 64
    assert first != changed


def test_hash_mismatch_stops_before_joblib_deserialization(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model_path = tmp_path / "candidate.joblib"
    model_path.write_bytes(b"untrusted model bytes")
    load = Mock(side_effect=AssertionError("joblib.load must not run"))
    monkeypatch.setattr(quality_models.joblib, "load", load)

    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        quality_models._load_bundle(model_path, "0" * 64, "test")

    load.assert_not_called()


def test_loader_deserializes_the_same_bytes_that_were_hashed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    trusted_buffer = io.BytesIO()
    joblib.dump({"identity": "trusted"}, trusted_buffer)
    trusted_bytes = trusted_buffer.getvalue()
    replacement_buffer = io.BytesIO()
    joblib.dump({"identity": "replacement"}, replacement_buffer)
    replacement_bytes = replacement_buffer.getvalue()
    model_path = tmp_path / "candidate.joblib"
    model_path.write_bytes(trusted_bytes)
    original_read_bytes = Path.read_bytes

    def replace_after_read(path: Path) -> bytes:
        observed = original_read_bytes(path)
        path.write_bytes(replacement_bytes)
        return observed

    monkeypatch.setattr(Path, "read_bytes", replace_after_read)
    expected_hash = hashlib.sha256(trusted_bytes).hexdigest()

    loaded = quality_models._load_bundle(
        model_path,
        expected_hash,
        "test",
    )

    assert loaded == {"identity": "trusted"}
    assert joblib.load(model_path) == {"identity": "replacement"}


def test_serialized_pipeline_preserves_prediction(tmp_path: Path) -> None:
    features = np.array([[0.0], [1.0], [2.0], [3.0]])
    truth = np.array([0.0, 10.0, 20.0, 30.0])
    pipeline = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("regressor", HistGradientBoostingRegressor(random_state=20260911)),
        ]
    ).fit(features, truth)
    model_path = tmp_path / "candidate.joblib"
    joblib.dump({"pipeline": pipeline}, model_path)
    expected = pipeline.predict(np.array([[1.5]]))

    loaded = joblib.load(model_path)
    observed = loaded["pipeline"].predict(np.array([[1.5]]))

    np.testing.assert_array_equal(observed, expected)


def test_fitted_median_imputer_accepts_nan_missing_value_parameter() -> None:
    imputer = SimpleImputer(
        strategy="median",
        keep_empty_features=True,
    ).fit(np.array([[1.0], [np.nan]]))

    quality_models._validate_median_imputer(imputer, "test")
