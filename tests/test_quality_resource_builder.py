from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from scripts import benchmark_contamination_methods as fna_benchmark
from scripts import benchmark_contamination_methods_faa as faa_benchmark
from scripts._contamination_features import (
    COMPLETENESS_MODE,
    CONTAMINATION_FEATURE_PROFILE,
    compute_query_features,
    contamination_feature_vector,
    validate_source_split_contract,
)
from src.bin.build_novelty_completeness_resources import (
    build_training_df,
    build_resource_bundle,
    make_validation_plan,
    reference_resources_for_ids,
    simulate_samples,
    train_models,
)
from src.core.contamination_scoring import CONTAMINATION_MODEL_FEATURES
from src.core.summarize.order_completeness import _normalize_order_completeness

ORDER_MARKERS = {"TestOrder": ["m1", "m2", "m3", "m4"]}


def _profiles():
    labels = {}
    profiles = {}
    for family in ("FamilyA", "FamilyB", "FamilyC"):
        for index in range(2):
            genome_id = f"{family}_{index}"
            labels[genome_id] = {
                "domain": "NCLDV",
                "phylum": "TestPhylum",
                "class": "TestClass",
                "order": "TestOrder",
                "family": family,
                "genus": "",
                "species": "",
            }
            profiles[genome_id] = {marker: 1 for marker in ORDER_MARKERS["TestOrder"]}
    return labels, profiles


def _write_compact_resources(resources: Path, *, include_profiles: bool = True) -> None:
    labels, profiles = _profiles()
    (resources / "markers").mkdir(parents=True)
    (resources / "parquet" / "labels").mkdir(parents=True)
    (resources / "markers" / "order_completeness.tab").write_text(
        "Order\tAverage_Percent\tStd_Percent\tOrthogroups\n"
        "TestOrder\t100\t0\tm1,m2,m3,m4\n"
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "label_id": genome_id,
                    "taxonomy": "|".join(labels[genome_id].values()),
                }
                for genome_id in sorted(labels)
            ]
        ),
        resources / "parquet" / "labels" / "labels.parquet",
    )
    if include_profiles:
        pq.write_table(
            pa.Table.from_pylist(
                [
                    {
                        "marker": marker,
                        "seq_id": f"{genome_id}|{marker}_protein",
                        "sequence": "M" * 20,
                    }
                    for genome_id in sorted(profiles)
                    for marker in ORDER_MARKERS["TestOrder"]
                ]
            ),
            resources / "parquet" / "faa.parquet",
        )


def test_compact_resource_builder_writes_nonempty_staged_bundle(tmp_path: Path) -> None:
    resources = tmp_path / "resources"
    output = tmp_path / "built"
    _write_compact_resources(resources)

    groups, models = build_resource_bundle(resources, output, smoke=True)

    assert groups > 0
    assert models == 1
    assert (output / "completeness" / "model.joblib").stat().st_size > 0
    assert (output / "completeness" / "tiers.tsv").stat().st_size > 0
    assert (output / "markers" / "marker_conservation.tsv").stat().st_size > 0
    assert not (resources / "labels.tsv").exists()
    assert not (resources / "database").exists()
    assert not (resources / ".gvclass_cache").exists()

    metadata = pd.read_csv(output / "completeness" / "model_metadata.tsv", sep="\t")
    row = metadata.loc[metadata["order_name"] == "TestOrder"].iloc[0]
    assert row["training_target"] == "marker_dropout_fraction_of_source_profile"
    assert not bool(row["target_is_nucleotide_completeness"])
    assert row["model_reliability"] == "advisory_only"
    assert pd.notna(row["marker_dropout_final_test_r2"])
    assert row["test_r2"] == row["marker_dropout_final_test_r2"]
    assert row["test_metric_scope"] == (
        "marker_dropout_only_not_nucleotide_completeness"
    )
    assert row["marker_panel_ascertainment"] == (
        "fixed_input_panel_predates_validation_holdouts"
    )

    splits = pd.read_csv(output / "completeness" / "validation_splits.tsv", sep="\t")
    final = splits[splits["stage"] == "final_evaluation"]
    development = set(final.loc[final["role"] == "development", "genome_id"])
    final_test = set(final.loc[final["role"] == "final_test", "genome_id"])
    assert development.isdisjoint(final_test)
    selection = splits[splits["stage"] == "candidate_selection"]
    assert selection["fold"].nunique() == int(row["n_selection_folds"])

    labels, profiles = _profiles()
    _, _, tiers, baselines, strategy1_baseline = reference_resources_for_ids(
        "TestOrder",
        profiles,
        list(development),
        labels,
        ORDER_MARKERS,
    )
    training_df = build_training_df(
        labels,
        {"TestOrder": {genome_id: profiles[genome_id] for genome_id in development}},
        ORDER_MARKERS,
        tiers,
        baselines,
        {"TestOrder": strategy1_baseline},
    )
    runtime_table = pd.read_csv(output / "markers" / "order_completeness.tab", sep="\t")
    sample = training_df.iloc[0]
    runtime_score, _, _ = _normalize_order_completeness(
        runtime_table, float(sample["raw_score"]), "TestOrder"
    )
    assert runtime_score == pytest.approx(float(sample["strategy1_score"]))
    assert runtime_table.iloc[0]["Baseline_Source"] == (
        "development_only_reference_profiles"
    )


def test_zero_profiles_fail_before_output_and_leave_source_intact(
    tmp_path: Path,
) -> None:
    resources = tmp_path / "resources"
    output = tmp_path / "built"
    _write_compact_resources(resources, include_profiles=False)
    before = {
        path.relative_to(resources): path.read_bytes()
        for path in resources.rglob("*")
        if path.is_file()
    }

    with pytest.raises(RuntimeError, match="No reference marker profiles"):
        build_resource_bundle(resources, output, smoke=True)

    after = {
        path.relative_to(resources): path.read_bytes()
        for path in resources.rglob("*")
        if path.is_file()
    }
    assert after == before
    assert not output.exists()


def test_grouped_folds_exclude_duplicate_sources_and_held_family_statistics() -> None:
    labels, profiles = _profiles()
    plan = make_validation_plan(profiles, labels)
    assert plan is not None
    assert plan["split_unit"] == "family"

    development = set(plan["development_ids"])
    final_test = set(plan["final_test_ids"])
    assert development.isdisjoint(final_test)
    _, development_family_profiles, _, _, _ = reference_resources_for_ids(
        "TestOrder",
        profiles,
        list(development),
        labels,
        ORDER_MARKERS,
    )
    final_families = {labels[genome_id]["family"] for genome_id in final_test}
    assert final_families.isdisjoint(
        {family for _, family in development_family_profiles}
    )
    for fold_index, fold in enumerate(plan["selection_folds"]):
        training_ids = set(fold["training_ids"])
        validation_ids = set(fold["validation_ids"])
        assert training_ids.isdisjoint(validation_ids)
        assert final_test.isdisjoint(training_ids | validation_ids)

        training_samples = simulate_samples(
            labels,
            {"TestOrder": {gid: profiles[gid] for gid in training_ids}},
            ORDER_MARKERS,
            seed=fold_index,
        )
        validation_samples = simulate_samples(
            labels,
            {"TestOrder": {gid: profiles[gid] for gid in validation_ids}},
            ORDER_MARKERS,
            seed=fold_index + 100,
        )
        assert {row["source_genome_id"] for row in training_samples}.isdisjoint(
            {row["source_genome_id"] for row in validation_samples}
        )

        _, family_profiles, _, _, _ = reference_resources_for_ids(
            "TestOrder",
            profiles,
            list(training_ids),
            labels,
            ORDER_MARKERS,
        )
        held_families = {labels[genome_id]["family"] for genome_id in validation_ids}
        derived_families = {family for _, family in family_profiles}
        assert held_families.isdisjoint(derived_families)

    unbalanced_labels = {}
    unbalanced_profiles = {}
    for family, size in (("FamilyA", 2), ("FamilyB", 2), ("FamilyZ", 6)):
        for index in range(size):
            genome_id = f"{family}_{index}"
            unbalanced_labels[genome_id] = {"family": family}
            unbalanced_profiles[genome_id] = {"m1": 1}
    balanced_plan = make_validation_plan(unbalanced_profiles, unbalanced_labels)
    assert balanced_plan is not None
    assert balanced_plan["final_held_out_group"] != "FamilyZ"
    assert len(balanced_plan["final_test_ids"]) == 2


def test_insufficient_source_groups_are_unsupported_not_random_split() -> None:
    labels, profiles = _profiles()
    profiles = dict(list(profiles.items())[:3])
    labels = {genome_id: labels[genome_id] for genome_id in profiles}

    models, metadata, tiers, *_ = train_models(
        labels, {"TestOrder": profiles}, ORDER_MARKERS, smoke=True
    )

    assert models == {}
    assert ("order", "TestOrder", "TestOrder") in tiers
    row = metadata.iloc[0]
    assert row["training_status"] == "unsupported"
    assert row["validation_type"] == "unsupported_grouped_validation"


class _FakeSummarizer:
    completeness_mode = COMPLETENESS_MODE
    sensitive_mode = True

    def summarize_query_full(self, query_id, query_output_dir, tree_nn, mode_fast=True):
        result = {
            feature: float(index)
            for index, feature in enumerate(CONTAMINATION_MODEL_FEATURES)
        }
        result["contamination_feature_version"] = CONTAMINATION_FEATURE_PROFILE
        return result


def test_benchmark_features_match_runtime_vector_and_settings(tmp_path: Path) -> None:
    summary = compute_query_features(
        "query",
        tmp_path,
        tmp_path,
        summarizer=_FakeSummarizer(),
        sensitive_mode=True,
        completeness_mode=COMPLETENESS_MODE,
    )

    vector = contamination_feature_vector(summary)
    assert list(vector) == CONTAMINATION_MODEL_FEATURES
    assert fna_benchmark.METHOD3_FEATURES == CONTAMINATION_MODEL_FEATURES
    assert faa_benchmark.METHOD3_FEATURES == CONTAMINATION_MODEL_FEATURES
    assert len(vector) == 21
    assert summary["benchmark_completeness_mode"] == "novelty-aware"
    assert summary["benchmark_sensitive_mode"] is True
    assert summary["benchmark_mode_fast"] is True
    assert summary["contamination_feature_profile"] == CONTAMINATION_FEATURE_PROFILE


@pytest.mark.parametrize("benchmark", [fna_benchmark, faa_benchmark])
def test_benchmark_run_profile_matches_feature_recompute(
    benchmark, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    results = tmp_path / "results"
    calls = []

    def fake_run(command, check):
        calls.append(command)
        results.mkdir()
        (results / "gvclass_summary.tsv").write_text("query\n")

    monkeypatch.setattr(benchmark.subprocess, "run", fake_run)
    benchmark.run_gvclass(tmp_path / "inputs", results, workers=2, force=False)

    assert calls[0][calls[0].index("--completeness-mode") + 1] == COMPLETENESS_MODE
    assert "--mode-fast" in calls[0]
    assert "--sensitive" in calls[0]
    profile = json.loads((results / "benchmark_run_profile.json").read_text())
    assert profile["completeness_mode"] == COMPLETENESS_MODE
    assert profile["sensitive_mode"] is True
    assert profile["mode_fast"] is True
    assert profile["contamination_feature_profile"] == CONTAMINATION_FEATURE_PROFILE
    assert profile["contamination_model_feature_names"] == (
        CONTAMINATION_MODEL_FEATURES
    )

    benchmark.run_gvclass(tmp_path / "inputs", results, workers=2, force=False)
    assert len(calls) == 1
    (results / "benchmark_run_profile.json").unlink()
    with pytest.raises(RuntimeError, match="no matching run profile"):
        benchmark.run_gvclass(tmp_path / "inputs", results, workers=2, force=False)


def test_benchmark_export_fits_train_only_and_records_final_test(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rows = []
    for split, source, target in (
        ("train", "train_a", 0.0),
        ("train", "train_b", 20.0),
        ("test", "final_c", 90.0),
    ):
        row = {
            "split": split,
            "base_accession": source,
            "true_contamination_pct": target,
            "contamination_present": int(target > 0),
        }
        row.update({feature: 0.0 for feature in CONTAMINATION_MODEL_FEATURES})
        rows.append(row)
    feature_df = pd.DataFrame(rows)

    class RecordingModel:
        def fit(self, X, y):
            self.fit_targets = list(y)
            return self

        def predict(self, X):
            return pd.Series([0.0] * len(X))

    model = RecordingModel()
    captured = {}
    monkeypatch.setattr(faa_benchmark, "create_model_candidate", lambda name: model)
    monkeypatch.setattr(
        faa_benchmark.joblib,
        "dump",
        lambda bundle, path: captured.update(bundle=bundle, path=path),
    )

    faa_benchmark.export_model(feature_df, str(tmp_path / "model.joblib"), "test")

    assert model.fit_targets == [0.0, 20.0]
    assert captured["bundle"]["fit_source_ids"] == ["train_a", "train_b"]
    assert captured["bundle"]["held_out_source_ids"] == ["final_c"]
    assert captured["bundle"]["feature_names"] == CONTAMINATION_MODEL_FEATURES

    with pytest.raises(ValueError, match="Missing source split columns"):
        validate_source_split_contract(feature_df.drop(columns=["split"]))
