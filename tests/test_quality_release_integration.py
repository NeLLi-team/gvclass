"""Release contracts for explicit model selection and unavailable evidence."""

import csv
import logging
import math
from unittest.mock import Mock

import pytest

from src.bin import gvclass_cli
from src.core.quality_models import CandidateQualityPredictions, QualityPrediction
from src.core.summarize_full import _apply_experimental_predictions
from src.pipeline import parallel_runner, query_processing_engine, summary_writer
from src.utils.database_manager import DatabaseManager
from src.utils.error_handling import ProcessingError

LOGGER = logging.getLogger(__name__)


@pytest.fixture
def candidate_summary():
    return {
        "query": "q1",
        "quality_models": "experimental",
        "quality_model_profile": "candidate-a",
        "completeness_query_support": "supported",
        "estimated_completeness_strategy": "experimental_family_union_hgb_v1",
        "completeness_model_group": "development_panel",
        "completeness_model_reliability": "development_only",
        "contamination_query_support": "supported",
        "estimated_contamination_strategy": "experimental_composition_regional_et_v1",
        "contamination_model_reliability": "development_only",
    }


def test_experimental_amount_does_not_assign_contamination_class():
    summary = {"contamination_flag_v1": "none"}
    predictions = CandidateQualityPredictions(
        completeness=QualityPrediction(72.5, "supported", ""),
        contamination=QualityPrediction(25.0, "supported", ""),
    )

    _apply_experimental_predictions(summary, predictions)

    assert summary["estimated_completeness"] == 72.5
    assert summary["estimated_contamination"] == 25.0
    assert summary["contamination_type"] == "uncertain"
    assert summary["contamination_flag_v1"] == "none"
    assert (
        summary["contamination_prediction_range_status"] == "uncalibrated_experimental"
    )


def test_unavailable_candidate_is_missing_and_keeps_reason():
    predictions = CandidateQualityPredictions(
        completeness=QualityPrediction(
            None, "unsupported_input", "nucleotide required"
        ),
        contamination=QualityPrediction(
            None, "unsupported_input", "nucleotide required"
        ),
    )
    summary = {}

    _apply_experimental_predictions(summary, predictions)

    assert math.isnan(summary["estimated_completeness"])
    assert math.isnan(summary["estimated_contamination"])
    assert summary["contamination_query_support"] == "unsupported_input"
    assert summary["experimental_contamination_detail"] == "nucleotide required"


@pytest.mark.parametrize(
    ("models", "profile", "complete"),
    [
        ("experimental", "candidate-a", True),
        ("experimental", "candidate-b", False),
        ("bundled", "bundled", False),
    ],
)
def test_resume_requires_saved_model_profile(
    tmp_path, candidate_summary, models, profile, complete
):
    saved = summary_writer.write_individual_final_summary_file(
        tmp_path, "q1", candidate_summary, LOGGER
    )
    (tmp_path / "q1.SUCCESS").touch()

    with saved.open(encoding="utf-8") as handle:
        fields = csv.DictReader(handle, delimiter="\t").fieldnames
    assert len(fields) == len(set(fields))
    assert (
        parallel_runner._query_is_resume_complete(
            tmp_path, "q1", LOGGER, models, profile
        )
        is complete
    )


@pytest.mark.parametrize(
    ("mode_fast", "sensitive_mode"), [(True, True), (False, False)]
)
def test_experimental_pipeline_rejects_unsupported_search_mode(
    mode_fast, sensitive_mode
):
    with pytest.raises(ValueError, match="require extended mode and sensitive"):
        parallel_runner.gvclass_flow(
            "unused-input",
            "unused-output",
            quality_models="experimental",
            mode_fast=mode_fast,
            sensitive_mode=sensitive_mode,
        )


def test_cli_installs_effective_pin_instead_of_zenodo_latest(tmp_path, monkeypatch):
    pinned = {
        "version": "v2.0.0",
        "url": "https://example.org/pinned.tar.gz",
        "filename": "pinned.tar.gz",
        "sha256": "a" * 64,
    }
    monkeypatch.setattr(
        gvclass_cli,
        "inspect_local_database_state",
        Mock(return_value={"exists": False, "complete": False, "version": "missing"}),
    )
    monkeypatch.setattr(
        DatabaseManager,
        "get_latest_database_source",
        Mock(return_value={"version": "v9.0.0"}),
    )
    monkeypatch.setattr(
        DatabaseManager, "get_effective_database_source", Mock(return_value=pinned)
    )
    install = Mock(return_value=tmp_path)
    monkeypatch.setattr(DatabaseManager, "setup_database", install)

    assert gvclass_cli.check_and_setup_database(
        tmp_path, {}, gvclass_cli.CliOutput(plain_output=True)
    )

    assert install.call_args.kwargs["preferred_source"] == pinned


def test_partial_hmm_set_fails_without_combined_resource(tmp_path):
    hmm_dir = tmp_path / "hmm"
    hmm_dir.mkdir()
    (hmm_dir / "busco69.hmm").write_text("partial", encoding="utf-8")

    with pytest.raises(ProcessingError, match="incomplete individual HMM set"):
        query_processing_engine._resolve_hmm_files(tmp_path)


def test_partial_hmm_set_uses_complete_combined_resource(tmp_path):
    hmm_dir = tmp_path / "hmm"
    hmm_dir.mkdir()
    (hmm_dir / "busco69.hmm").write_text("partial", encoding="utf-8")
    combined = hmm_dir / "combined.hmm"
    combined.write_text("combined", encoding="utf-8")

    assert query_processing_engine._resolve_hmm_files(tmp_path) == [str(combined)]


def test_marker_failure_propagates_from_thread_worker(tmp_path, monkeypatch):
    monkeypatch.setattr(
        query_processing_engine,
        "_process_single_marker",
        Mock(side_effect=RuntimeError("similarity search failed")),
    )

    with pytest.raises(RuntimeError, match="similarity search failed"):
        query_processing_engine._run_marker_tasks(
            query_name="q1",
            marker_files=[("marker", tmp_path / "query.faa")],
            database_path=tmp_path,
            query_output_dir=tmp_path,
            tree_method="veryfasttree",
            max_parallel=1,
            threads_per_marker=1,
            logger=LOGGER,
        )


def test_summary_failure_cannot_become_successful_fallback(tmp_path, monkeypatch):
    monkeypatch.setattr(
        query_processing_engine,
        "FullSummarizer",
        Mock(side_effect=ValueError("invalid required labels")),
    )

    with pytest.raises(ValueError, match="invalid required labels"):
        query_processing_engine._generate_summary_data(
            "q1", tmp_path, tmp_path, {}, False, "novelty-aware", True, LOGGER
        )


def test_tree_analysis_uses_only_current_attempt_paths(tmp_path, monkeypatch):
    tree_dir = tmp_path / "queryrefs_genetrees"
    tree_dir.mkdir()
    current = tree_dir / "current.treefile"
    current.write_text("(q1:1,reference:1);", encoding="utf-8")
    stale = tree_dir / "stale.treefile"
    stale.write_text("(q1:1,old_reference:1);", encoding="utf-8")
    analyze = Mock(return_value={})
    monkeypatch.setattr(query_processing_engine, "_run_tree_analysis", analyze)

    query_processing_engine._analyze_marker_trees(
        "q1", tmp_path, tmp_path, {"current": {"tree": current}}, LOGGER
    )

    assert analyze.call_args.args[4] == [current]


def test_retry_preserves_old_evidence_and_starts_an_empty_directory(tmp_path):
    old_work = tmp_path / "q1"
    old_blast = old_work / "blastp_out"
    old_blast.mkdir(parents=True)
    (old_blast / "old.m8").write_text("previous evidence", encoding="utf-8")

    query_name, current_work = query_processing_engine._initialize_query(
        tmp_path / "q1.fna", tmp_path
    )

    assert query_name == "q1"
    assert current_work == old_work
    assert list(current_work.iterdir()) == []
    retained = list(tmp_path.glob("q1.previous.*"))
    assert len(retained) == 1
    assert (retained[0] / "blastp_out/old.m8").read_text(
        encoding="utf-8"
    ) == "previous evidence"


def test_canonical_gene_outputs_exist_before_quality_scoring(tmp_path, monkeypatch):
    gene_dir = tmp_path / "gene_calling/q1_reformatted"
    gene_dir.mkdir(parents=True)
    nucleotide = gene_dir / "q1_reformatted.fna"
    nucleotide.write_text(">q1_1\nACGTACGT\n", encoding="utf-8")
    protein = gene_dir / "q1_reformatted.faa"
    protein.write_text(">q1_1_1\nMK\n", encoding="utf-8")
    gff = gene_dir / "q1_reformatted.gff"
    gff.write_text("current coordinates\n", encoding="utf-8")
    prepared = query_processing_engine.PreparedQueryInput(
        protein_file=protein,
        reformatted_file=nucleotide,
        best_code=11,
        is_nucleotide=True,
    )
    monkeypatch.setattr(
        query_processing_engine, "_prepare_query_input", Mock(return_value=prepared)
    )
    monkeypatch.setattr(query_processing_engine, "_run_hmm_and_blast", Mock())
    monkeypatch.setattr(
        query_processing_engine, "_process_markers", Mock(return_value={})
    )

    def inspect_current_files(*_args: object, **_kwargs: object) -> dict[str, object]:
        assert (tmp_path / "query_gff/q1.gff").read_bytes() == gff.read_bytes()
        assert (tmp_path / "query_fna/q1.fna").read_bytes() == nucleotide.read_bytes()
        return {}

    summarize = Mock(side_effect=inspect_current_files)
    monkeypatch.setattr(query_processing_engine, "_generate_summary_data", summarize)

    query_processing_engine._run_query_stages(
        query_file=tmp_path / "q1.fna",
        query_name="q1",
        query_output_dir=tmp_path,
        database_path=tmp_path,
        genetic_codes=[11],
        tree_method="veryfasttree",
        mode_fast=False,
        completeness_mode="novelty-aware",
        sensitive_mode=True,
        threads=1,
        logger=LOGGER,
        quality_models="experimental",
    )

    summarize.assert_called_once()


@pytest.mark.parametrize(
    "writer",
    [
        summary_writer.write_individual_summary_file,
        summary_writer.write_individual_final_summary_file,
    ],
)
def test_required_summary_write_failure_raises(tmp_path, candidate_summary, writer):
    with pytest.raises(FileNotFoundError):
        writer(tmp_path / "missing", "q1", candidate_summary, LOGGER)
