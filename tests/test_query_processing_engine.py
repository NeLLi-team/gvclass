"""Stages of the per-query processing engine."""

from __future__ import annotations

import logging
from pathlib import Path
from unittest.mock import Mock

import pytest

from src.pipeline import query_processing_engine
from src.pipeline.query_processing_engine import _initialize_query
from src.utils.error_handling import ProcessingError
from tests.pipeline_helpers import make_config

LOGGER = logging.getLogger(__name__)


def test_query_initialization_rejects_a_working_directory_that_holds_the_input(
    tmp_path: Path,
) -> None:
    query_dir = tmp_path / ".gvclass_work" / "sample"
    query_dir.mkdir(parents=True)
    query = query_dir / "sample.fna"
    query.write_text(">contig1\nACGT\n", encoding="utf-8")

    with pytest.raises(ProcessingError, match="holds the input file"):
        _initialize_query(query, tmp_path)

    assert query.read_text(encoding="utf-8") == ">contig1\nACGT\n"


def test_missing_combined_hmm_fails_the_marker_search(tmp_path):
    protein = tmp_path / "query.faa"
    protein.write_text(">q|p_1\nMKVLAA\n", encoding="utf-8")
    (tmp_path / "db" / "hmm").mkdir(parents=True)
    (tmp_path / "db" / "hmm" / "busco69.hmm").write_text("partial", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match=r"combined\.hmm"):
        query_processing_engine._run_hmm_and_blast(
            "q",
            tmp_path / "out",
            protein,
            make_config(tmp_path),
            1,
            LOGGER,
        )

    assert not (tmp_path / "out" / "hmmout" / "models.out").exists()


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
            query_output_dir=tmp_path,
            config=make_config(tmp_path),
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
            "q1", tmp_path, make_config(tmp_path, mode_fast=False), {}, LOGGER
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
        config=make_config(tmp_path, mode_fast=False),
        threads=1,
        logger=LOGGER,
    )

    summarize.assert_called_once()
