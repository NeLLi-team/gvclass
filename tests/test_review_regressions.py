"""Regression tests for defects found by the independent review of the 2.1.0 launcher.

The launcher exists twice, as the ``gvclass`` script and as
``src.bin.gvclass_apptainer``; each behaviour is checked on both copies.
"""

from __future__ import annotations

import logging
import runpy
import sys
from pathlib import Path
from typing import Any
from unittest.mock import Mock

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.bin import gvclass_apptainer, gvclass_cli
from src.core.contamination_scoring import ContaminationScorer
from src.core.reformat import reformat_sequences
from src.pipeline import query_processing_engine
from src.pipeline.parallel_runner import _resolve_worker_distribution
from src.pipeline.query_processing_engine import _initialize_query
from src.utils.error_handling import ProcessingError

REPO_ROOT = Path(__file__).resolve().parents[1]
LOGGER = logging.getLogger(__name__)


@pytest.fixture(params=["module", "script"])
def launcher(request: pytest.FixtureRequest) -> dict[str, Any]:
    if request.param == "module":
        return vars(gvclass_apptainer)
    return runpy.run_path(str(REPO_ROOT / "gvclass"), run_name="gvclass_launcher")


def test_image_cache_names_differ_for_urls_with_the_same_basename(
    launcher: dict[str, Any],
) -> None:
    name = launcher["cache_filename_for_image"]

    first = name("https://host-a.example/v1/gvclass.sif")
    second = name("https://host-b.example/v2/gvclass.sif")

    assert first != second
    assert first.endswith("_gvclass.sif")
    assert first == name("https://host-a.example/v1/gvclass.sif")


def test_library_image_cache_name_keeps_the_tag(launcher: dict[str, Any]) -> None:
    name = launcher["cache_filename_for_image"](
        "library://group/collection/gvclass:2.1.0"
    )

    assert name == "group_collection_gvclass_2.1.0.sif"


def test_container_run_lets_no_sensitive_win_like_the_native_cli(
    launcher: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: list[list[str]] = []
    run_container = launcher["run_container"]
    monkeypatch.setitem(run_container.__globals__, "ensure_local_image", lambda i: i)
    monkeypatch.setattr("subprocess.call", lambda cmd: captured.append(list(cmd)) or 0)
    query = tmp_path / "query"
    query.mkdir()

    run_container(
        query,
        tmp_path / "out",
        2,
        "/tmp/image.sif",
        tmp_path / "cache",
        sensitive_mode=True,
        no_sensitive=True,
    )

    assert "--no-sensitive" in captured[0]
    assert "--sensitive" not in captured[0]


def test_default_output_name_for_the_current_directory(
    launcher: dict[str, Any], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = tmp_path / "genomes"
    query.mkdir()
    monkeypatch.chdir(query)

    output = launcher["resolve_output_path"](Path(), None, None)

    assert output == tmp_path / "genomes_results"


def test_repo_detection_ignores_an_unrelated_src_directory(tmp_path: Path) -> None:
    script = runpy.run_path(str(REPO_ROOT / "gvclass"), run_name="gvclass_launcher")
    (tmp_path / "src").mkdir()

    assert script["is_gvclass_repo"](tmp_path) is False
    assert script["is_gvclass_repo"](REPO_ROOT) is True


def test_query_initialization_rejects_a_working_directory_that_holds_the_input(
    tmp_path: Path,
) -> None:
    query_dir = tmp_path / "sample"
    query_dir.mkdir()
    query = query_dir / "sample.fna"
    query.write_text(">contig1\nACGT\n", encoding="utf-8")

    with pytest.raises(ProcessingError, match="holds the input file"):
        _initialize_query(query, tmp_path)

    assert query.read_text(encoding="utf-8") == ">contig1\nACGT\n"


def test_query_initialization_keeps_earlier_work_beside_the_new_directory(
    tmp_path: Path,
) -> None:
    query = tmp_path / "input" / "sample.fna"
    query.parent.mkdir()
    query.write_text(">contig1\nACGT\n", encoding="utf-8")
    earlier = tmp_path / "out" / "sample"
    earlier.mkdir(parents=True)
    (earlier / "old.txt").write_text("earlier attempt", encoding="utf-8")

    name, work = _initialize_query(query, tmp_path / "out")

    retained = list((tmp_path / "out").glob("sample.previous.*"))
    assert (name, work) == ("sample", earlier)
    assert list(work.iterdir()) == []
    assert [(path / "old.txt").read_text(encoding="utf-8") for path in retained] == [
        "earlier attempt"
    ]


@pytest.mark.parametrize(
    ("total_threads", "threads_per_worker", "expected"),
    [(2, 8, (1, 2)), (16, 4, (4, 4)), (1, 1, (1, 1))],
)
def test_explicit_threads_per_worker_never_exceeds_the_total(
    total_threads: int, threads_per_worker: int, expected: tuple[int, int]
) -> None:
    plan = _resolve_worker_distribution(8, total_threads, None, threads_per_worker)

    assert plan == expected


def test_single_contig_keeps_one_identity_for_nucleotide_and_protein_records(
    tmp_path: Path,
) -> None:
    fna = tmp_path / "parent__ctg_reformatted.fna"
    records = reformat_sequences(
        [SeqRecord(Seq("ACGT" * 6000), id="ctg")], "parent__ctg.fna"
    )
    SeqIO.write(records, fna, "fasta")
    query_ids = [f"parent__ctg|ctg_{index}" for index in range(1, 6)]
    faa = tmp_path / "parent__ctg_reformatted.faa"
    faa.write_text(
        "".join(f">{name}\n{'M' * 100}\n" for name in query_ids), encoding="utf-8"
    )
    blast = tmp_path / "blastp_out"
    blast.mkdir()
    rows = [
        f"{name}\tBAC__ref|p1\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100"
        for name in query_ids
    ]
    (blast / "OG1.m8").write_text("\n".join(rows) + "\n", encoding="utf-8")
    scorer = ContaminationScorer.__new__(ContaminationScorer)
    scorer.labels = {"BAC__ref": {"domain": "BAC", "order": "OrderA", "family": "FamA"}}
    tree_neighbors = {"OG1": {name: {"BAC__ref|p1": 0.1} for name in query_ids}}

    result = scorer.collect_contig_features(fna, faa, blast, "", "", tree_neighbors)

    assert result["contig_count"] == 1
    assert result["suspicious_bp_fraction"] == 100.0
    assert result["cellular_coherent_bp_fraction"] == 100.0


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
    )

    summarize.assert_called_once()


def test_container_mode_without_a_query_is_an_error(
    launcher: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    entry = launcher.get("_standalone_apptainer_main") or launcher["main"]

    with pytest.raises(SystemExit) as raised:
        entry(["-o", "results", "--image", "gvclass.sif"])

    assert raised.value.code == 2
    assert "query path is required" in capsys.readouterr().err


def test_container_mode_without_arguments_prints_help(
    launcher: dict[str, Any], capsys: pytest.CaptureFixture[str]
) -> None:
    entry = launcher.get("_standalone_apptainer_main") or launcher["main"]

    assert entry([]) == 0
    assert "usage:" in capsys.readouterr().out


def test_abbreviated_container_flag_runs_the_container(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    called: list[list[str]] = []
    monkeypatch.setattr(
        gvclass_apptainer, "main", lambda argv: called.append(list(argv)) or 0
    )
    monkeypatch.setattr(sys, "argv", ["gvclass", "queries", "--imag", "x.sif"])

    assert gvclass_cli.main() == 0
    assert called == [["queries", "--imag", "x.sif"]]


def test_unrecognized_native_arguments_still_fail(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sys, "argv", ["gvclass", "queries", "--no-such-option"])

    with pytest.raises(SystemExit) as raised:
        gvclass_cli.main()

    assert raised.value.code == 2
