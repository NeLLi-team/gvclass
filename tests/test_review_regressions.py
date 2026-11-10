"""Regression tests for defects found by the independent review of 2.1.0."""

from __future__ import annotations

import logging
from pathlib import Path
from types import ModuleType
from unittest.mock import Mock

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.bin import gvclass_cli
from src.core.contamination_scoring import ContaminationScorer
from src.core.labels import LabelIndex
from src.core.reformat import reformat_sequences
from src.core.summarize.paths_stats import QueryEvidence
from src.pipeline import query_processing_engine
from src.pipeline.query_processing_engine import _initialize_query
from src.pipeline.run_config import plan_workers
from src.utils.error_handling import ProcessingError
from tests.launcher_helpers import load_launcher
from tests.pipeline_helpers import make_config

REPO_ROOT = Path(__file__).resolve().parents[1]
LOGGER = logging.getLogger(__name__)


@pytest.fixture
def launcher() -> ModuleType:
    return load_launcher()


def test_image_cache_names_differ_for_urls_with_the_same_basename(
    launcher: ModuleType,
) -> None:
    name = launcher.cache_filename_for_image

    first = name("https://host-a.example/v1/gvclass.sif")
    second = name("https://host-b.example/v2/gvclass.sif")

    assert first != second
    assert first.endswith("_gvclass.sif")
    assert first == name("https://host-a.example/v1/gvclass.sif")


def test_image_cache_name_keeps_the_readable_tag(launcher: ModuleType) -> None:
    name = launcher.cache_filename_for_image

    oras = name("oras://ghcr.io/nelli-team/gvclass:2.1.1")
    library = name("library://group/collection/gvclass:2.1.0")

    assert oras.endswith("_ghcr.io_nelli-team_gvclass_2.1.1.sif")
    assert library.endswith("_group_collection_gvclass_2.1.0.sif")


def test_image_cache_names_differ_for_one_path_under_different_schemes(
    launcher: ModuleType,
) -> None:
    name = launcher.cache_filename_for_image
    addresses = [
        "oras://ghcr.io/nelli-team/gvclass:2.1.1",
        "docker://ghcr.io/nelli-team/gvclass:2.1.1",
        "library://ghcr.io/nelli-team/gvclass:2.1.1",
        "oras://ghcr.io/nelli-team/gvclass_2.1.1",
        "oras://ghcr.io/nelli-team/gvclass:latest",
    ]

    names = [name(address) for address in addresses]

    assert len(set(names)) == len(addresses)
    assert names[0] == name(addresses[0])


def test_default_output_name_for_the_current_directory(
    launcher: ModuleType, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = tmp_path / "genomes"
    query.mkdir()
    monkeypatch.chdir(query)

    output = launcher.default_output_dir(Path())

    assert output == Path.cwd() / "genomes_results"


def test_repo_detection_ignores_an_unrelated_src_directory(
    launcher: ModuleType, tmp_path: Path
) -> None:
    (tmp_path / "src").mkdir()

    assert launcher.is_gvclass_repo(tmp_path) is False
    assert launcher.is_gvclass_repo(REPO_ROOT) is True


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


@pytest.mark.parametrize(
    ("total_threads", "threads_per_worker", "expected"),
    [(2, 8, (1, 2)), (16, 4, (4, 4)), (1, 1, (1, 1))],
)
def test_explicit_threads_per_worker_never_exceeds_the_total(
    total_threads: int, threads_per_worker: int, expected: tuple[int, int]
) -> None:
    plan = plan_workers(8, total_threads, None, threads_per_worker)

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
    scorer.labels = LabelIndex.from_taxonomy(
        {"BAC__ref": ["BAC", "", "", "OrderA", "FamA"]}
    )
    tree_neighbors = {"OG1": {name: {"BAC__ref|p1": 0.1} for name in query_ids}}

    result = scorer.collect_contig_features(
        fna,
        faa,
        scorer.scan_blast_dir(blast),
        QueryEvidence({}, {}, {}),
        "",
        "",
        tree_neighbors,
    )

    assert result["contig_count"] == 1
    assert result["suspicious_bp_fraction"] == 100.0
    assert result["cellular_coherent_bp_fraction"] == 100.0


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


def test_container_mode_without_a_query_is_an_error(
    launcher: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as raised:
        launcher.container_main(["-o", "results", "--image", "gvclass.sif"])

    assert raised.value.code == 2
    assert "query path is required" in capsys.readouterr().err


def test_container_mode_without_arguments_prints_help(
    launcher: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    assert launcher.container_main([]) == 0
    assert "usage:" in capsys.readouterr().out


def test_unrecognized_native_arguments_still_fail() -> None:
    with pytest.raises(SystemExit) as raised:
        gvclass_cli.main(["queries", "--no-such-option"])

    assert raised.value.code == 2
