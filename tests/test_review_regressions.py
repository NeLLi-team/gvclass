"""Regression tests for defects found by the independent review of the 2.1.0 launcher.

The launcher exists twice, as the ``gvclass`` script and as
``src.bin.gvclass_apptainer``; each behaviour is checked on both copies.
"""

from __future__ import annotations

import runpy
from pathlib import Path
from typing import Any

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.bin import gvclass_apptainer
from src.core.contamination_scoring import ContaminationScorer
from src.core.reformat import reformat_sequences
from src.pipeline.parallel_runner import _resolve_worker_distribution
from src.pipeline.query_processing_engine import _initialize_query
from src.utils.error_handling import ProcessingError

REPO_ROOT = Path(__file__).resolve().parents[1]


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
