"""Requested species trees must not turn inference errors into successful runs."""

import json
import subprocess
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
from click.testing import CliRunner

from src.bin import gvclass_cli, gvclass_runner
from src.core import alignment
from src.core.species_tree import (
    handoff,
    neighbor_selection,
    orchestration,
    tree_inference,
)
from src.pipeline import parallel_runner, query_processing_engine
from src.pipeline.run_config import RunConfig
from tests.pipeline_helpers import (
    fake_setup,
    make_config,
    make_fake_process,
    write_queries,
)


@pytest.mark.parametrize("mode", ["fast", "ufboot"])
def test_iqtree_has_four_hour_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    monkeypatch.setattr(alignment, "_IQTREE_MODE", mode)
    run = Mock()
    monkeypatch.setattr(subprocess, "run", run)
    tree = tmp_path / "result.treefile"
    tree.write_text("(a,b);", encoding="utf-8")

    assert alignment.run_iqtree(tmp_path / "input.faa", tree) == "(a,b);"
    assert run.call_args.kwargs["timeout"] == 14400
    assert run.call_args.kwargs["check"] is True


@pytest.mark.parametrize(
    "error",
    [
        subprocess.TimeoutExpired("iqtree", 14400),
        subprocess.CalledProcessError(1, "iqtree"),
    ],
)
def test_iqtree_subprocess_errors_propagate(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, error: Exception
) -> None:
    monkeypatch.setattr(subprocess, "run", Mock(side_effect=error))
    with pytest.raises(type(error)):
        tree_inference.infer_species_tree(
            tmp_path / "in.faa", tmp_path / "out.treefile", tree_method="iqtree"
        )


@pytest.mark.parametrize("method", ["iqtree", "veryfasttree"])
def test_empty_species_tree_is_an_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, method: str
) -> None:
    monkeypatch.setattr(tree_inference, "run_iqtree", Mock(return_value=" "))
    monkeypatch.setattr(tree_inference, "run_veryfasttree", Mock(return_value=" "))
    with pytest.raises(RuntimeError, match="produced no tree"):
        tree_inference.infer_species_tree(
            tmp_path / "in.faa", tmp_path / "out.treefile", tree_method=method
        )


@pytest.fixture
def flow_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> RunConfig:
    config = make_config(tmp_path, threads=1)
    write_queries(config, ["query"])
    monkeypatch.setattr(parallel_runner, "validate_and_setup_task", fake_setup)
    monkeypatch.setattr(parallel_runner, "check_marker_panels", lambda _database: None)
    return config


@pytest.fixture
def prepared_flow(
    flow_config: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> tuple[RunConfig, Mock]:
    monkeypatch.setattr(
        query_processing_engine,
        "_run_query_stages",
        Mock(
            return_value=(SimpleNamespace(best_code=1), {"taxonomy_majority": "NCLDV"})
        ),
    )
    monkeypatch.setattr(orchestration, "select_panel", Mock(return_value=object()))
    archive = Mock(wraps=query_processing_engine._post_process_query)
    monkeypatch.setattr(query_processing_engine, "_post_process_query", archive)
    return flow_config, archive


@pytest.mark.parametrize(
    "error",
    [
        subprocess.TimeoutExpired("iqtree", 14400),
        subprocess.CalledProcessError(1, "iqtree"),
        RuntimeError("IQ-TREE produced no tree"),
    ],
)
def test_requested_per_query_failure_marks_query_and_run_failed(
    prepared_flow: tuple[RunConfig, Mock],
    monkeypatch: pytest.MonkeyPatch,
    error: Exception,
) -> None:
    config, archive = prepared_flow
    monkeypatch.setattr(orchestration, "_run_per_query", Mock(side_effect=error))

    with pytest.raises(RuntimeError, match="failed for 1 of 1 queries"):
        parallel_runner.gvclass_flow(replace(config, species_tree=True))

    status = json.loads(
        (config.output_dir / "run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "failed"
    assert status["queries"]["query"]["status"] == "failed"
    assert str(error) in status["queries"]["query"]["error"]
    archive.assert_not_called()
    assert (config.output_dir / ".gvclass_work" / "query").is_dir()


@pytest.mark.parametrize("stage", ["setup", "combined"])
def test_requested_run_level_species_failure_is_not_suppressed(
    flow_config: RunConfig, monkeypatch: pytest.MonkeyPatch, stage: str
) -> None:
    monkeypatch.setattr(parallel_runner, "process_single_query", make_fake_process())
    target = handoff if stage == "setup" else orchestration
    name = "init_handoff" if stage == "setup" else "_run_combined"
    monkeypatch.setattr(target, name, Mock(side_effect=RuntimeError("species failure")))

    with pytest.raises(RuntimeError, match="species failure"):
        parallel_runner.gvclass_flow(replace(flow_config, species_tree_combined=True))

    status = json.loads(
        (flow_config.output_dir / "run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "failed"


def test_ineligible_panel_is_still_skipped(tmp_path: Path) -> None:
    orchestration.run_per_query_species_tree(
        "query", tmp_path, tmp_path, {"taxonomy_majority": "unknown"}, tmp_path, 1
    )


def test_empty_neighbor_tree_is_an_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = tmp_path / "query.faa"
    reference = tmp_path / "reference.faa"
    query.write_text(">q\nMA\n", encoding="utf-8")
    reference.write_text(">a\nMA\n>b\nMA\n", encoding="utf-8")
    monkeypatch.setattr(
        neighbor_selection,
        "align_sequences_pyfamsa",
        Mock(return_value=[("q", "MA"), ("a", "MA"), ("b", "MA")]),
    )
    monkeypatch.setattr(neighbor_selection, "run_veryfasttree", Mock(return_value=" "))
    with pytest.raises(RuntimeError, match="produced no tree"):
        neighbor_selection.build_species_tree_gene_tree(
            query, ["a", "b"], reference, tmp_path / "out.treefile", threads=1
        )


def test_runner_exits_nonzero_for_requested_species_failure(
    flow_config: RunConfig, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = flow_config.output_dir
    monkeypatch.setattr(parallel_runner, "process_single_query", make_fake_process())
    monkeypatch.setattr(
        orchestration,
        "_run_combined",
        Mock(side_effect=RuntimeError("species inference failed")),
    )
    result = CliRunner().invoke(
        gvclass_runner.main,
        [
            str(flow_config.query_path),
            "-o",
            str(output),
            "-d",
            str(flow_config.database_path),
            "-t",
            "1",
            "--species-tree-combined",
        ],
    )
    assert result.exit_code == 1
    assert "species inference failed" in result.output
    assert "Pipeline completed successfully" not in result.output
    assert (output / "gvclass_summary.tsv").is_file()
    assert (output / "gvclass_summary.csv").is_file()
    status = json.loads((output / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] == "failed"


def test_native_cli_preserves_runner_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        gvclass_cli,
        "run_pipeline",
        Mock(return_value=(1, b"", b"species inference failed", None)),
    )
    output = Mock()
    result = gvclass_cli.execute_pipeline_command(
        SimpleNamespace(resume=False, verbose=False),
        SimpleNamespace(n_queries=1, n_skipped=0, output_dir=tmp_path),
        tmp_path,
        ["runner"],
        {},
        output,
    )
    assert result == 1
    output.line.assert_called_once_with(
        "\nPipeline failed with exit code: 1", key="error"
    )


def test_iqtree_preserves_dotted_query_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def write_iqtree_output(command: list[str], **_kwargs: object) -> None:
        prefix = command[command.index("-pre") + 1]
        Path(f"{prefix}.treefile").write_text("(a,b);", encoding="utf-8")

    monkeypatch.setattr(subprocess, "run", write_iqtree_output)
    tree = tmp_path / "bin.12.treefile"

    assert alignment.run_iqtree(tmp_path / "input.faa", tree) == "(a,b);"
    assert tree.read_text(encoding="utf-8") == "(a,b);"


def test_healthy_species_fixture_completes_with_real_archive(
    prepared_flow: tuple[RunConfig, Mock], monkeypatch: pytest.MonkeyPatch
) -> None:
    config, archive = prepared_flow
    output = config.output_dir
    monkeypatch.setattr(orchestration, "_run_per_query", Mock())

    parallel_runner.gvclass_flow(replace(config, species_tree=True))

    status = json.loads((output / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] == "complete"
    assert status["queries"]["query"]["status"] == "complete"
    assert (output / "query.tar.gz").is_file()
    archive.assert_called_once()
    assert not (output / ".gvclass_work").exists()


def test_missing_iqtree_output_is_an_inference_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(subprocess, "run", Mock())
    tree = tmp_path / "absent.treefile"

    with pytest.raises(RuntimeError, match="IQ-TREE produced no tree"):
        tree_inference.infer_species_tree(
            tmp_path / "input.faa", tree, tree_method="iqtree"
        )

    assert not tree.exists()
