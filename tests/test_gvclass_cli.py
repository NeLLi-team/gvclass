import os
from pathlib import Path

import pytest

from src.bin.gvclass_cli import (
    CliOutput,
    ContigInput,
    PipelineContext,
    RESOURCE_CACHE_ENV,
    WorkerPlan,
    build_pipeline_command,
    configure_resource_cache,
    inspect_local_database_state,
    load_config,
    prompt_yes_no,
    resolve_contigs_min_length,
    resolve_input_min_length,
    resolve_resource_cache_path,
    resolve_validation_min_length,
    should_offer_database_update,
)
from src.utils.database_manager import DatabaseManager
from src.utils.resource_store import ResourceStore


def test_inspect_local_database_state_requires_database_directory(tmp_path):
    db_path = tmp_path / "resources"
    db_path.mkdir()
    for relative_path in DatabaseManager.REQUIRED_FILES:
        target = db_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test\n")

    state = inspect_local_database_state(db_path)

    assert state["exists"] is True
    assert state["complete"] is False
    assert "database" in state["missing_files"]


def test_prompt_yes_no_accepts_enter_as_default_yes(monkeypatch):
    monkeypatch.setattr("builtins.input", lambda _: "")

    assert prompt_yes_no("Update database?", default=True) is True


def test_should_offer_database_update_uses_version_comparison():
    assert should_offer_database_update("v1.2.2", {"version": "v1.4.0"}) is True
    assert should_offer_database_update("v1.4.0", {"version": "v1.4.0"}) is False
    assert should_offer_database_update("v2.0.0", {"version": "v1.7.1"}) is False


def _load(tmp_path: Path, text: str):
    cfg = tmp_path / "gvclass_config.yaml"
    cfg.write_text(text)
    return load_config(str(cfg), tmp_path, CliOutput(plain_output=True))


def test_load_config_partial_pipeline_uses_defaults(tmp_path):
    # Only one pipeline key set; the rest must fall back to defaults.
    config = _load(tmp_path, "pipeline:\n  threads: 4\n")
    assert config["pipeline"]["threads"] == 4
    assert config["pipeline"]["tree_method"] == "veryfasttree"
    assert config["pipeline"]["mode_fast"] is True
    assert config["pipeline"]["completeness_mode"] == "novelty-aware"
    # The whole database section must still be present.
    assert config["database"]["download_version"] == "v2.0.0"
    assert config["database"]["cache_path"] is None


def test_load_config_missing_database_uses_defaults(tmp_path):
    config = _load(tmp_path, "pipeline:\n  tree_method: iqtree\n")
    assert config["pipeline"]["tree_method"] == "iqtree"
    assert "database" in config
    assert config["database"]["download_url"].startswith("https://")


def test_load_config_partial_database_preserves_download_fields(tmp_path):
    config = _load(tmp_path, "database:\n  path: /custom/db\n")
    assert config["database"]["path"] == "/custom/db"
    assert config["database"]["cache_path"] is None
    # Setting only path must not drop the download fields needed for setup-db.
    assert config["database"]["download_url"].startswith("https://")
    assert config["database"]["download_sha256"]
    assert config["pipeline"]["sensitive_mode"] is True


def test_resolve_resource_cache_path_relative_to_database(tmp_path, monkeypatch):
    monkeypatch.delenv(RESOURCE_CACHE_ENV, raising=False)
    db_path = tmp_path / "resources"
    config = {"database": {"cache_path": ".gvclass_cache"}}

    assert (
        resolve_resource_cache_path(config, db_path)
        == (db_path / ".gvclass_cache").resolve()
    )


def test_configure_resource_cache_respects_env_override(tmp_path, monkeypatch):
    monkeypatch.setenv(RESOURCE_CACHE_ENV, "/env/cache")
    config = {"database": {"cache_path": "configured_cache"}}

    configure_resource_cache(config, tmp_path / "resources")

    assert os.environ[RESOURCE_CACHE_ENV] == "/env/cache"


def test_configure_resource_cache_sets_configured_path_for_resource_store(
    tmp_path, monkeypatch
):
    monkeypatch.delenv(RESOURCE_CACHE_ENV, raising=False)
    db_path = tmp_path / "resources"
    config = {"database": {"cache_path": "configured_cache"}}
    expected = (db_path / "configured_cache").resolve()

    configure_resource_cache(config, db_path)

    assert os.environ[RESOURCE_CACHE_ENV] == str(expected)
    assert ResourceStore(db_path)._cache_root() == expected


def test_load_config_empty_file_returns_defaults(tmp_path):
    config = _load(tmp_path, "# only a comment\n")
    assert config["pipeline"]["threads"] == 16
    assert config["database"]["download_version"] == "v2.0.0"
    assert config["quality"]["min_length"] == 20000


def test_resolve_input_min_length_uses_cli_over_config(tmp_path):
    import argparse

    config = _load(tmp_path, "quality:\n  min_length: 30000\n")
    args = argparse.Namespace(min_length=12000)

    assert resolve_input_min_length(args, config) == 12000


def test_resolve_contigs_min_length_uses_cli_over_config(tmp_path):
    import argparse

    config = _load(tmp_path, "pipeline:\n  contigs_min_length: 30000\n")
    args = argparse.Namespace(contigs_min_length=5000)

    assert resolve_contigs_min_length(args, config) == 5000


def test_resolve_min_length_rejects_negative_cli_values(tmp_path):
    import argparse

    config = _load(tmp_path, "quality:\n  min_length: 20000\n")

    with pytest.raises(ValueError, match="--min-length"):
        resolve_input_min_length(argparse.Namespace(min_length=-1), config)

    with pytest.raises(ValueError, match="--contigs-min-length"):
        resolve_contigs_min_length(
            argparse.Namespace(contigs_min_length=-1),
            config,
        )


def test_resolve_validation_min_length_uses_contig_floor_in_contigs_mode():
    assert (
        resolve_validation_min_length(
            ContigInput(enabled=True, input_file=Path("input.fna")),
            input_min_length=30000,
            contigs_min_length=5000,
        )
        == 5000
    )
    assert (
        resolve_validation_min_length(
            ContigInput(enabled=False),
            input_min_length=30000,
            contigs_min_length=5000,
        )
        == 30000
    )


def test_build_pipeline_command_forwards_effective_min_length(tmp_path):
    import argparse

    args = argparse.Namespace(
        cluster_type="local",
        cluster_queue=None,
        cluster_project=None,
        cluster_walltime="04:00:00",
        verbose=False,
        resume=False,
        allow_short=False,
        species_tree=False,
        species_tree_combined=False,
        species_tree_trim="witchi",
    )
    context = PipelineContext(
        query_dir=tmp_path / "queries",
        active_query_dir=tmp_path / "split_contigs",
        output_dir=tmp_path / "out",
        database=tmp_path / "resources",
        threads=8,
        tree_method="veryfasttree",
        iqtree_mode="fast",
        mode_fast=True,
        completeness_mode="novelty-aware",
        sensitive_mode=True,
        contigs=ContigInput(enabled=True, input_file=tmp_path / "input.fna"),
        input_min_length=20000,
        contigs_min_length=5000,
        validation_min_length=5000,
        temp_contigs_dir=tmp_path / "split_contigs",
        n_queries=2,
        n_input_files=1,
        n_skipped=0,
        workers=WorkerPlan(workers=2, threads_per_worker=4),
    )

    cmd = build_pipeline_command(
        args=args,
        context=context,
        python_exe="python",
        pixi_cmd=None,
        use_pixi_run=False,
    )

    assert cmd[cmd.index("--min-length") + 1] == "5000"
