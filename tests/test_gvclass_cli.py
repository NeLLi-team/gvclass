import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import yaml

from src.bin import gvclass_cli
from src.bin.gvclass_cli import (
    BUILTIN_CONFIG_DEFAULTS,
    RESOURCE_CACHE_ENV,
    THREAD_ENV_VARS,
    ConfigError,
    ContigInput,
    build_parser,
    configure_resource_cache,
    default_output_dir,
    find_config_file,
    inspect_local_database_state,
    load_config,
    prompt_yes_no,
    resolve_completeness_mode,
    resolve_contigs_min_length,
    resolve_input_min_length,
    resolve_mode_fast,
    resolve_resource_cache_path,
    resolve_sensitive_mode,
    resolve_threads,
    resolve_validation_min_length,
    should_offer_database_update,
    thread_environment,
)
from src.pipeline.run_config import plan_workers
from src.pipeline.run_status import completed_query_names
from src.utils.database_manager import DatabaseManager
from src.utils.resource_store import ResourceStore
from tests.launcher_helpers import load_launcher

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_inspect_local_database_state_accepts_compact_resources(tmp_path):
    db_path = tmp_path / "resources"
    db_path.mkdir()
    for relative_path in DatabaseManager.REQUIRED_FILES:
        target = db_path / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("test\n")

    state = inspect_local_database_state(db_path)

    assert state["exists"] is True
    assert state["complete"] is True
    assert not (db_path / "database").exists()
    assert state["missing_files"] == []


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
    return load_config(cfg)


def test_load_config_partial_pipeline_uses_defaults(tmp_path):
    # Only one pipeline key set; the rest must fall back to defaults. Use a
    # value that differs from the default so this proves the override applies.
    config = _load(tmp_path, "pipeline:\n  threads: 12\n")
    assert config["pipeline"]["threads"] == 12
    assert config["pipeline"]["tree_method"] == "veryfasttree"
    assert config["pipeline"]["mode_fast"] is True
    assert config["pipeline"]["completeness_mode"] == "novelty-aware"
    # The whole database section must still be present.
    assert config["database"]["download_version"] == "v2.1.0"
    assert config["database"]["cache_path"] is None


def test_load_config_missing_database_uses_defaults(tmp_path):
    config = _load(tmp_path, "pipeline:\n  tree_method: iqtree\n")
    assert config["pipeline"]["tree_method"] == "iqtree"
    assert "database" in config
    assert config["database"]["download_version"] == "v2.1.0"


def test_load_config_partial_database_preserves_download_fields(tmp_path):
    config = _load(tmp_path, "database:\n  path: /custom/db\n")
    assert config["database"]["path"] == "/custom/db"
    assert config["database"]["cache_path"] is None
    # Setting only path must keep the download version needed for setup-db; the
    # source itself is built in (the archive is stored in parts).
    assert config["database"]["download_version"] == "v2.1.0"
    assert config["database"]["download_url"] is None
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
    assert config["pipeline"]["threads"] == 4
    assert config["database"]["download_version"] == "v2.1.0"
    assert config["quality"]["min_length"] == 20000


def test_builtin_defaults_match_shipped_config_yaml():
    """The shipped YAML is the single source of defaults; the CLI literal is
    only the fallback for installs without ``config/``. This keeps the literal
    from becoming a second source of truth."""
    shipped = yaml.safe_load((REPO_ROOT / "config" / "gvclass_config.yaml").read_text())

    assert shipped == BUILTIN_CONFIG_DEFAULTS


def _auto_plan(n_queries: int, threads: int):
    return plan_workers(n_queries, threads)


def test_worker_plan_never_oversubscribes():
    """Assert the invariants rather than the numbers: never hand out more threads
    than the run was given, never start more workers than there are queries."""
    for n_queries in (1, 2, 3, 4, 5, 8, 16, 50, 100, 200, 500):
        for threads in (1, 2, 3, 4, 8, 16, 32, 64, 128):
            workers, per_worker = _auto_plan(n_queries, threads)

            assert workers >= 1 and per_worker >= 1
            assert workers <= n_queries, (n_queries, threads)
            assert workers * per_worker <= threads, (n_queries, threads)


def test_worker_plan_honours_explicit_overrides():
    assert plan_workers(50, 64, 8, 3) == (8, 3)
    assert plan_workers(50, 64, 8, None) == (8, 8)
    assert plan_workers(50, 64, None, 3)[1] == 3


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


def test_parser_rejects_quality_models_option():
    with pytest.raises(SystemExit):
        build_parser().parse_args(
            ["example", "-o", "out", "--quality-models", "bundled"]
        )


def test_empty_download_keys_leave_the_built_in_database_source(tmp_path):
    config = _load(tmp_path, "database:\n  download_url:\n  download_sha256:\n")

    assert gvclass_cli.resolve_database_download_source(config) is None


def test_main_exits_2_for_an_explicit_config_that_does_not_exist(
    tmp_path, monkeypatch, capsys
):
    missing = tmp_path / "missing.yaml"
    monkeypatch.setattr(
        gvclass_cli, "check_and_setup_database", Mock(side_effect=AssertionError)
    )

    code = gvclass_cli.main([str(tmp_path), "-c", str(missing)])

    assert code == 2
    assert str(missing) in capsys.readouterr().err


def test_find_config_file_reads_the_working_directory_then_the_repository(
    tmp_path, monkeypatch
):
    repo = tmp_path / "repo"
    (repo / "config").mkdir(parents=True)
    shipped = repo / "config" / "gvclass_config.yaml"
    shipped.write_text("pipeline: {}\n")
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    local = work / "gvclass_config.yaml"

    assert find_config_file(None, repo) == shipped
    local.write_text("pipeline: {}\n")
    assert find_config_file(None, repo).resolve() == local.resolve()
    shipped.unlink()
    local.unlink()
    assert find_config_file(None, repo) is None


@pytest.mark.parametrize(
    ("text", "section"),
    [
        ("pipeline:\n", "pipeline"),
        ("database: []\n", "database"),
        ("quality: 5\n", "quality"),
    ],
)
def test_load_config_names_a_section_that_is_not_a_mapping(tmp_path, text, section):
    with pytest.raises(ConfigError, match=section):
        _load(tmp_path, text)


def test_load_config_without_pyyaml_is_an_error(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "yaml", None)

    with pytest.raises(ConfigError, match="PyYAML"):
        _load(tmp_path, "pipeline: {}\n")


@pytest.mark.parametrize(
    "option",
    [
        "--cluster-type",
        "--cluster-queue",
        "--cluster-project",
        "--cluster-walltime",
        "--apptainer",
        "--image",
    ],
)
def test_parser_rejects_removed_and_container_options(option):
    argv = ["queries", option] if option == "--apptainer" else ["queries", option, "x"]

    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(argv)

    assert raised.value.code == 2


@pytest.mark.parametrize(
    "flags",
    [
        ["-t", "0"],
        ["-t", "many"],
        ["-j", "0"],
        ["--threads-per-worker", "0"],
        ["--threads-per-worker", "-2"],
    ],
)
def test_parser_rejects_thread_counts_below_one(flags):
    with pytest.raises(SystemExit) as raised:
        build_parser().parse_args(["queries", *flags])

    assert raised.value.code == 2


def test_resolve_threads_uses_the_flag_then_the_config(tmp_path):
    config = _load(tmp_path, "pipeline:\n  threads: 6\n")

    assert resolve_threads(build_parser().parse_args(["q"]), config) == 6
    assert resolve_threads(build_parser().parse_args(["q", "-t", "3"]), config) == 3


def test_resolve_threads_rejects_a_bad_config_value(tmp_path):
    config = _load(tmp_path, "pipeline:\n  threads: 0\n")

    with pytest.raises(ValueError, match=r"pipeline\.threads"):
        resolve_threads(build_parser().parse_args(["q"]), config)


@pytest.mark.parametrize(
    "flags",
    [
        ["--extended", "--mode-fast"],
        ["--mode-fast", "--extended"],
        ["-e", "--mode-fast"],
        ["--mode-fast", "-e"],
    ],
)
def test_extended_wins_over_mode_fast_in_either_order(flags):
    args = build_parser().parse_args(["q", *flags])

    assert resolve_mode_fast(args, BUILTIN_CONFIG_DEFAULTS) is False


@pytest.mark.parametrize(
    "flags", [["--sensitive", "--no-sensitive"], ["--no-sensitive", "--sensitive"]]
)
def test_no_sensitive_wins_over_sensitive_in_either_order(flags):
    args = build_parser().parse_args(["q", *flags])

    assert resolve_sensitive_mode(args, BUILTIN_CONFIG_DEFAULTS) is False


@pytest.mark.parametrize(
    ("inherited", "expected"),
    [(None, "4"), ("2", "2"), ("4", "4"), ("8", "4"), ("many", "4"), ("0", "4")],
)
def test_thread_environment_keeps_only_a_smaller_inherited_integer(inherited, expected):
    environ = {} if inherited is None else dict.fromkeys(THREAD_ENV_VARS, inherited)

    result = thread_environment(environ, 4)

    assert {result[name] for name in THREAD_ENV_VARS} == {expected}


def test_default_output_dir_is_the_working_directory_plus_the_input_name(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    folder = tmp_path / "bins"
    folder.mkdir()
    assembly = tmp_path / "assembly.fna"
    assembly.write_text(">c\nACGT\n", encoding="utf-8")

    assert default_output_dir(folder) == Path.cwd() / "bins_results"
    assert default_output_dir(assembly) == Path.cwd() / "assembly_results"


def test_native_and_container_runs_choose_the_same_default_output_dir(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    folder = tmp_path / "bins"
    folder.mkdir()
    assembly = tmp_path / "assembly.fna"
    assembly.write_text(">c\nACGT\n", encoding="utf-8")
    launcher = load_launcher()

    assert launcher.default_output_dir(folder) == default_output_dir(folder)
    assert launcher.default_output_dir(assembly) == default_output_dir(assembly)


def _record(query: str, status: str, archive: dict | None) -> str:
    return json.dumps(
        {
            "query": query,
            "status": status,
            "error": None,
            "input_name": f"{query}.faa",
            "input_sha256": "0" * 64,
            "settings": {},
            "database_version": "v2.1.0",
            "software_version": "v2.1.1",
            "finished_at": "2026-11-09T00:00:00+00:00",
            "archive": archive,
            "summary": None,
        }
    )


def test_completed_query_names_lists_complete_records_with_intact_archives(tmp_path):
    def archive(name: str, size: int) -> dict:
        return {"filename": f"{name}.tar.gz", "bytes": size, "sha256": "x"}

    (tmp_path / "a.tar.gz").write_bytes(b"12345")
    (tmp_path / "b.tar.gz").write_bytes(b"123")
    lines = [
        _record("a", "complete", archive("a", 5)),
        _record("b", "complete", archive("b", 99)),
        _record("c", "failed", None),
        _record("d", "complete", archive("d", 1)),
        _record("a", "failed", None),
        _record("a", "complete", archive("a", 5)),
    ]
    (tmp_path / "results.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")

    assert completed_query_names(tmp_path) == {"a"}


def test_importing_the_cli_loads_no_scientific_library():
    code = (
        "import sys; sys.path.insert(0, sys.argv[1]); import src.bin.gvclass_cli; "
        "print(sorted(m for m in ('Bio', 'numpy', 'pandas') if m in sys.modules))"
    )

    result = subprocess.run(
        [sys.executable, "-c", code, str(REPO_ROOT)],
        capture_output=True,
        text=True,
        check=True,
    )

    assert result.stdout.strip() == "[]"


def test_resolve_completeness_mode_honors_explicit_legacy() -> None:
    args = SimpleNamespace(completeness_mode=None)
    config = {"pipeline": {"completeness_mode": "legacy"}}
    assert resolve_completeness_mode(args, config) == "legacy"


def test_resolve_completeness_mode_cli_flag_wins_over_config() -> None:
    args = SimpleNamespace(completeness_mode="legacy")
    config = {"pipeline": {"completeness_mode": "novelty-aware"}}
    assert resolve_completeness_mode(args, config) == "legacy"


def test_unrecognized_native_arguments_still_fail() -> None:
    with pytest.raises(SystemExit) as raised:
        gvclass_cli.main(["queries", "--no-such-option"])

    assert raised.value.code == 2
