#!/usr/bin/env python
"""GVClass command line: resolve the options into a run configuration and run it.

``build_parser`` is the one definition of the pipeline options. ``main`` resolves them
with the configuration file into a ``RunConfig``, writes it to ``run_config.json`` in
the output directory and starts ``python -m src.bin.gvclass_runner`` on that file.
"""

from __future__ import annotations

import argparse
import copy
import os
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import deque
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager, suppress
from dataclasses import dataclass, replace
from pathlib import Path
from types import FrameType
from typing import Any, Optional

from src.__version__ import __version__ as _GVCLASS_VERSION
from src.bin.cli_display import print_banner, print_configuration
from src.bin.progress_monitor import ResourceMonitor
from src.pipeline.run_config import RunConfig, plan_workers

SOFTWARE_VERSION = f"v{_GVCLASS_VERSION}"
PLAIN_OUTPUT_ENV = "GVCLASS_PLAIN_OUTPUT"
DATABASE_PATH_ENV = "GVCLASS_DB"
RESOURCE_CACHE_ENV = "GVCLASS_RESOURCE_CACHE"
DEFAULT_INPUT_MIN_LENGTH = 20000
DEFAULT_CONTIGS_MIN_LENGTH = 10000
CONFIG_FILENAME = "gvclass_config.yaml"
RUN_CONFIG_FILENAME = "run_config.json"
#: The numeric libraries that start their own threads.
THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)
#: A failed run prints this many of the last lines of the runner output.
RUNNER_TAIL_LINES = 200
MONITOR_THREAD_NAME = "gvclass-memory-monitor"
#: Seconds that the runner and its descendants get to stop after SIGTERM before they
#: are killed.
STOP_GRACE_S = 10.0
STOP_POLL_S = 0.05
MAX_ISSUES_SHOWN = 5


class CliOutput:
    """Console output helper with optional plain-text mode."""

    ICONS = {
        "success": "✅",
        "warning": "⚠️",
        "error": "❌",
        "database": "📦",
        "download": "📥",
        "contigs": "📂",
        "validate": "🔍",
        "summary": "📊",
        "runtime": "⏱️",
        "memory": "💾",
    }

    def __init__(self, plain_output: bool = False):
        self.plain_output = plain_output

    def _prefix(self, key: str) -> str:
        if self.plain_output:
            return ""
        icon = self.ICONS.get(key, "")
        return f"{icon} " if icon else ""

    def line(self, message: str = "", key: Optional[str] = None, file=None):
        prefix = self._prefix(key) if key else ""
        print(f"{prefix}{message}", file=file)


class ConfigError(ValueError):
    """A configuration file that GVClass cannot use."""


@dataclass
class ContigInput:
    """Resolved contig splitting input."""

    enabled: bool
    input_file: Optional[Path] = None
    input_dir: Optional[Path] = None


@dataclass(frozen=True, slots=True)
class PipelineContext:
    """One resolved invocation: the run configuration and what the banner shows.

    ``config.query_path`` is the folder of the queries that the runner reads; with
    ``--contigs`` it is the temporary folder of the split contigs.
    """

    config: RunConfig
    config_file: Path | None
    contigs: ContigInput
    n_queries: int
    n_input_files: int
    n_skipped: int
    workers: int
    threads_per_worker: int

    @property
    def n_remaining(self) -> int:
        """Return the queries that the runner will process."""
        return self.n_queries - self.n_skipped


@dataclass(frozen=True, slots=True)
class RunnerOutcome:
    """How the runner process ended.

    ``tail`` holds the last lines of its output that were not progress events; it is
    empty in verbose mode, where the output goes straight to the terminal.
    """

    return_code: int
    tail: tuple[str, ...]
    interrupted: bool
    peak_memory_mb: float


def positive_int(text: str) -> int:
    """Parse ``text`` as an integer of at least 1, for argparse."""
    try:
        value = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected an integer, got {text!r}") from None
    if value < 1:
        raise argparse.ArgumentTypeError(f"expected at least 1, got {value}")
    return value


def _add_core_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "query_dir", nargs="?", help="Query directory containing .fna or .faa files"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output directory (default: <name>_results in the current directory)",
    )
    parser.add_argument(
        "-c",
        "--config",
        help=(
            f"Configuration file (default: ./{CONFIG_FILENAME}, then "
            f"config/{CONFIG_FILENAME} in the repository, then built-in defaults)"
        ),
    )
    parser.add_argument(
        "-d",
        "--database",
        help=f"Database path (overrides {DATABASE_PATH_ENV} and config)",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        help="Number of threads (overrides config)",
    )
    parser.add_argument(
        "-j", "--max-workers", type=positive_int, help="Maximum parallel workers"
    )
    parser.add_argument(
        "--threads-per-worker", type=positive_int, help="Threads per worker"
    )
    parser.add_argument(
        "--tree-method",
        choices=["veryfasttree", "iqtree", "fasttree"],
        help="Tree building method (overrides config); default veryfasttree",
    )
    parser.add_argument(
        "--iqtree-mode",
        choices=["fast", "ufboot"],
        help="IQ-TREE species-tree search mode (overrides config)",
    )
    parser.add_argument(
        "--resource-cache-dir",
        help="Directory for Parquet materialization cache",
    )


def _add_mode_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--mode-fast", action="store_true", help="Enable fast mode (overrides config)"
    )
    parser.add_argument(
        "-e",
        "--extended",
        action="store_true",
        help="Extended mode: build trees for all markers (slower but more comprehensive)",
    )
    parser.add_argument(
        "--sensitive",
        action="store_true",
        help="Force sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
    )
    parser.add_argument(
        "--no-sensitive",
        action="store_true",
        help="Disable sensitive HMM mode: use model GA cutoffs when available",
    )
    parser.add_argument(
        "--completeness-mode",
        choices=["legacy", "novelty-aware"],
        help="Completeness estimator to surface as the primary estimate",
    )
    parser.add_argument(
        "--plain-output",
        action="store_true",
        help="Disable emojis and ANSI colors in CLI output",
    )
    parser.add_argument(
        "-C",
        "--contigs",
        action="store_true",
        help="Split FNA file(s) into individual contigs for analysis (accepts file or directory)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        help=(
            "Minimum total nucleotide length (bp) for each input FNA. "
            "Overrides quality.min_length. Use 0 to disable the length floor."
        ),
    )
    parser.add_argument(
        "--contigs-min-length",
        type=int,
        help=(
            "Minimum contig length (bp) when splitting inputs with --contigs. "
            "Overrides pipeline.contigs_min_length."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from previous run, skipping completed queries",
    )
    parser.add_argument(
        "--allow-short",
        action="store_true",
        help=(
            "Bypass the nucleotide length floor. Use for exploratory runs "
            "on short inputs."
        ),
    )
    parser.add_argument(
        "--no-screen",
        action="store_true",
        help=(
            "Skip the contig screen that adds the screen_* summary columns "
            "(overrides pipeline.contig_screen)."
        ),
    )
    parser.add_argument(
        "--species-tree",
        action="store_true",
        help=(
            "Build one supermatrix species tree per query (NCLDV GVOG8 by "
            "default; PPV/MIRUS routed by taxonomy) and assign tree-placement "
            "taxonomy (writes out/species_tree/<query>/ and summary columns)."
        ),
    )
    parser.add_argument(
        "--species-tree-combined",
        action="store_true",
        help=(
            "Implies --species-tree; additionally build one combined species "
            "tree over all queries at once (writes out/species_tree/combined.*)."
        ),
    )
    parser.add_argument(
        "--species-tree-trim",
        choices=["witchi", "pytrimal", "none"],
        default="witchi",
        help=(
            "Supermatrix column trimming for --species-tree: witchi (rigorous, "
            "default), pytrimal (fast), or none (unpruned)."
        ),
    )


def _add_misc_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "--version", action="store_true", help="Show version information"
    )


def build_parser() -> argparse.ArgumentParser:
    """Return the parser that defines every pipeline option."""
    parser = argparse.ArgumentParser(
        description="Run GVClass pipeline with clean output",
        epilog=(
            "Container mode is handled by the gvclass launcher before this parser "
            "runs: --apptainer runs the Apptainer image, and --image IMAGE selects "
            "an image address or .sif path (it implies --apptainer)."
        ),
    )
    _add_core_arguments(parser)
    _add_mode_arguments(parser)
    _add_misc_arguments(parser)
    return parser


def resolve_plain_output(args) -> bool:
    if args.plain_output:
        return True
    env_value = os.environ.get(PLAIN_OUTPUT_ENV, "").strip().lower()
    return env_value in {"1", "true", "yes", "on"}


#: Fallback for installs without ``config/gvclass_config.yaml``, which is the
#: single source of these defaults. A test pins the two together.
BUILTIN_CONFIG_DEFAULTS = {
    "database": {
        "path": "resources",
        "cache_path": None,
        "download_url": None,
        "download_version": "v2.1.0",
        "download_sha256": None,
    },
    "pipeline": {
        "tree_method": "veryfasttree",
        "iqtree_mode": "fast",
        "mode_fast": True,
        "completeness_mode": "novelty-aware",
        "sensitive_mode": True,
        "contig_screen": True,
        "contigs_min_length": DEFAULT_CONTIGS_MIN_LENGTH,
        "threads": 4,
    },
    "quality": {
        "min_length": DEFAULT_INPUT_MIN_LENGTH,
    },
}


def find_config_file(config_file: str | None, repo_dir: Path) -> Path | None:
    """Return the configuration file to read, or ``None`` for the built-in defaults.

    An explicit ``config_file`` must exist. Without one, the file is
    ``gvclass_config.yaml`` in the working directory, else ``config/`` of the
    repository.

    Raises:
        ConfigError: If ``config_file`` does not exist.
    """
    if config_file is not None:
        path = Path(config_file).expanduser()
        if not path.is_file():
            raise ConfigError(f"configuration file not found: {path}")
        return path
    for candidate in (
        Path.cwd() / CONFIG_FILENAME,
        repo_dir / "config" / CONFIG_FILENAME,
    ):
        if candidate.is_file():
            return candidate
    return None


def load_config(config_path: Path | None) -> dict[str, Any]:
    """Read ``config_path`` over the built-in defaults; ``None`` gives the defaults.

    Raises:
        ConfigError: If PyYAML is missing, or the file is not a mapping of sections,
            or a known section is empty or not a mapping.
    """
    defaults = copy.deepcopy(BUILTIN_CONFIG_DEFAULTS)
    if config_path is None:
        return defaults
    try:
        import yaml  # noqa: PLC0415  # --help and --version work without PyYAML.
    except ModuleNotFoundError as exc:
        raise ConfigError(
            f"PyYAML is needed to read {config_path}; "
            "install it in the GVClass environment"
        ) from exc
    with config_path.open(encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if loaded is None:
        return defaults
    if not isinstance(loaded, dict):
        raise ConfigError(
            f"{config_path}: expected sections of keys, got {type(loaded).__name__}"
        )
    for section in BUILTIN_CONFIG_DEFAULTS:
        if section in loaded and not isinstance(loaded[section], dict):
            raise ConfigError(
                f"{config_path}: section '{section}' must hold keys, "
                f"got {loaded[section]!r}"
            )
    return _deep_merge_config(defaults, loaded)


def _deep_merge_config(base: dict, overlay: dict) -> dict:
    """Recursively merge ``overlay`` onto a copy of ``base``.

    Nested dicts merge key-by-key so a config that omits a whole section (e.g.
    ``pipeline``) or only sets one key still inherits every default. Scalars and
    lists in ``overlay`` overwrite the corresponding ``base`` value.
    """
    merged = dict(base)
    for key, value in overlay.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge_config(merged[key], value)
        else:
            merged[key] = value
    return merged


def resolve_database_setting(database_override: Optional[str], config) -> str:
    if database_override:
        return database_override

    env_database = os.environ.get(DATABASE_PATH_ENV, "").strip()
    if env_database:
        return env_database

    return config["database"]["path"]


def resolve_database_path(args, config, repo_dir: Path) -> Path:
    database = Path(resolve_database_setting(args.database, config)).expanduser()
    if not database.is_absolute():
        database = repo_dir / database
    return database.resolve()


def resolve_resource_cache_path(config, database: Path) -> Optional[Path]:
    """Return configured resource cache path, or None for ResourceStore default."""
    env_cache = os.environ.get(RESOURCE_CACHE_ENV, "").strip()
    if env_cache:
        return None

    cache_setting = config.get("database", {}).get("cache_path")
    if cache_setting is None:
        return None

    cache_setting = str(cache_setting).strip()
    if not cache_setting:
        return None

    cache_path = Path(cache_setting).expanduser()
    if not cache_path.is_absolute():
        cache_path = database / cache_path
    return cache_path.resolve()


def configure_resource_cache(config, database: Path) -> None:
    """Expose config-file cache setting to ResourceStore via its env hook."""
    cache_path = resolve_resource_cache_path(config, database)
    if cache_path is not None:
        os.environ[RESOURCE_CACHE_ENV] = str(cache_path)


def resolve_database_download_source(config) -> Optional[dict]:
    database_config = config.get("database", {})
    url = str(database_config.get("download_url") or "").strip()
    if not url:
        return None

    source = {"url": url}
    version = str(database_config.get("download_version", "")).strip()
    if version:
        source["version"] = version

    sha256 = database_config.get("download_sha256")
    if sha256 is not None:
        sha256 = str(sha256).strip()
        if sha256:
            source["sha256"] = sha256

    return source


def resolve_tree_method(args, config) -> str:
    return args.tree_method or config["pipeline"]["tree_method"]


def resolve_iqtree_mode(args, config) -> str:
    return args.iqtree_mode or config["pipeline"]["iqtree_mode"]


def resolve_mode_fast(args, config) -> bool:
    if args.extended:
        return False
    if args.mode_fast:
        return True
    return config["pipeline"]["mode_fast"]


def resolve_sensitive_mode(args, config) -> bool:
    if getattr(args, "no_sensitive", False):
        return False
    if getattr(args, "sensitive", False):
        return True
    return bool(config["pipeline"]["sensitive_mode"])


def resolve_contig_screen(args: argparse.Namespace, config: dict[str, Any]) -> bool:
    """Return True if the run screens its bins: ``--no-screen`` wins over the config."""
    if getattr(args, "no_screen", False):
        return False
    return bool(config["pipeline"]["contig_screen"])


def resolve_completeness_mode(args, config) -> str:
    return args.completeness_mode or config["pipeline"]["completeness_mode"]


def _resolve_int(raw_value: object, setting_name: str, *, minimum: int) -> int:
    try:
        value = int(raw_value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Invalid {setting_name} value: {raw_value!r} (must be an integer)"
        ) from exc
    if value < minimum:
        raise ValueError(
            f"Invalid {setting_name} value: {value} (must be >= {minimum})"
        )
    return value


def resolve_threads(args: argparse.Namespace, config: Mapping[str, Any]) -> int:
    """Return ``-t``, else ``pipeline.threads`` of the configuration."""
    if args.threads is not None:
        return args.threads
    return _resolve_int(config["pipeline"]["threads"], "pipeline.threads", minimum=1)


def resolve_input_min_length(args, config) -> int:
    if args.min_length is not None:
        return _resolve_int(args.min_length, "--min-length", minimum=0)
    return _resolve_int(
        config["quality"]["min_length"], "quality.min_length", minimum=0
    )


def resolve_contigs_min_length(args, config) -> int:
    if args.contigs_min_length is not None:
        return _resolve_int(args.contigs_min_length, "--contigs-min-length", minimum=0)
    return _resolve_int(
        config["pipeline"]["contigs_min_length"],
        "pipeline.contigs_min_length",
        minimum=0,
    )


def resolve_validation_min_length(
    contigs: ContigInput,
    input_min_length: int,
    contigs_min_length: int,
) -> int:
    return contigs_min_length if contigs.enabled else input_min_length


def default_output_dir(input_path: Path) -> Path:
    """Return ``<name>_results`` in the working directory.

    ``<name>`` is the name of the input directory, or the file name without its
    suffix for a single file. The launcher gives the same answer in container mode.
    """
    resolved = input_path.resolve()
    name = resolved.stem if resolved.is_file() else resolved.name
    return Path.cwd() / f"{name}_results"


def resolve_output_dir(args: argparse.Namespace, input_path: Path) -> Path:
    """Return ``-o`` as an absolute path, else the default output directory."""
    if args.output_dir:
        return Path(args.output_dir).resolve()
    return default_output_dir(input_path)


def resolve_contigs_input(query_dir_abs: Path, enabled: bool) -> ContigInput:
    if not enabled:
        return ContigInput(enabled=False)

    if query_dir_abs.is_file():
        if query_dir_abs.suffix.lower() not in [".fna", ".fa", ".fasta", ".fas"]:
            raise ValueError(
                f"--contigs requires FNA file(s) (.fna, .fa, .fasta, .fas): {query_dir_abs}"
            )
        return ContigInput(enabled=True, input_file=query_dir_abs)

    if query_dir_abs.is_dir():
        return ContigInput(enabled=True, input_dir=query_dir_abs)

    raise ValueError(f"--contigs input not found: {query_dir_abs}")


def split_contig_inputs(
    contigs: ContigInput, output: CliOutput, min_length: int, workdir: Path
) -> int:
    """Write the contigs of the input into ``workdir``; return the input file count."""
    from src.utils.contig_splitter import split_contigs, split_contigs_from_directory

    if contigs.input_file:
        prefix = contigs.input_file.stem
        output.line(
            f"Splitting contigs from: {contigs.input_file.name} (min length: {min_length} bp)",
            key="contigs",
        )
        _, n_contigs = split_contigs(
            contigs.input_file,
            prefix=prefix,
            min_length=min_length,
            output_dir=workdir,
        )
        output.line(f"Found {n_contigs} contigs >= {min_length} bp", key="success")
        return 1

    output.line(
        f"Splitting contigs from directory: {contigs.input_dir} (min length: {min_length} bp)",
        key="contigs",
    )
    _, n_contigs, n_files = split_contigs_from_directory(
        contigs.input_dir,
        min_length=min_length,
        output_dir=workdir,
    )
    output.line(
        f"Found {n_contigs} contigs >= {min_length} bp from {n_files} files",
        key="success",
    )
    return n_files
def inspect_local_database_state(db_path: Path) -> dict:
    from src.utils.database_manager import DatabaseManager

    if not db_path.exists():
        return {
            "exists": False,
            "complete": False,
            "version": "not installed",
            "missing_files": list(DatabaseManager.REQUIRED_FILES),
        }

    missing_files = DatabaseManager._check_missing_files(db_path)
    version = DatabaseManager.get_database_version(db_path)
    return {
        "exists": True,
        "complete": not missing_files,
        "version": version,
        "has_release_version": DatabaseManager.parse_version(version) is not None,
        "missing_files": missing_files,
    }


def prompt_yes_no(prompt: str, default: bool = True) -> bool:
    suffix = "[Y/n]" if default else "[y/N]"
    while True:
        try:
            response = input(f"{prompt} {suffix}: ").strip().lower()
        except EOFError:
            return default

        if not response:
            return default
        if response in {"y", "yes"}:
            return True
        if response in {"n", "no"}:
            return False

        print("Please answer yes or no.")


def is_interactive_terminal() -> bool:
    """Tell whether both stdin and stdout are terminals, so a prompt can be shown."""
    return sys.stdin.isatty() and sys.stdout.isatty()


def should_offer_database_update(
    local_version: str, latest_source: Optional[dict]
) -> bool:
    """Tell whether the pinned source is newer than a comparable installed version.

    A database without a release version number is a custom one; it is never
    offered an update.
    """
    from src.utils.database_manager import DatabaseManager

    if not latest_source:
        return False
    if DatabaseManager.parse_version(local_version) is None:
        return False
    return DatabaseManager.is_newer_version(latest_source["version"], local_version)


def confirm_database_update(db_state: dict, source: dict, output: CliOutput) -> bool:
    """Apply the update policy to a complete database; True means replace it.

    A complete database is replaced only after a "yes" at a terminal prompt
    (default No). Without a terminal, one warning names the update command and
    the run goes on with the installed database.
    """
    installed, pinned = db_state["version"], source["version"]
    if not db_state["has_release_version"]:
        output.line(
            f"Installed database has no release version ({installed}); "
            "treating it as a custom database and not checking for updates.",
            key="database",
        )
        return False
    if not should_offer_database_update(installed, source):
        return False
    if is_interactive_terminal():
        return prompt_yes_no(
            f"Update database from {installed} to {pinned}?", default=False
        )
    output.line(
        f"Installed database {installed} is older than the pinned {pinned}; "
        "continuing with the installed database. "
        "To update, run: pixi run setup-db --update",
        key="warning",
        file=sys.stderr,
    )
    return False


def check_and_setup_database(
    db_path: Path, config: dict, output: CliOutput, *, update: bool = False
) -> bool:
    """Install a missing or incomplete database; replace a complete one on consent.

    ``update`` is the explicit consent of ``setup-db --update``.
    """
    from src.utils.database_manager import DatabaseManager

    preferred_source = resolve_database_download_source(config)
    db_state = inspect_local_database_state(db_path)
    source = DatabaseManager.get_effective_database_source(preferred_source)

    if db_state["exists"] and db_state["complete"]:
        output.line(f"Database found at: {db_path}", key="success")
        output.line(
            f"Installed database version: {db_state['version']}", key="database"
        )
        output.line(f"Configured database version: {source['version']}", key="download")
        if not update and not confirm_database_update(db_state, source, output):
            return True
        output.line(
            f"Updating database from {db_state['version']} to {source['version']}...",
            key="download",
        )
        return install_database(
            DatabaseManager.setup_database, db_path, source, output, force=True
        )

    if db_state["exists"]:
        output.line(
            f"Database exists but missing files: {db_state['missing_files']}",
            key="warning",
        )
        output.line(
            f"Installed database version: {db_state['version']}", key="database"
        )
    else:
        output.line(f"\nDatabase not found at {db_path}", key="database")

    latest_source = DatabaseManager.get_latest_database_source(preferred_source)
    if latest_source:
        output.line(
            f"Latest database version on Zenodo: {latest_source['version']}",
            key="download",
        )
    output.line("Setting up database...", key="download")
    return install_database(
        DatabaseManager.setup_database, db_path, source, output, force=False
    )


def install_database(
    setup_database: Callable[..., Path],
    db_path: Path,
    source: dict,
    output: CliOutput,
    *,
    force: bool,
) -> bool:
    """Run the installer and report the outcome; False means it failed."""
    try:
        setup_database(str(db_path), preferred_source=source, force=force)
    except Exception as exc:
        output.line(f"\nFailed to setup database: {exc}", key="error")
        output.line("You can try manually with: pixi run setup-db")
        return False
    output.line("Database setup complete!", key="success")
    return True

def discover_queries(query_dir: Path) -> list[Path]:
    """Return the query files of ``query_dir``, a folder that must hold some."""
    from src.utils.input_validation import (  # noqa: PLC0415  # Biopython loads only for a run.
        discover_query_files,
    )

    if not query_dir.is_dir():
        raise ValueError(
            f"Query path is not a directory: {query_dir} "
            "(use --contigs for a single FASTA file)"
        )
    query_files = discover_query_files(query_dir)
    if not query_files:
        raise ValueError(f"No .fna or .faa files found in {query_dir}")
    return query_files


def count_completed_queries(output_dir: Path, query_files: Sequence[Path]) -> int:
    """Count the queries that ``--resume`` would skip, from ``results.jsonl``.

    The count needs the archive in place. The runner also compares the input, the
    settings and the versions, so it can skip fewer queries than this preview shows.
    """
    from src.pipeline.run_status import (  # noqa: PLC0415  # Pipeline code loads only for a run.
        completed_query_names,
    )

    return len(completed_query_names(output_dir) & {path.stem for path in query_files})


def resolve_run_config(
    args: argparse.Namespace,
    config: dict[str, Any],
    repo_dir: Path,
    contigs: ContigInput,
) -> RunConfig:
    """Resolve the options and the configuration into one run configuration.

    ``query_path`` is the input that the user named. With ``--contigs``, ``run_cli``
    replaces it by the folder of the split contigs.
    """
    input_path = Path(args.query_dir).resolve()
    contigs_min_length = resolve_contigs_min_length(args, config)
    return RunConfig(
        query_path=input_path,
        output_dir=resolve_output_dir(args, input_path),
        database_path=resolve_database_path(args, config, repo_dir),
        threads=resolve_threads(args, config),
        max_workers=args.max_workers,
        threads_per_worker=args.threads_per_worker,
        tree_method=resolve_tree_method(args, config),
        iqtree_mode=resolve_iqtree_mode(args, config),
        mode_fast=resolve_mode_fast(args, config),
        completeness_mode=resolve_completeness_mode(args, config),
        sensitive_mode=resolve_sensitive_mode(args, config),
        resume=args.resume,
        allow_short=args.allow_short,
        min_sequence_length=resolve_validation_min_length(
            contigs, resolve_input_min_length(args, config), contigs_min_length
        ),
        species_tree=args.species_tree,
        species_tree_combined=args.species_tree_combined,
        species_tree_trim=args.species_tree_trim,
        contig_screen=resolve_contig_screen(args, config),
        contigs_mode=contigs.enabled,
        contigs_min_length=contigs_min_length,
    )


def apply_resource_cache(
    args: argparse.Namespace, config: dict[str, Any], database: Path
) -> None:
    """Pass ``--resource-cache-dir``, else the configured cache, to the runner."""
    if args.resource_cache_dir:
        os.environ[RESOURCE_CACHE_ENV] = str(
            Path(args.resource_cache_dir).expanduser().resolve()
        )
    else:
        configure_resource_cache(config, database)


def plan_run(
    args: argparse.Namespace,
    run_config: RunConfig,
    query_dir: Path,
    *,
    config_file: Path | None,
    contigs: ContigInput,
    n_input_files: int,
) -> PipelineContext:
    """Count the queries in ``query_dir`` and plan the workers for the banner."""
    query_files = discover_queries(query_dir)
    n_queries = len(query_files)
    n_skipped = (
        count_completed_queries(run_config.output_dir, query_files)
        if args.resume
        else 0
    )
    workers, threads_per_worker = plan_workers(
        max(1, n_queries - n_skipped),
        run_config.threads,
        run_config.max_workers,
        run_config.threads_per_worker,
    )
    return PipelineContext(
        config=replace(run_config, query_path=query_dir),
        config_file=config_file,
        contigs=contigs,
        n_queries=n_queries,
        n_input_files=n_input_files,
        n_skipped=n_skipped,
        workers=workers,
        threads_per_worker=threads_per_worker,
    )


def thread_environment(
    environ: Mapping[str, str], per_marker_threads: int
) -> dict[str, str]:
    """Return ``environ`` with the thread variables of the numeric libraries set.

    Each gets ``per_marker_threads``, unless the inherited value is a smaller
    integer: a user who limits the threads further keeps that limit. Without the
    clamp, each inner library (pyhmmer, pyswrd, pyfamsa, NumPy/BLAS) would start up
    to ``threads_per_worker`` threads of its own and oversubscribe a large host.
    """
    env = dict(environ)
    for name in THREAD_ENV_VARS:
        inherited = env.get(name, "")
        if not (inherited.isdecimal() and 0 < int(inherited) < per_marker_threads):
            env[name] = str(per_marker_threads)
    return env


def build_runtime_env(repo_dir: Path, per_marker_threads: int) -> dict[str, str]:
    """Build the environment of the runner process."""
    env = thread_environment(os.environ, per_marker_threads)
    existing_pythonpath = env.get("PYTHONPATH", "")
    repo_str = str(repo_dir)
    if existing_pythonpath:
        if repo_str not in existing_pythonpath.split(os.pathsep):
            env["PYTHONPATH"] = os.pathsep.join([repo_str, existing_pythonpath])
    else:
        env["PYTHONPATH"] = repo_str
    env["PYTHONWARNINGS"] = "ignore"
    return env


@contextmanager
def _sigterm_as_keyboard_interrupt() -> Iterator[None]:
    """Make SIGTERM raise ``KeyboardInterrupt``, so it takes the path of Ctrl-C."""

    def raise_interrupt(_signum: int, _frame: FrameType | None) -> None:
        raise KeyboardInterrupt

    previous = signal.signal(signal.SIGTERM, raise_interrupt)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


@contextmanager
def _contig_workdir(contigs: ContigInput) -> Iterator[Path | None]:
    """Give a ``--contigs`` run its temporary folder, removed on every way out."""
    if not contigs.enabled:
        yield None
        return
    with tempfile.TemporaryDirectory(prefix="gvclass_contigs_") as name:
        yield Path(name)


def _process_group_is_alive(process: subprocess.Popen) -> bool:
    """Tell whether the process group that ``process`` leads has a member left."""
    process.poll()  # reaps an exited leader, which would otherwise count as a member
    try:
        os.killpg(process.pid, 0)
    except ProcessLookupError:
        return False
    return True


def _stop_process_group(process: subprocess.Popen) -> None:
    """End the whole process group of ``process``: SIGTERM, a grace period, SIGKILL.

    The grace period waits for every member of the group, not only the leader: a
    descendant that ignores SIGTERM outlives a leader that exits, and it is killed
    when the period ends.
    """
    with suppress(ProcessLookupError):
        os.killpg(process.pid, signal.SIGTERM)
    deadline = time.monotonic() + STOP_GRACE_S
    while _process_group_is_alive(process) and time.monotonic() < deadline:
        time.sleep(STOP_POLL_S)
    if _process_group_is_alive(process):
        with suppress(ProcessLookupError):  # the last member may exit just now
            os.killpg(process.pid, signal.SIGKILL)
    process.wait()


def _collect_output(
    process: subprocess.Popen, monitor: ResourceMonitor, tail: deque[str]
) -> None:
    """Read the runner output to its end; progress events go to the monitor."""
    if process.stdout is None:
        return
    for line in process.stdout:
        if not monitor.handle_progress_line(line):
            tail.append(line)


def run_runner(
    command: Sequence[str],
    *,
    env: Mapping[str, str],
    cwd: Path,
    output_dir: Path,
    queries_to_process: int,
    verbose: bool,
) -> RunnerOutcome:
    """Run the runner process and return how it ended.

    The runner has its own process group. Ctrl-C or SIGTERM ends that group, also
    when the runner itself has already exited, and returns an interrupted outcome.
    In the default mode the output goes through a pipe, and the outcome keeps its
    last non-progress lines; the full log is in ``gvclass.log``. Memory is sampled
    on a timer by a daemon thread.
    """
    child_env = dict(env) if verbose else {**env, "GVCLASS_PROGRESS_EVENTS": "1"}
    process = subprocess.Popen(
        list(command),
        env=child_env,
        cwd=cwd,
        stdout=None if verbose else subprocess.PIPE,
        stderr=None if verbose else subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        bufsize=1,
        start_new_session=True,
    )
    monitor = ResourceMonitor(output_dir, queries_to_process, process.pid)
    sampler = threading.Thread(
        target=monitor.monitor_memory, name=MONITOR_THREAD_NAME, daemon=True
    )
    sampler.start()
    tail: deque[str] = deque(maxlen=RUNNER_TAIL_LINES)
    interrupted = False
    try:
        with _sigterm_as_keyboard_interrupt():
            _collect_output(process, monitor, tail)
            process.wait()
    except KeyboardInterrupt:
        interrupted = True
    finally:
        if interrupted or process.poll() is None:
            _stop_process_group(process)
        monitor.stop()
        sampler.join(timeout=2.0)
    monitor.finish(process.returncode == 0)
    return RunnerOutcome(
        return_code=process.returncode,
        tail=tuple(tail),
        interrupted=interrupted,
        peak_memory_mb=monitor.peak_memory_mb,
    )


def validate_outputs(output_dir: Path, output: CliOutput, *, verbose: bool) -> bool:
    """Validate the pipeline outputs, print the issues, and say whether all is well."""
    from src.utils.output_validator import (  # noqa: PLC0415  # Pipeline code loads only for a run.
        validate_pipeline_outputs,
    )

    output.line("Validating pipeline outputs...", key="validate")
    validation = validate_pipeline_outputs(output_dir=output_dir, verbose=verbose)
    if validation["success"]:
        output.line("All outputs validated successfully", key="success")
        return True
    issues = validation["issues"]
    output.line(f"Validation found {len(issues)} issue(s):", key="warning")
    for issue in issues[:MAX_ISSUES_SHOWN]:
        print(f"  - {issue}")
    if len(issues) > MAX_ISSUES_SHOWN:
        print(f"  ... and {len(issues) - MAX_ISSUES_SHOWN} more")
    return False


def report_failure(outcome: RunnerOutcome, output: CliOutput) -> None:
    """Say how the runner failed and show the end of its output."""
    if outcome.interrupted:
        output.line("\nPipeline interrupted by user", key="warning")
        return
    output.line(f"\nPipeline failed with exit code: {outcome.return_code}", key="error")
    if outcome.tail:
        output.line(
            f"Last {len(outcome.tail)} lines of the runner output "
            "(the full log is gvclass.log in the output directory):"
        )
        print("".join(outcome.tail), end="")


def report_runtime(
    elapsed_s: float, peak_memory_mb: float, output_dir: Path, output: CliOutput
) -> None:
    """Print the run time, the peak memory and where the results are."""
    minutes, seconds = divmod(int(elapsed_s), 60)
    output.line(f"Runtime: {minutes}m {seconds}s", key="runtime")
    output.line(f"Peak memory usage: {peak_memory_mb:.0f} MB", key="memory")
    output.line(f"Results saved to: {output_dir}", key="summary")


def describe_contig_screen(config: RunConfig) -> str:
    """Say what the run does about the contig screen, for the startup lines."""
    # Imported here so that --help does not load numpy.
    from src.contig_screen.layout import bundle_present  # noqa: PLC0415

    if not config.contig_screen:
        return "disabled"
    if config.contigs_mode:
        return "not applicable in --contigs mode"
    try:
        installed = bundle_present(config.database_path)
    except OSError:
        return "add-on manifest not readable (screen_status will be failed)"
    if not installed:
        return "add-on not installed (install it with: gvclass screen --install)"
    return "enabled"


def display_startup_configuration(context: PipelineContext, output: CliOutput) -> None:
    """Print the banner and the settings of the run."""
    print_banner(output.plain_output, SOFTWARE_VERSION)
    print_configuration(context)
    output.line(f"Completeness mode: {context.config.completeness_mode}")
    output.line(f"Contig screen: {describe_contig_screen(context.config)}")


def execute_pipeline(
    args: argparse.Namespace,
    context: PipelineContext,
    repo_dir: Path,
    output: CliOutput,
) -> int:
    """Write the run configuration, run the runner on it and check the outputs.

    Return 0 only if the runner exited 0 and the output validation found no issue.
    """
    config = context.config
    config.output_dir.mkdir(parents=True, exist_ok=True)
    config_path = config.output_dir / RUN_CONFIG_FILENAME
    config_path.write_text(config.to_json(), encoding="utf-8")
    display_startup_configuration(context, output)

    command = [sys.executable, "-m", "src.bin.gvclass_runner", str(config_path)]
    if args.verbose:
        command.append("--verbose")
    # Half of the threads of a worker mirrors the per-marker thread cap that the
    # engine applies (see _calculate_marker_threading).
    env = build_runtime_env(repo_dir, max(1, context.threads_per_worker // 2))
    started = time.monotonic()
    outcome = run_runner(
        command,
        env=env,
        cwd=repo_dir,
        output_dir=config.output_dir,
        queries_to_process=context.n_remaining,
        verbose=args.verbose,
    )
    if outcome.interrupted or outcome.return_code != 0:
        report_failure(outcome, output)
        return 1

    output.line("\nPipeline completed successfully!", key="success")
    valid = validate_outputs(config.output_dir, output, verbose=args.verbose)
    report_runtime(
        time.monotonic() - started, outcome.peak_memory_mb, config.output_dir, output
    )
    return 0 if valid else 1


def run_cli(args: argparse.Namespace, repo_dir: Path, output: CliOutput) -> int:
    """Resolve the invocation, set up the database and run; return the exit status."""
    config_file = find_config_file(args.config, repo_dir)
    config = load_config(config_file)
    contigs = resolve_contigs_input(Path(args.query_dir).resolve(), args.contigs)
    run_config = resolve_run_config(args, config, repo_dir, contigs)
    apply_resource_cache(args, config, run_config.database_path)
    if contigs.enabled and args.min_length is not None:
        output.line(
            "Warning: --min-length is ignored in --contigs mode; "
            "use --contigs-min-length for split contigs.",
            key="warning",
            file=sys.stderr,
        )
    with _sigterm_as_keyboard_interrupt(), _contig_workdir(contigs) as workdir:
        if not check_and_setup_database(run_config.database_path, config, output):
            output.line("Database setup failed", key="error")
            return 1
        n_input_files = 0
        if workdir is not None:
            n_input_files = split_contig_inputs(
                contigs, output, run_config.contigs_min_length, workdir
            )
        context = plan_run(
            args,
            run_config,
            workdir or run_config.query_path,
            config_file=config_file,
            contigs=contigs,
            n_input_files=n_input_files,
        )
        return execute_pipeline(args, context, repo_dir, output)


def show_version(args: argparse.Namespace, repo_dir: Path) -> None:
    """Print the software version and the installed and latest database versions."""
    from src.utils.database_manager import (  # noqa: PLC0415  # Pipeline code loads only when asked.
        DatabaseManager,
    )

    config = load_config(find_config_file(args.config, repo_dir))
    db_path = resolve_database_path(args, config, repo_dir)

    print("GVClass Pipeline")
    print(f"  Software version: {SOFTWARE_VERSION}")
    if db_path.exists():
        print(f"  Database version: {DatabaseManager.get_database_version(db_path)}")
    else:
        print("  Database version: not installed")

    latest_source = DatabaseManager.get_latest_database_source(
        resolve_database_download_source(config)
    )
    if latest_source:
        print(f"  Latest Zenodo database: {latest_source['version']}")


def _run_screen(argv: list[str]) -> int:
    """Run ``gvclass screen``; the screens import torch only when they score."""
    from src.bin.gvclass_screen import main as screen_main  # noqa: PLC0415

    return screen_main(argv)


def _run_command(args: argparse.Namespace, repo_dir: Path, output: CliOutput) -> int:
    """Show the version, or run the pipeline."""
    if args.version:
        show_version(args, repo_dir)
        return 0
    return run_cli(args, repo_dir, output)


def main(argv: Sequence[str] | None = None, repo_dir: Path | None = None) -> int:
    """Run the GVClass command line and return the exit status.

    Container mode is the launcher's business: it decides before this module loads.
    """
    given = list(sys.argv[1:] if argv is None else argv)
    if given[:1] == ["screen"]:
        return _run_screen(given[1:])

    parser = build_parser()
    args = parser.parse_args(given)
    if not args.version and not args.query_dir:
        parser.print_help()
        return 0

    active_repo = repo_dir or Path(__file__).resolve().parents[2]
    output = CliOutput(plain_output=resolve_plain_output(args))
    try:
        return _run_command(args, active_repo, output)
    except ConfigError as exc:
        output.line(f"error: {exc}", key="error", file=sys.stderr)
        return 2
    except ValueError as exc:
        output.line(str(exc), key="error", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        output.line("\nPipeline interrupted by user", key="warning")
        return 1
    except Exception as exc:
        output.line(f"Pipeline failed: {exc}", key="error")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
