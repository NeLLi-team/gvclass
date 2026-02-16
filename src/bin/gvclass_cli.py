#!/usr/bin/env python
"""GVClass CLI runtime implementation."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from src.bin.cli_display import print_banner, print_configuration
from src.bin.progress_monitor import ResourceMonitor


SOFTWARE_VERSION = "v1.2.1"
PLAIN_OUTPUT_ENV = "GVCLASS_PLAIN_OUTPUT"
DATABASE_PATH_ENV = "GVCLASS_DB"
TWO_DECIMAL_SUMMARY_COLUMNS = {
    "avgdist",
    "order_dup",
    "order_completeness",
    "weighted_order_completeness",
    "order_weighted_completeness",
    "order_confidence_score",
    "gvog8_dup",
    "vp_df",
    "mirus_df",
    "cellular_dup",
    "GCperc",
    "CODINGperc",
}
NUMERIC_VALUE_PATTERN = re.compile(r"^[+-]?\d+(?:\.\d+)?$")


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


@dataclass
class ContigInput:
    """Resolved contig splitting input."""

    enabled: bool
    input_file: Optional[Path] = None
    input_dir: Optional[Path] = None


@dataclass
class WorkerPlan:
    """Parallel worker allocation."""

    workers: int
    threads_per_worker: int


@dataclass
class PipelineContext:
    """Resolved runtime context for one CLI invocation."""

    query_dir: Path
    active_query_dir: Path
    output_dir: Path
    database: Path
    threads: int
    tree_method: str
    mode_fast: bool
    sensitive_mode: bool
    contigs: ContigInput
    temp_contigs_dir: Optional[Path]
    n_queries: int
    n_input_files: int
    n_skipped: int
    workers: WorkerPlan


def count_query_files(directory: Path) -> int:
    if not directory.exists():
        return 0
    fna_files = list(directory.glob("*.fna"))
    faa_files = list(directory.glob("*.faa"))
    return len(fna_files) + len(faa_files)


def _add_core_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "query_dir", nargs="?", help="Query directory containing .fna or .faa files"
    )
    parser.add_argument(
        "-o", "--output-dir", help="Output directory (default: <query_dir>_results)"
    )
    parser.add_argument(
        "-c",
        "--config",
        default="config/gvclass_config.yaml",
        help="Configuration file (default: config/gvclass_config.yaml)",
    )
    parser.add_argument(
        "-d",
        "--database",
        help=f"Database path (overrides {DATABASE_PATH_ENV} and config)",
    )
    parser.add_argument(
        "-t", "--threads", type=int, help="Number of threads (overrides config)"
    )
    parser.add_argument("-j", "--max-workers", type=int, help="Maximum parallel workers")
    parser.add_argument("--threads-per-worker", type=int, help="Threads per worker")
    parser.add_argument(
        "--tree-method",
        choices=["fasttree", "iqtree"],
        help="Tree building method (overrides config)",
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
        help="Sensitive HMM mode: use E-value 1e-5 instead of GA cutoffs",
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
        "--resume",
        action="store_true",
        help="Resume from previous run, skipping completed queries",
    )


def _add_cluster_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cluster-type",
        default="local",
        choices=["local", "slurm", "pbs", "sge"],
        help="Type of cluster to use",
    )
    parser.add_argument("--cluster-queue", help="Queue/partition for cluster jobs")
    parser.add_argument("--cluster-project", help="Project/account for cluster billing")
    parser.add_argument(
        "--cluster-walltime", default="04:00:00", help="Walltime for cluster jobs"
    )


def _add_misc_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--version", action="store_true", help="Show version information")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run GVClass pipeline with clean output")
    _add_core_arguments(parser)
    _add_mode_arguments(parser)
    _add_cluster_arguments(parser)
    _add_misc_arguments(parser)
    return parser


def resolve_plain_output(args) -> bool:
    if args.plain_output:
        return True
    env_value = os.environ.get(PLAIN_OUTPUT_ENV, "").strip().lower()
    return env_value in {"1", "true", "yes", "on"}


def load_config(config_file: str, repo_dir: Path, output: CliOutput):
    config_path = Path(config_file)
    default_config = {
        "database": {
            "path": str(repo_dir / "resources"),
            "download_url": "https://zenodo.org/records/18662446/files/resources_v1_2_1.tar.gz?download=1",
            "expected_size": 1565,
        },
        "pipeline": {
            "tree_method": "fasttree",
            "mode_fast": True,
            "sensitive_mode": False,
            "threads": 16,
            "output_pattern": "{query_dir}_results",
        },
    }

    if not config_path.is_absolute():
        config_path = Path.cwd() / config_file
        if not config_path.exists():
            config_path = Path.cwd() / "gvclass_config.yaml"
        if not config_path.exists():
            config_path = repo_dir / config_file

    if not config_path.exists():
        return default_config

    try:
        import yaml
    except ModuleNotFoundError:
        output.line(
            "Warning: PyYAML is not available; using built-in defaults instead "
            f"of {config_path}.",
            file=sys.stderr,
        )
        return default_config

    with open(config_path, "r") as handle:
        return yaml.safe_load(handle)


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


def resolve_tree_method(args, config) -> str:
    if args.tree_method:
        return args.tree_method
    return config["pipeline"].get("tree_method", "fasttree")


def resolve_mode_fast(args, config) -> bool:
    if args.extended:
        return False
    if args.mode_fast:
        return True
    return config["pipeline"].get("mode_fast", True)


def resolve_sensitive_mode(args, config) -> bool:
    return args.sensitive or config["pipeline"].get("sensitive_mode", False)


def resolve_output_dir(args, query_dir: Path, config) -> Path:
    if args.output_dir:
        return Path(args.output_dir).resolve()

    output_pattern = config["pipeline"].get("output_pattern", "{query_dir}_results")
    query_base = query_dir.name
    return Path(output_pattern.format(query_dir=query_base)).resolve()


def resolve_contigs_input(query_dir_abs: Path, enabled: bool) -> ContigInput:
    if not enabled:
        return ContigInput(enabled=False)

    if query_dir_abs.is_file():
        if query_dir_abs.suffix.lower() not in [".fna", ".fa", ".fasta"]:
            raise ValueError(
                f"--contigs requires FNA file(s) (.fna, .fa, .fasta): {query_dir_abs}"
            )
        return ContigInput(enabled=True, input_file=query_dir_abs)

    if query_dir_abs.is_dir():
        return ContigInput(enabled=True, input_dir=query_dir_abs)

    raise ValueError(f"--contigs input not found: {query_dir_abs}")


def split_contig_inputs(contigs: ContigInput, output: CliOutput):
    if not contigs.enabled:
        return None, 0, 0

    from src.utils.contig_splitter import split_contigs, split_contigs_from_directory

    if contigs.input_file:
        prefix = contigs.input_file.stem
        output.line(f"Splitting contigs from: {contigs.input_file.name}", key="contigs")
        temp_dir, n_contigs = split_contigs(contigs.input_file, prefix=prefix)
        output.line(f"Found {n_contigs} contigs", key="success")
        return temp_dir, n_contigs, 1

    output.line(f"Splitting contigs from directory: {contigs.input_dir}", key="contigs")
    temp_dir, n_contigs, n_files = split_contigs_from_directory(contigs.input_dir)
    output.line(f"Found {n_contigs} contigs from {n_files} files", key="success")
    return temp_dir, n_contigs, n_files


def check_and_setup_database(db_path: Path, output: CliOutput) -> bool:
    required_checks = [
        db_path / "models" / "combined.hmm",
        db_path / "database",
        db_path / "gvclassFeb26_labels.tsv",
    ]
    if db_path.exists() and all(check.exists() for check in required_checks):
        output.line(f"Database found at: {db_path}", key="success")
        return True

    output.line(f"\nDatabase not found at {db_path}", key="database")
    output.line("Setting up database...", key="download")

    try:
        from src.utils.database_manager import DatabaseManager

        DatabaseManager.setup_database(str(db_path))
        output.line("Database setup complete!", key="success")
        return True
    except Exception as exc:
        output.line(f"\nFailed to setup database: {exc}", key="error")
        output.line("You can try manually with: pixi run setup-db")
        return False


def count_resume_skips(output_dir: Path) -> int:
    summary_files = set(
        path.name.replace(".summary.tab", "")
        for path in output_dir.glob("*.summary.tab")
    )
    tar_files = set(path.name.replace(".tar.gz", "") for path in output_dir.glob("*.tar.gz"))
    return len(summary_files & tar_files)


def calculate_worker_plan(args, n_queries: int, threads: int) -> WorkerPlan:
    if args.max_workers and args.threads_per_worker:
        return WorkerPlan(workers=args.max_workers, threads_per_worker=args.threads_per_worker)

    if args.max_workers:
        workers = args.max_workers
        return WorkerPlan(workers=workers, threads_per_worker=max(1, threads // workers))

    if args.threads_per_worker:
        workers = min(n_queries, max(1, threads // args.threads_per_worker))
        return WorkerPlan(workers=workers, threads_per_worker=args.threads_per_worker)

    if n_queries <= 4 and threads >= n_queries * 2:
        workers = n_queries
        return WorkerPlan(workers=workers, threads_per_worker=max(1, threads // workers))

    threads_per_worker = 4
    workers = min(n_queries, max(1, threads // threads_per_worker))
    return WorkerPlan(workers=workers, threads_per_worker=threads_per_worker)


def resolve_runtime_python(repo_dir: Path):
    pixi_python = repo_dir / ".pixi" / "envs" / "default" / "bin" / "python"
    python_exe = str(pixi_python) if pixi_python.exists() else sys.executable
    pixi_cmd = shutil.which("pixi")
    use_pixi_run = (not pixi_python.exists()) and bool(pixi_cmd)
    return python_exe, pixi_cmd, use_pixi_run


def check_missing_dependencies(repo_dir: Path, python_exe: str, pixi_cmd: Optional[str], use_pixi_run: bool):
    missing = []
    for dep in ["click"]:
        if use_pixi_run and pixi_cmd:
            result = subprocess.run(
                [pixi_cmd, "run", "python", "-c", f"import {dep}"],
                capture_output=True,
                text=True,
                cwd=repo_dir,
            )
        else:
            result = subprocess.run(
                [python_exe, "-c", f"import {dep}"], capture_output=True, text=True
            )
        if result.returncode != 0:
            missing.append(dep)
    return missing


def ensure_dependencies(repo_dir: Path, python_exe: str, pixi_cmd: Optional[str], use_pixi_run: bool, output: CliOutput):
    missing = check_missing_dependencies(repo_dir, python_exe, pixi_cmd, use_pixi_run)
    if missing and pixi_cmd:
        output.line("Installing dependencies with pixi...", key="database")
        subprocess.run([pixi_cmd, "install"], cwd=repo_dir)
        missing = check_missing_dependencies(repo_dir, python_exe, pixi_cmd, use_pixi_run)

    if not missing:
        return True

    output.line(
        "Missing Python dependencies in the runtime environment: " + ", ".join(missing),
        key="error",
    )
    output.line("Run `pixi install` or use `pixi run gvclass ...`.")
    return False


def _base_prefect_command(
    python_exe: str, pixi_cmd: Optional[str], use_pixi_run: bool
) -> list[str]:
    if use_pixi_run and pixi_cmd:
        return [pixi_cmd, "run", "python", "-m", "src.bin.gvclass_prefect"]
    return [python_exe, "-m", "src.bin.gvclass_prefect"]


def _append_optional_pipeline_flags(
    cmd: list[str], args, mode_fast: bool, sensitive_mode: bool
) -> None:
    cmd.append("--mode-fast" if mode_fast else "--extended")
    if sensitive_mode:
        cmd.append("--sensitive")
    if args.cluster_queue:
        cmd.extend(["--cluster-queue", args.cluster_queue])
    if args.cluster_project:
        cmd.extend(["--cluster-project", args.cluster_project])
    if args.cluster_walltime:
        cmd.extend(["--cluster-walltime", args.cluster_walltime])
    if args.verbose:
        cmd.append("--verbose")
    if args.resume:
        cmd.append("--resume")


def build_pipeline_command(
    args,
    context: PipelineContext,
    python_exe: str,
    pixi_cmd: Optional[str],
    use_pixi_run: bool,
) -> list[str]:
    cmd = _base_prefect_command(python_exe, pixi_cmd, use_pixi_run)
    cmd.extend(
        [
            str(context.active_query_dir),
            "-o",
            str(context.output_dir),
            "-d",
            str(context.database),
            "-t",
            str(context.threads),
            "-j",
            str(context.workers.workers),
            "--threads-per-worker",
            str(context.workers.threads_per_worker),
            "--tree-method",
            context.tree_method,
            "--cluster-type",
            args.cluster_type,
        ]
    )
    _append_optional_pipeline_flags(cmd, args, context.mode_fast, context.sensitive_mode)
    return cmd


def build_runtime_env(repo_dir: Path):
    env = os.environ.copy()
    existing_pythonpath = env.get("PYTHONPATH", "")
    repo_str = str(repo_dir)
    if existing_pythonpath:
        if repo_str not in existing_pythonpath.split(os.pathsep):
            env["PYTHONPATH"] = os.pathsep.join([repo_str, existing_pythonpath])
    else:
        env["PYTHONPATH"] = repo_str
    env["PREFECT_API_URL"] = ""
    env["PREFECT_LOGGING_LEVEL"] = "ERROR"
    env["PYTHONWARNINGS"] = "ignore"
    env["PREFECT_SERVER_STARTUP_TIMEOUT"] = "300"
    env["PREFECT_API_DATABASE_TIMEOUT"] = "60"
    return env


def combine_summary_files(output_dir: Path, output: CliOutput):
    summary_files = list(output_dir.glob("*.summary.tab"))
    if not summary_files:
        output.line("No individual summary files found", key="warning")
        return

    combined_rows = []
    headers = None
    for summary_file in sorted(summary_files):
        with open(summary_file, "r", newline="") as handle:
            rows = list(csv.reader(handle, delimiter="\t"))
        if not rows:
            continue

        if headers is None and rows[0]:
            headers = rows[0]
        elif headers is not None and rows[0] != headers:
            output.line(
                f"Header mismatch in {summary_file.name}; using first header",
                key="warning",
            )

        for row in rows[1:]:
            if row:
                combined_rows.append(row)

    if not (combined_rows and headers):
        output.line("No summary data found to combine", key="warning")
        return

    combined_tsv = output_dir / "gvclass_summary.tsv"
    combined_csv = output_dir / "gvclass_summary.csv"
    with open(combined_tsv, "w") as tsv_handle, open(combined_csv, "w", newline="") as csv_handle:
        writer = csv.writer(csv_handle)
        tsv_handle.write("\t".join(headers) + "\n")
        writer.writerow(headers)
        for row in combined_rows:
            formatted_row = [
                _format_summary_cell(header, value)
                for header, value in zip(headers, row)
            ]
            tsv_handle.write("\t".join(formatted_row) + "\n")
            writer.writerow(formatted_row)

    output.line(f"Combined summary written to: {combined_tsv}", key="success")
    output.line(f"CSV summary written to: {combined_csv}", key="success")


def _format_summary_cell(header: str, value: str) -> str:
    if header not in TWO_DECIMAL_SUMMARY_COLUMNS:
        return value
    if not value or not NUMERIC_VALUE_PATTERN.match(value):
        return value
    return f"{float(value):.2f}"


def cleanup_temp_dir(temp_dir: Optional[Path]):
    if temp_dir and Path(temp_dir).exists():
        shutil.rmtree(temp_dir)


def run_pipeline(
    cmd,
    env,
    repo_dir: Path,
    output_dir: Path,
    queries_to_process: int,
    verbose: bool,
):
    if verbose:
        process = subprocess.Popen(cmd, env=env, cwd=repo_dir)
    else:
        process = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=repo_dir,
        )

    existing_completed = len(list(output_dir.glob("*.summary.tab")))
    monitor = ResourceMonitor(output_dir, queries_to_process, process.pid)
    monitor.current_query = existing_completed
    monitor_thread = threading.Thread(target=monitor.monitor)
    monitor_thread.start()

    stdout_data, stderr_data = process.communicate()
    return_code = process.returncode
    monitor.stop()
    monitor_thread.join()
    print("\r" + " " * 80 + "\r", end="")
    return return_code, stdout_data, stderr_data, monitor


def show_version(config_path: str, database_override: Optional[str], repo_dir: Path):
    from src.utils.database_manager import DatabaseManager

    output = CliOutput(plain_output=True)
    config = load_config(config_path, repo_dir, output)
    db_path = Path(resolve_database_setting(database_override, config)).expanduser()
    if not db_path.is_absolute():
        db_path = repo_dir / db_path
    db_path = db_path.resolve()

    print("GVClass Pipeline")
    print(f"  Software version: {SOFTWARE_VERSION}")
    if db_path.exists():
        print(f"  Database version: {DatabaseManager.get_database_version(db_path)}")
    else:
        print("  Database version: not installed")


def resolve_pipeline_context(args, repo_dir: Path, output: CliOutput) -> PipelineContext:
    config = load_config(args.config, repo_dir, output)
    query_dir = Path(args.query_dir).resolve()

    contigs = resolve_contigs_input(query_dir, args.contigs)
    database = resolve_database_path(args, config, repo_dir)
    threads = args.threads if args.threads else config["pipeline"].get("threads", 16)
    tree_method = resolve_tree_method(args, config)
    mode_fast = resolve_mode_fast(args, config)
    sensitive_mode = resolve_sensitive_mode(args, config)
    output_dir = resolve_output_dir(args, query_dir, config)

    if not check_and_setup_database(database, output):
        raise RuntimeError("Database setup failed")

    temp_contigs_dir, _, n_input_files = split_contig_inputs(contigs, output)
    active_query_dir = temp_contigs_dir if temp_contigs_dir else query_dir

    n_queries = count_query_files(active_query_dir)
    if n_queries == 0:
        raise ValueError(f"No .fna or .faa files found in {active_query_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    n_skipped = count_resume_skips(output_dir) if args.resume else 0
    workers = calculate_worker_plan(args, n_queries, threads)

    return PipelineContext(
        query_dir=query_dir,
        active_query_dir=active_query_dir,
        output_dir=output_dir,
        database=database,
        threads=threads,
        tree_method=tree_method,
        mode_fast=mode_fast,
        sensitive_mode=sensitive_mode,
        contigs=contigs,
        temp_contigs_dir=temp_contigs_dir,
        n_queries=n_queries,
        n_input_files=n_input_files,
        n_skipped=n_skipped,
        workers=workers,
    )


def _emit_subprocess_output(data: Optional[bytes]) -> None:
    if data:
        print(data.decode(errors="replace").strip())


def handle_failed_pipeline(
    return_code: int,
    stdout_data: Optional[bytes],
    stderr_data: Optional[bytes],
    verbose: bool,
    output: CliOutput,
) -> int:
    output.line(f"\nPipeline failed with exit code: {return_code}", key="error")
    if not verbose:
        _emit_subprocess_output(stderr_data)
        _emit_subprocess_output(stdout_data)
    return 1


def validate_and_summarize_outputs(
    output_dir: Path, verbose: bool, output: CliOutput
) -> None:
    from src.utils.output_validator import validate_pipeline_outputs

    output.line("Validating pipeline outputs...", key="validate")
    validation = validate_pipeline_outputs(output_dir=output_dir, verbose=verbose)
    if validation["success"]:
        output.line("All outputs validated successfully", key="success")
    else:
        output.line(f"Validation found {len(validation['issues'])} issue(s):", key="warning")
        for issue in validation["issues"][:5]:
            print(f"  - {issue}")

    output.line("Generating combined summary...", key="summary")
    combine_summary_files(output_dir, output)


def report_runtime(
    start_time: datetime, monitor: ResourceMonitor, output_dir: Path, output: CliOutput
) -> None:
    elapsed = datetime.now() - start_time
    minutes = int(elapsed.total_seconds() / 60)
    seconds = int(elapsed.total_seconds() % 60)
    output.line(f"Runtime: {minutes}m {seconds}s", key="runtime")
    output.line(f"Peak memory usage: {monitor.peak_memory_mb:.0f} MB", key="memory")
    output.line(f"Results saved to: {output_dir}", key="summary")


def display_startup_configuration(args, context: PipelineContext, output: CliOutput) -> None:
    print_banner(output.plain_output, SOFTWARE_VERSION)
    print_configuration(
        args=args,
        query_dir_abs=context.active_query_dir,
        output_dir=context.output_dir,
        database=context.database,
        n_queries=context.n_queries,
        n_input_files=context.n_input_files,
        n_skipped=context.n_skipped,
        mode_fast=context.mode_fast,
        sensitive_mode=context.sensitive_mode,
        threads=context.threads,
        workers=context.workers,
        contigs=context.contigs,
    )


def build_runtime_command(args, context: PipelineContext, repo_dir: Path, output: CliOutput):
    python_exe, pixi_cmd, use_pixi_run = resolve_runtime_python(repo_dir)
    deps_ok = ensure_dependencies(repo_dir, python_exe, pixi_cmd, use_pixi_run, output)
    if not deps_ok:
        return None, None
    cmd = build_pipeline_command(
        args=args,
        context=context,
        python_exe=python_exe,
        pixi_cmd=pixi_cmd,
        use_pixi_run=use_pixi_run,
    )
    return cmd, build_runtime_env(repo_dir)


def execute_pipeline_command(
    args, context: PipelineContext, repo_dir: Path, cmd, env, output: CliOutput
) -> int:
    start_time = datetime.now()
    queries_to_process = (
        context.n_queries - context.n_skipped if args.resume else context.n_queries
    )

    try:
        return_code, stdout_data, stderr_data, monitor = run_pipeline(
            cmd=cmd,
            env=env,
            repo_dir=repo_dir,
            output_dir=context.output_dir,
            queries_to_process=queries_to_process,
            verbose=args.verbose,
        )
    except KeyboardInterrupt:
        output.line("\nPipeline interrupted by user", key="warning")
        return 1

    if return_code != 0:
        return handle_failed_pipeline(
            return_code=return_code,
            stdout_data=stdout_data,
            stderr_data=stderr_data,
            verbose=args.verbose,
            output=output,
        )

    output.line("\nPipeline completed successfully!", key="success")
    validate_and_summarize_outputs(context.output_dir, args.verbose, output)
    report_runtime(start_time, monitor, context.output_dir, output)
    return 0


def resolve_cli_context(args, repo_dir: Path, output: CliOutput) -> Optional[PipelineContext]:
    try:
        return resolve_pipeline_context(args, repo_dir, output)
    except ValueError as exc:
        output.line(f"ERROR: {exc}", key="error", file=sys.stderr)
        return None
    except RuntimeError as exc:
        output.line(str(exc), key="error")
        return None


def run_cli(args, repo_dir: Path, output: CliOutput) -> int:
    context = resolve_cli_context(args, repo_dir, output)
    if context is None:
        return 1

    display_startup_configuration(args, context, output)
    cmd, env = build_runtime_command(args, context, repo_dir, output)
    if not cmd:
        cleanup_temp_dir(context.temp_contigs_dir)
        return 1

    try:
        return execute_pipeline_command(args, context, repo_dir, cmd, env, output)
    finally:
        cleanup_temp_dir(context.temp_contigs_dir)

    # Unreachable; retained to satisfy static analysis in some tooling.
    return 0


def main(repo_dir: Optional[Path] = None) -> int:
    parser = build_parser()
    args = parser.parse_args()

    active_repo = repo_dir if repo_dir else Path(__file__).resolve().parents[2]
    if args.version:
        show_version(args.config, args.database, active_repo)
        return 0

    if not args.query_dir:
        parser.print_help()
        return 0

    plain_output = resolve_plain_output(args)
    output = CliOutput(plain_output=plain_output)

    try:
        return run_cli(args, active_repo, output)
    except ValueError as exc:
        output.line(str(exc), key="error", file=sys.stderr)
        return 1
    except Exception as exc:
        output.line(f"Pipeline failed: {exc}", key="error")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
