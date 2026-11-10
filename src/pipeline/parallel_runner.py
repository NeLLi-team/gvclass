"""GVClass pipeline execution engine."""

from __future__ import annotations

import contextlib
import logging
import traceback
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from typing import Any

from src.config.marker_sets import missing_panel_models
from src.core.alignment import configure_iqtree
from src.core.hmm_search import COMBINED_HMM
from src.pipeline.contig_screen_step import ScreenOptions, add_screen_columns
from src.pipeline.progress_events import emit_progress_event
from src.pipeline.query_processing_engine import WORK_DIRNAME, process_single_query
from src.pipeline.run_config import RunConfig, plan_workers
from src.pipeline.run_status import (
    COMPLETE,
    FAILED,
    INVALID_INPUT,
    QueryRecord,
    RunRecorder,
    resumable_queries,
)
from src.pipeline.summary_writer import (
    archive_final_summary_extended_files,
    write_final_summary_extended_files,
    write_final_summary_files,
)
from src.utils.common import sha256_file
from src.utils.database_manager import DatabaseManager
from src.utils.error_handling import FileError, ProcessingError, ValidationError
from src.utils.input_validation import (
    RESERVED_QUERY_STEMS,
    InputValidator,
    discover_query_files,
    is_safe_query_stem,
)

logger = logging.getLogger("gvclass_runner")


@dataclass(frozen=True, slots=True)
class RunInputs:
    """The checked inputs of a run.

    ``rejected_inputs`` pairs each invalid file with the reason; the other files
    are ``query_files``. ``database_path`` is the resolved, complete database.
    """

    query_files: tuple[Path, ...]
    rejected_inputs: tuple[tuple[Path, str], ...]
    database_path: Path


def validate_and_setup_task(config: RunConfig) -> RunInputs:
    """Validate inputs and setup directories.

    Each query file is validated on its own. Valid files become ``query_files``;
    each invalid file becomes a ``(path, error text)`` pair in ``rejected_inputs``
    so one bad file does not stop the batch. Colliding or reserved query names
    and a directory without any valid file are errors for the whole run, and so is a
    query name that is empty, starts with a dot or holds a path separator.
    """
    logger.info(
        "Validating inputs - Query: %s, Output: %s",
        config.query_path,
        config.output_dir,
    )
    query_path = InputValidator.validate_directory(config.query_path, must_exist=True)
    discovered_files = discover_query_files(query_path)
    if not discovered_files:
        raise ValidationError(
            f"No valid sequence files found in query directory: {config.query_path}",
            field="query_directory",
            value=str(config.query_path),
        )
    _reject_unsafe_stems(discovered_files)
    _reject_colliding_stems(discovered_files)
    _reject_reserved_stems(discovered_files)

    query_files: list[Path] = []
    rejected_inputs: list[tuple[Path, str]] = []
    for query_file in discovered_files:
        error = _validation_error(
            query_file,
            allow_short=config.allow_short,
            min_sequence_length=config.min_sequence_length,
        )
        if error is None:
            query_files.append(query_file)
            continue
        logger.error("Rejected input %s: %s", query_file.name, error)
        rejected_inputs.append((query_file, error))
    if not query_files:
        details = "; ".join(f"{path.stem}: {error}" for path, error in rejected_inputs)
        raise ValidationError(
            f"No valid query file in {config.query_path}. {details}",
            field="query_directory",
            value=str(config.query_path),
        )

    config.output_dir.mkdir(parents=True, exist_ok=True)
    database_path = DatabaseManager.setup_database(config.database_path)
    logger.info("Using database at: %s", database_path)
    logger.info(
        "Found %d query files, %d rejected", len(query_files), len(rejected_inputs)
    )
    return RunInputs(tuple(query_files), tuple(rejected_inputs), database_path)


def _validation_error(
    query_file: Path, *, allow_short: bool, min_sequence_length: int | None
) -> str | None:
    """Return why ``query_file`` is invalid, or ``None`` when it is valid."""
    try:
        InputValidator.validate_sequence_file(
            query_file,
            allow_short=allow_short,
            min_sequence_length=min_sequence_length,
        )
    except (ValidationError, FileError) as exc:
        return exc.message
    return None


def _reject_unsafe_stems(query_files: list[Path]) -> None:
    unsafe = sorted(
        query_file.name
        for query_file in query_files
        if not is_safe_query_stem(query_file.stem)
    )
    if unsafe:
        raise ValidationError(
            "These query filenames cannot name a work directory: a name before the "
            "extension must not be empty, start with a dot or hold a slash. "
            f"Rename them: {', '.join(unsafe)}"
        )


def _reject_colliding_stems(query_files: list[Path]) -> None:
    files_by_stem: dict[str, list[Path]] = {}
    for query_file in query_files:
        files_by_stem.setdefault(query_file.stem, []).append(query_file)
    collisions = {
        stem: files for stem, files in files_by_stem.items() if len(files) > 1
    }
    if collisions:
        details = "; ".join(
            f"{stem}: {', '.join(path.name for path in files)}"
            for stem, files in sorted(collisions.items())
        )
        raise ValidationError(
            "Query filenames must have unique stems because stems identify output "
            f"files. Rename these colliding inputs: {details}"
        )


def _reject_reserved_stems(query_files: list[Path]) -> None:
    reserved = sorted(
        query_file.name
        for query_file in query_files
        if query_file.stem in RESERVED_QUERY_STEMS
    )
    if reserved:
        raise ValidationError(
            "These query filenames use a reserved name that GVClass needs for its "
            f"own outputs. Rename them: {', '.join(reserved)}"
        )


def check_marker_panels(database_path: Path) -> None:
    """Raise if a marker panel names a model that the database lacks.

    A panel that no model matches would count zero markers for every query without
    any message, so the run stops at the start instead.
    """
    models: set[str] = set()
    with (database_path / COMBINED_HMM).open(encoding="utf-8") as handle:
        models.update(
            line.split(maxsplit=1)[1].strip()
            for line in handle
            if line.startswith("NAME ")
        )
    missing = missing_panel_models(models)
    if missing:
        details = "; ".join(
            f"{panel}: {', '.join(names)}" for panel, names in sorted(missing.items())
        )
        raise ProcessingError(
            f"The database {database_path} lacks models that the marker panels "
            f"name ({details}). Use a database of the same release as GVClass.",
            step="marker_panels",
        )


def gvclass_flow(config: RunConfig) -> Path:
    """Run the GVClass pipeline for every query of ``config`` and return the summary."""
    configure_iqtree(config.iqtree_mode)
    logger.info("Validating inputs and setting up directories")
    inputs = validate_and_setup_task(config)
    config = replace(config, database_path=inputs.database_path)
    check_marker_panels(config.database_path)

    database_version = DatabaseManager.get_database_version(config.database_path)
    skipped = (
        resumable_queries(config, inputs.query_files, database_version)
        if config.resume
        else {}
    )
    to_run = [path for path in inputs.query_files if path.stem not in skipped]
    logger.info("Found %d queries to process", len(to_run))
    if config.resume:
        logger.info(
            "Resume mode: skipped %d completed queries, %d remaining",
            len(skipped),
            len(to_run),
        )

    recorder = RunRecorder(
        config,
        database_version,
        [path.stem for path in inputs.query_files]
        + [path.stem for path, _ in inputs.rejected_inputs],
        skipped,
    )
    recorder.start()
    try:
        for path, error in inputs.rejected_inputs:
            recorder.record_failed(path.stem, error, INVALID_INPUT)
        results = _process_new_queries(to_run, config, skipped, recorder)
        all_results = [
            *_skipped_results(skipped),
            *results,
            *(
                {"query": path.stem, "status": INVALID_INPUT, "error": error}
                for path, error in inputs.rejected_inputs
            ),
        ]
        return _summarize_run(config, inputs, to_run, results, all_results, recorder)
    except BaseException:
        recorder.finish(FAILED)
        raise


def _skipped_results(skipped: dict[str, QueryRecord]) -> list[dict[str, Any]]:
    """Return the results of the queries that resume reuses, from their records."""
    return [
        {"query": name, "status": COMPLETE, "summary_data": dict(record.summary or {})}
        for name, record in skipped.items()
    ]


def _process_new_queries(
    to_run: list[Path],
    config: RunConfig,
    skipped: dict[str, QueryRecord],
    recorder: RunRecorder,
) -> list[dict[str, Any]]:
    """Process the queries that resume did not skip, in a pool of threads."""
    if not to_run:
        logger.info("No queries remaining after the resume filter; nothing to process.")
        return []
    n_workers, threads_per_worker = plan_workers(
        len(to_run), config.threads, config.max_workers, config.threads_per_worker
    )
    logger.info(
        "Parallelization strategy: %d workers x %d threads",
        n_workers,
        threads_per_worker,
    )
    if config.species_tree:
        _init_species_tree(config, to_run, skipped)
    logger.info("Starting parallel query processing")
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(
            executor.map(
                partial(
                    _process_one,
                    config=config,
                    threads=threads_per_worker,
                    recorder=recorder,
                ),
                to_run,
            )
        )
    # A query that failed leaves its folder, so the folder stays then.
    with contextlib.suppress(OSError):
        (config.output_dir / WORK_DIRNAME).rmdir()
    completed = sum(1 for result in results if result["status"] == COMPLETE)
    logger.info(
        "All queries complete. Successfully processed %d/%d queries",
        completed,
        len(results),
    )
    return results


def _init_species_tree(
    config: RunConfig, to_run: list[Path], skipped: dict[str, QueryRecord]
) -> None:
    """Clear earlier species-tree scratch and outputs and record this run's queries.

    The per-query hook and the opt-in combined step then read exactly this run's
    sidecars. Setup failures must fail the requested run.
    """
    if config.resume and skipped:
        logger.warning(
            "%d resume-skipped queries will have no species-tree placement this run "
            "(per-query trees + species_tree_* columns are produced only for queries "
            "processed this run); re-run those without --resume to place them.",
            len(skipped),
        )
    from src.core.species_tree.handoff import (  # noqa: PLC0415  # Load optional species-tree code only when requested.
        init_handoff,
    )

    init_handoff(
        config.output_dir, [path.stem for path in to_run], resume=config.resume
    )


def _process_one(
    query_file: Path, config: RunConfig, threads: int, recorder: RunRecorder
) -> dict[str, Any]:
    """Process one query, record how it finished, and return its result."""
    name = query_file.stem
    recorder.record_query_started(name)
    try:
        input_sha256 = sha256_file(query_file)
        result = process_single_query(query_file, config, threads)
        _require_unchanged_input(query_file, input_sha256)
        recorder.record_complete(query_file, result, input_sha256)
    except Exception as exc:  # noqa: BLE001  # one query's failure is recorded; the others go on
        logger.error("Query %s failed: %s", name, exc)
        logger.error("Traceback: %s", traceback.format_exc())
        recorder.record_failed(name, str(exc))
        emit_progress_event(
            {
                "event": "query_failed",
                "query": name,
                "progress": 100,
                "stage": "failed",
            }
        )
        return {"query": name, "status": FAILED, "error": str(exc)}
    logger.info("Completed: %s - %s", name, result["status"])
    return result


def _require_unchanged_input(query_file: Path, input_sha256: str) -> None:
    """Raise if ``query_file`` no longer has the SHA-256 it had when the query began."""
    if sha256_file(query_file) != input_sha256:
        raise ProcessingError(
            f"The input file {query_file.name} changed while GVClass processed it, "
            "so the results do not match the file. Run the query again.",
            step="input_fingerprint",
            input_file=str(query_file),
        )


def _summarize_run(
    config: RunConfig,
    inputs: RunInputs,
    to_run: list[Path],
    results: list[dict[str, Any]],
    all_results: list[dict[str, Any]],
    recorder: RunRecorder,
) -> Path:
    """Screen the bins, write the summary tables and close the run."""
    screen_options = ScreenOptions(
        enabled=config.contig_screen,
        contigs_mode=config.contigs_mode,
        resume=config.resume,
        plan_workers=partial(
            plan_workers,
            total_threads=config.threads,
            max_workers=config.max_workers,
            threads_per_worker=config.threads_per_worker,
        ),
    )
    logger.info("Screening bins for contaminating contigs")
    add_screen_columns(
        all_results,
        {path.stem: path for path in inputs.query_files},
        config.output_dir,
        config.database_path,
        screen_options,
        logger,
    )
    summary_file = create_final_summary_task(all_results, config.output_dir)

    # Opt-in combined tree from this run's sidecars (--species-tree-combined). The
    # per-query hook already owns the summary columns, so this writes only the
    # additional <out>/species_tree/ combined artifacts. Failures fail the run.
    # Resume-skipped queries have no sidecar this run and stay out: a documented
    # limitation.
    if config.species_tree_combined and to_run:
        from src.core.species_tree.orchestration import (  # noqa: PLC0415  # Load optional species-tree code only when requested.
            run_combined_species_tree,
        )

        run_combined_species_tree(
            config.output_dir,
            config.database_path,
            results,
            config.threads,
            config.species_tree_trim,
            config.tree_method,
        )
    elif config.species_tree_combined:
        logger.info(
            "All queries resume-skipped; the combined species tree is built only "
            "from this run's queries, so it is not regenerated. Re-run without "
            "--resume to rebuild out/species_tree/combined.*"
        )

    failed_results = [result for result in all_results if result["status"] != COMPLETE]
    recorder.finish(FAILED if failed_results else COMPLETE)
    logger.info("Pipeline completed! Results in: %s", config.output_dir)
    logger.info("Summary written to: %s", summary_file)
    if failed_results:
        raise RuntimeError(
            f"GVClass failed for {len(failed_results)} of {len(all_results)} "
            f"queries; see gvclass_failed_queries.tsv and run_status.json in "
            f"{config.output_dir}."
        )
    return summary_file


def create_final_summary_task(results: list[dict[str, Any]], output_dir: Path) -> Path:
    """Create the final summary file combining all results."""
    logger.info("Creating final summary for %d queries", len(results))
    summary_tsv = write_final_summary_files(results, output_dir)
    extended_tsv = write_final_summary_extended_files(results, output_dir)
    extended_archive = archive_final_summary_extended_files(output_dir)
    logger.info("Summary written to: %s", summary_tsv)
    logger.info("CSV summary written to: %s", output_dir / "gvclass_summary.csv")
    logger.info(
        "Extended diagnostics archived in: %s", extended_archive or extended_tsv
    )
    return summary_tsv
