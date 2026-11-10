"""Output validation for the GVClass pipeline.

A run is valid when ``results.jsonl`` holds a ``complete`` last record for every
query of the run (the queries in ``run_status.json``) and each record's archive
exists with the recorded size and holds a summary and final summary with a data row
for its query. A missing or empty member is an issue, not a warning. A record of a
query that the run does not have, such as an input removed since an earlier run, is
not read.
"""

import csv
import io
import logging
import tarfile
from pathlib import Path
from typing import Any

from src.pipeline.run_status import (
    COMPLETE,
    RESULTS_FILENAME,
    RUN_STATUS_FILENAME,
    QueryRecord,
    archive_is_intact,
    read_results,
    read_run_query_names,
)

logger = logging.getLogger(__name__)

REQUIRED_COLUMNS = ("query", "taxonomy_majority")


def validate_pipeline_outputs(
    output_dir: Path, *, verbose: bool = False
) -> dict[str, Any]:
    """Validate that pipeline outputs are complete and properly formatted.

    Args:
        output_dir: Path to the output directory containing pipeline results
        verbose: Whether to log detailed validation information

    Returns:
        Dictionary with validation results:
        - 'success': Boolean indicating if validation passed
        - 'issues': List of issue descriptions (empty if success=True)
        - 'warnings': List of non-critical warnings
    """
    issues: list[str] = []
    warnings: list[str] = []

    if verbose:
        logger.info("Validating outputs in: %s", output_dir)

    if not output_dir.exists():
        issues.append(f"Output directory does not exist: {output_dir}")
        return {"success": False, "issues": issues, "warnings": warnings}

    query_names = read_run_query_names(output_dir)
    if query_names is None:
        issues.append(f"{RUN_STATUS_FILENAME} is missing or unreadable")
        query_names = []
    records = read_results(output_dir)
    for name in query_names:
        record = records.get(name)
        if record is None:
            issues.append(f"Query {name} has no record in {RESULTS_FILENAME}")
            continue
        _validate_record(output_dir, name, record, issues)

    combined_summary = output_dir / "gvclass_summary.tsv"
    if combined_summary.exists():
        _validate_combined_summary(combined_summary, issues, warnings)

    if verbose:
        _log_outcome(len(query_names), issues, warnings)
    return {"success": not issues, "issues": issues, "warnings": warnings}


def _log_outcome(n_queries: int, issues: list[str], warnings: list[str]) -> None:
    logger.info("Checked %d query record(s)", n_queries)
    if issues:
        logger.warning("Validation found %s issue(s)", len(issues))
    else:
        logger.info("Validation successful")
    if warnings:
        logger.warning("%s warning(s) found", len(warnings))


def _validate_record(
    output_dir: Path, name: str, record: QueryRecord, issues: list[str]
) -> None:
    """Add an issue for each way the last record of a query is not a good result."""
    if record.status != COMPLETE:
        issues.append(f"Query {name} is {record.status}: {record.error or 'no error'}")
        return
    if not archive_is_intact(output_dir, record):
        issues.append(f"Archive of query {name} is missing or has changed in size")
        return
    _validate_archive_summaries(output_dir / f"{name}.tar.gz", name, issues)


def _validate_archive_summaries(
    archive_file: Path, query_name: str, issues: list[str]
) -> None:
    """Add an issue for each summary member that is missing or lacks a query row."""
    members = (
        f"{query_name}/{query_name}.summary.tab",
        f"{query_name}/{query_name}.final_summary.tsv",
    )
    try:
        with tarfile.open(archive_file, "r:gz") as archive:
            for member in members:
                _validate_member(archive, member, query_name, issues)
    except (OSError, tarfile.TarError) as exc:
        issues.append(f"Invalid query archive {archive_file.name}: {exc}")


def _validate_member(
    archive: tarfile.TarFile, member: str, query_name: str, issues: list[str]
) -> None:
    try:
        extracted = archive.extractfile(member)
    except KeyError:
        issues.append(f"Archive member {member} is missing")
        return
    if extracted is None:
        issues.append(f"Archive member {member} is not a file")
        return
    text = extracted.read().decode("utf-8", errors="replace")
    table = list(csv.reader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE))
    header = table[0] if table else []
    rows = [row for row in table[1:] if any(cell.strip() for cell in row)]
    if not rows:
        issues.append(f"Archive member {member} is empty or has no data row")
        return
    missing = [column for column in REQUIRED_COLUMNS if column not in header]
    if missing:
        issues.append(
            f"Archive member {member} lacks required columns: {', '.join(missing)}"
        )
        return
    query_column = header.index("query")
    if not any(
        len(row) > query_column and row[query_column] == query_name for row in rows
    ):
        issues.append(f"Archive member {member} has no row for query {query_name}")


def _validate_combined_summary(
    combined_summary: Path, issues: list[str], warnings: list[str]
) -> None:
    """Add issues for an empty combined summary or one without the key columns."""
    lines = combined_summary.read_text(encoding="utf-8").splitlines()
    if not lines:
        issues.append("Combined summary file is empty")
        return
    missing = [
        column for column in REQUIRED_COLUMNS if column not in lines[0].split("\t")
    ]
    if missing:
        issues.append(
            f"Combined summary missing required columns: {', '.join(missing)}"
        )
        return
    if len(lines) == 1:
        warnings.append("Combined summary has header but no data rows")
