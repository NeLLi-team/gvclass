"""
Output validation utilities for GVClass pipeline.

This module provides functions to validate the outputs of the GVClass pipeline,
ensuring all expected files are present and properly formatted.
"""

from pathlib import Path
from typing import Dict, List, Any
import logging

logger = logging.getLogger(__name__)


def validate_pipeline_outputs(
    output_dir: Path, verbose: bool = False
) -> Dict[str, Any]:
    """
    Validate that pipeline outputs are complete and properly formatted.

    Args:
        output_dir: Path to the output directory containing pipeline results
        verbose: Whether to log detailed validation information

    Returns:
        Dictionary with validation results:
        - 'success': Boolean indicating if validation passed
        - 'issues': List of issue descriptions (empty if success=True)
        - 'warnings': List of non-critical warnings
    """
    issues = []
    warnings = []

    if verbose:
        logger.info(f"Validating outputs in: {output_dir}")

    # Check if output directory exists
    if not output_dir.exists():
        issues.append(f"Output directory does not exist: {output_dir}")
        return {"success": False, "issues": issues, "warnings": warnings}

    # Check for summary files (*.summary.tab)
    summary_files = list(output_dir.glob("*.summary.tab"))
    if not summary_files:
        issues.append("No summary.tab files found in output directory")
    else:
        if verbose:
            logger.info(f"Found {len(summary_files)} summary file(s)")

        # Validate each summary file
        for summary_file in summary_files:
            if not _validate_summary_file(summary_file, issues, warnings, verbose):
                issues.append(f"Invalid summary file format: {summary_file.name}")

    # Check for compressed result archives (*.tar.gz)
    archive_files = list(output_dir.glob("*.tar.gz"))
    if not archive_files and summary_files:
        warnings.append(
            "No compressed archives (.tar.gz) found - intermediate files may not be archived"
        )

    # Check for the combined summary file (gvclass_summary.tsv)
    combined_summary = output_dir / "gvclass_summary.tsv"
    if combined_summary.exists():
        if verbose:
            logger.info("Found combined summary file: gvclass_summary.tsv")
        if not _validate_combined_summary(combined_summary, issues, warnings, verbose):
            issues.append("Invalid combined summary file format")
    else:
        # It's okay if combined summary doesn't exist yet - it might be created after validation
        if verbose:
            logger.info(
                "Combined summary file not yet created (will be generated after validation)"
            )

    # Check for query subdirectories
    query_dirs = [d for d in output_dir.iterdir() if d.is_dir()]
    if query_dirs and verbose:
        logger.info(f"Found {len(query_dirs)} query subdirectory(ies)")
        for query_dir in query_dirs:
            _validate_query_directory(query_dir, issues, warnings, verbose)

    success = len(issues) == 0

    if verbose:
        if success:
            logger.info("✅ Validation successful")
        else:
            logger.warning(f"❌ Validation found {len(issues)} issue(s)")
        if warnings:
            logger.warning(f"⚠️  {len(warnings)} warning(s) found")

    return {"success": success, "issues": issues, "warnings": warnings}


def _validate_summary_file(
    summary_file: Path, issues: List[str], warnings: List[str], verbose: bool
) -> bool:
    """
    Validate individual summary.tab file format.

    Returns True if file is valid, False otherwise.
    """
    try:
        with open(summary_file, "r") as f:
            lines = f.readlines()

        if len(lines) < 2:
            warnings.append(
                f"{summary_file.name} has fewer than 2 lines (header + data)"
            )
            return False

        # Check header line contains expected columns
        header = lines[0].strip().split("\t")
        required_columns = ["query", "taxonomy_majority", "taxonomy_strict"]

        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            issues.append(
                f"{summary_file.name} missing required columns: {', '.join(missing_columns)}"
            )
            return False

        if verbose:
            logger.debug(f"✓ {summary_file.name} has valid format")

        return True

    except Exception as e:
        issues.append(f"Error reading {summary_file.name}: {str(e)}")
        return False


def _validate_combined_summary(
    combined_summary: Path, issues: List[str], warnings: List[str], verbose: bool
) -> bool:
    """
    Validate combined gvclass_summary.tsv file.

    Returns True if file is valid, False otherwise.
    """
    try:
        with open(combined_summary, "r") as f:
            lines = f.readlines()

        if len(lines) < 1:
            issues.append("Combined summary file is empty")
            return False

        # Check header
        header = lines[0].strip().split("\t")
        required_columns = ["query", "taxonomy_majority", "taxonomy_strict"]

        missing_columns = [col for col in required_columns if col not in header]
        if missing_columns:
            issues.append(
                f"Combined summary missing required columns: {', '.join(missing_columns)}"
            )
            return False

        # Warn if no data rows
        if len(lines) == 1:
            warnings.append("Combined summary has header but no data rows")

        if verbose:
            logger.debug(f"✓ Combined summary has {len(lines)-1} data row(s)")

        return True

    except Exception as e:
        issues.append(f"Error reading combined summary: {str(e)}")
        return False


def _validate_query_directory(
    query_dir: Path, issues: List[str], warnings: List[str], verbose: bool
) -> None:
    """
    Validate contents of individual query subdirectory.
    """
    query_name = query_dir.name

    # Check for key output files
    expected_files = [f"{query_name}.summary.tab", f"{query_name}.tar.gz"]

    for expected_file in expected_files:
        file_path = query_dir.parent / expected_file
        if not file_path.exists():
            # Check if file is in the subdirectory instead
            alt_path = query_dir / expected_file
            if not alt_path.exists():
                warnings.append(
                    f"Expected file not found for {query_name}: {expected_file}"
                )

    # Check for marker results (if not in fast mode)
    marker_files = list(query_dir.glob("*.faa"))
    tree_files = list(query_dir.glob("*.tree"))

    if marker_files and not tree_files:
        warnings.append(f"Query {query_name} has alignment files but no tree files")

    if verbose and (marker_files or tree_files):
        logger.debug(
            f"  {query_name}: {len(marker_files)} alignments, {len(tree_files)} trees"
        )
