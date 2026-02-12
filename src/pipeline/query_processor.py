"""Compatibility wrapper for query-level GVClass processing."""

from pathlib import Path
from typing import List, Dict, Any

from src.pipeline.query_processing_engine import process_single_query


def run_query_processing(
    query_file: Path,
    output_base: Path,
    database_path: Path,
    genetic_codes: List[int],
    tree_method: str,
    mode_fast: bool,
    sensitive_mode: bool = False,
    threads: int = 4,
) -> Dict[str, Any]:
    """Process a single query through the entire pipeline."""
    return process_single_query(
        query_file=query_file,
        output_base=output_base,
        database_path=database_path,
        genetic_codes=genetic_codes,
        tree_method=tree_method,
        mode_fast=mode_fast,
        sensitive_mode=sensitive_mode,
        threads=threads,
    )
