"""
File and path helpers for the GVClass summary.

Stateless counterparts of the ``FullSummarizer`` methods that only read query
output files; the class keeps thin delegating wrappers for every name here.
"""

import logging
from pathlib import Path
from typing import Any, Dict, Set

import pandas as pd

from src.core.marker_extraction import (
    parse_hmm_output,
    parse_hmm_scores,
)

logger = logging.getLogger(__name__)


def _resolve_query_fna_path(query_id: str, query_output_dir: Path) -> Path:
    candidates = [
        query_output_dir / "query_fna" / f"{query_id}.fna",
        query_output_dir / f"{query_id}_reformatted.fna",
        query_output_dir
        / "gene_calling"
        / f"{query_id}_reformatted"
        / f"{query_id}_reformatted.fna",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return Path()


def _resolve_query_faa_path(query_id: str, query_output_dir: Path) -> Path:
    candidates = [
        query_output_dir / "query_faa" / f"{query_id}.faa",
        query_output_dir
        / "gene_calling"
        / f"{query_id}_reformatted"
        / f"{query_id}_reformatted.faa",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(query_output_dir.glob("gene_calling/**/*.faa"))
    return matches[0] if matches else Path()


def _read_stats_tsv(stats_tsv_file: Path) -> Dict[str, Any]:
    """Read base stats from the reformat stage TSV."""
    basic_stats = {}
    if not stats_tsv_file.exists():
        return basic_stats

    try:
        stats_df = pd.read_csv(stats_tsv_file, sep="\t")
        if not stats_df.empty:
            basic_stats = stats_df.iloc[0].to_dict()
    except Exception as e:
        logger.error(f"Error reading stats TSV: {e}")
    return basic_stats


def _read_stats_tab(stats_tab_file: Path, basic_stats: Dict[str, Any]) -> None:
    """Merge coding table/genecount fields from legacy stats.tab file."""
    if not stats_tab_file.exists():
        return

    try:
        with open(stats_tab_file, "r") as f:
            for line in f:
                _update_basic_stats_from_tab_line(line, basic_stats)
    except Exception as e:
        logger.error(f"Error reading stats tab: {e}")


def _update_basic_stats_from_tab_line(line: str, basic_stats: Dict[str, Any]) -> None:
    if line.startswith("ttable\t"):
        basic_stats["ttable"] = line.strip().split("\t")[1]
        return
    if line.startswith("genes\t"):
        basic_stats["genecount"] = int(line.strip().split("\t")[1])
        return
    if line.startswith("coding_density\t"):
        basic_stats["CODINGperc"] = float(line.strip().split("\t")[1])


def _load_basic_stats(query_id: str, query_output_dir: Path) -> Dict[str, Any]:
    """Load query summary statistics from available stats files."""
    stats_tsv_file = query_output_dir / "stats" / f"{query_id}_stats.tsv"
    stats_tab_file = query_output_dir / "stats" / f"{query_id}.stats.tab"

    basic_stats = _read_stats_tsv(stats_tsv_file)
    _read_stats_tab(stats_tab_file, basic_stats)
    return basic_stats


def _initialize_summary_result(
    query_id: str, basic_stats: Dict[str, Any]
) -> Dict[str, Any]:
    """Create the base result dictionary with normalized defaults."""
    ttable = basic_stats.get("ttable", "unknown")
    return {
        "query": query_id,
        "contigs": basic_stats.get("contigs", 0),
        "LENbp": basic_stats.get("LENbp", 0),
        "GCperc": basic_stats.get("GCperc", 0.0),
        "genecount": basic_stats.get("genecount", 0),
        "CODINGperc": basic_stats.get("CODINGperc", 0.0),
        "ttable": "codemeta" if ttable == "0" else ttable,
    }


def _load_marker_counts(counts_file: Path) -> Dict[str, int]:
    """Load marker counts from a models.counts file."""
    marker_counts = {}
    if not counts_file.exists():
        return marker_counts

    try:
        with open(counts_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) != 2:
                    continue
                try:
                    marker_counts[parts[0]] = int(parts[1])
                except ValueError:
                    logger.warning(
                        "Skipping malformed marker count '%s' in %s",
                        line.strip(),
                        counts_file,
                    )
    except Exception as e:
        logger.error(f"Error reading marker counts from {counts_file}: {e}")
    return marker_counts


def _load_marker_hits(hmmout_file: Path) -> Dict[str, Set[str]]:
    """Load per-marker protein hits from a filtered HMM output file."""
    if not hmmout_file.exists():
        return {}

    try:
        return parse_hmm_output(hmmout_file)
    except Exception as exc:
        logger.error("Error reading marker hits from %s: %s", hmmout_file, exc)
        return {}


def _load_marker_scores(hmmout_file: Path) -> Dict[str, Dict[str, float]]:
    """Load per-protein per-marker bit scores from a filtered HMM output file.

    Used to consolidate the MRYA panel against its marker groups (see
    :meth:`_add_marker_metrics`).
    """
    if not hmmout_file.exists():
        return {}

    try:
        return parse_hmm_scores(hmmout_file)
    except Exception as exc:
        logger.error("Error reading marker scores from %s: %s", hmmout_file, exc)
        return {}


def _candidate_type_for_reason(reason: str) -> str:
    mapping = {
        "cellular_markers": "cellular",
        "cellular_hits": "cellular",
        "phage_markers": "phage",
        "phage_hits": "phage",
        "viral_mixture": "mixed_viral",
    }
    return mapping.get(str(reason), "uncertain")
