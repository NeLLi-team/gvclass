"""File and path helpers for the GVClass summary.

Stateless helpers that only read query output files.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from src.core.marker_extraction import parse_hmm_scores

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class QueryEvidence:
    """Marker evidence of one query, parsed once from its ``hmmout`` files.

    Attributes:
        counts: Proteins per marker model (``models.counts``).
        hits: Proteins that hit each marker model (``models.out.filtered``).
        scores: Best bit score of each protein per marker model.
    """

    counts: dict[str, int]
    hits: dict[str, set[str]]
    scores: dict[str, dict[str, float]]


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


def _read_stats_tsv(stats_tsv_file: Path) -> dict[str, Any]:
    """Read base stats from the reformat stage TSV."""
    if not stats_tsv_file.exists():
        return {}
    stats_df = pd.read_csv(stats_tsv_file, sep="\t")
    if stats_df.empty:
        return {}
    return stats_df.iloc[0].to_dict()


def _read_stats_tab(stats_tab_file: Path, basic_stats: dict[str, Any]) -> None:
    """Merge coding table/genecount fields from legacy stats.tab file."""
    if not stats_tab_file.exists():
        return
    with stats_tab_file.open(encoding="utf-8") as handle:
        for line in handle:
            _update_basic_stats_from_tab_line(line, basic_stats)


def _update_basic_stats_from_tab_line(line: str, basic_stats: dict[str, Any]) -> None:
    if line.startswith("ttable\t"):
        basic_stats["ttable"] = line.strip().split("\t")[1]
        return
    if line.startswith("genes\t"):
        basic_stats["genecount"] = int(line.strip().split("\t")[1])
        return
    if line.startswith("coding_density\t"):
        basic_stats["CODINGperc"] = float(line.strip().split("\t")[1])


def _load_basic_stats(query_id: str, query_output_dir: Path) -> dict[str, Any]:
    """Load query summary statistics from available stats files."""
    stats_tsv_file = query_output_dir / "stats" / f"{query_id}_stats.tsv"
    stats_tab_file = query_output_dir / "stats" / f"{query_id}.stats.tab"

    basic_stats = _read_stats_tsv(stats_tsv_file)
    _read_stats_tab(stats_tab_file, basic_stats)
    return basic_stats


def _initialize_summary_result(
    query_id: str, basic_stats: dict[str, Any]
) -> dict[str, Any]:
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


def load_query_evidence(query_output_dir: Path) -> QueryEvidence:
    """Parse the marker tables of a query once.

    A missing table means the query has no marker evidence. The hits derive from
    the scores, so ``models.out.filtered`` is read a single time.
    """
    hmmout_dir = query_output_dir / "hmmout"
    counts = _load_marker_counts(hmmout_dir / "models.counts")
    filtered_table = hmmout_dir / "models.out.filtered"
    scores = parse_hmm_scores(filtered_table) if filtered_table.exists() else {}
    hits: dict[str, set[str]] = {}
    for protein, model_scores in scores.items():
        for model in model_scores:
            hits.setdefault(model, set()).add(protein)
    return QueryEvidence(counts=counts, hits=hits, scores=scores)


def _load_marker_counts(counts_file: Path) -> dict[str, int]:
    """Load marker counts from a models.counts file."""
    marker_counts: dict[str, int] = {}
    if not counts_file.exists():
        return marker_counts

    with counts_file.open(encoding="utf-8") as handle:
        for line in handle:
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
    return marker_counts


def _candidate_type_for_reason(reason: str) -> str:
    mapping = {
        "cellular_coherent_lineage": "cellular",
        "cellular_markers": "cellular",
        "cellular_hits": "cellular",
        "phage_markers": "phage",
        "phage_hits": "phage",
        "viral_mixture": "mixed_viral",
    }
    return mapping.get(str(reason), "uncertain")
