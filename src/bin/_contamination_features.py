"""Shared feature recomputation for the contamination benchmark/retraining scripts.

The contamination model's training features (``CONTAMINATION_MODEL_FEATURES``) are
no longer all present in ``gvclass_summary.tsv``: the final-summary overhaul
dropped the v1 rule-based signals, and the marker-panel standardization dropped the
cellular/phage marker counts. Scraping the summary table therefore silently zeros
those features, which would degrade any retrained model.

This helper instead recomputes the full feature dict directly from a query's
per-query output directory (``hmmout/``, ``stats/``, ``blastp_out/``, ``query_*``)
using the runtime :class:`~src.core.summarize_full.FullSummarizer` — the single
source of truth that the live pipeline uses to feed the model.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional

from src.core.summarize_full import FullSummarizer


def load_tree_nn_results(
    tree_nn_file: Path,
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Reconstruct ``tree_nn_results`` (marker -> query_protein -> {neighbor: dist})
    from a ``stats/<query>.tree_nn`` TSV written by the runtime."""
    results: Dict[str, Dict[str, Dict[str, float]]] = defaultdict(dict)
    if not tree_nn_file.exists():
        return results
    with tree_nn_file.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "marker":
                continue
            marker, query_protein, neighbor, distance = (
                parts[0],
                parts[1],
                parts[2],
                parts[3],
            )
            try:
                dist = float(distance)
            except ValueError:
                continue
            results[marker].setdefault(query_protein, {})[neighbor] = dist
    return results


def compute_query_features(
    query_id: str,
    query_output_dir: Path,
    database_path: Path,
    summarizer: Optional[FullSummarizer] = None,
    sensitive_mode: bool = True,
) -> Dict[str, object]:
    """Recompute the full feature dict for one query from its per-query output dir.

    Returns the complete ``summary_data`` dict, a superset of
    ``CONTAMINATION_MODEL_FEATURES`` (v1 rule-based signals, ``order_dup`` /
    ``gvog8_dup``, the cellular/phage marker counts, and the per-contig
    ``suspicious_*`` / ``cellular_coherent_*`` features). Pass a shared
    ``summarizer`` to avoid reloading the database for every query.
    """
    if summarizer is None:
        summarizer = FullSummarizer(database_path, sensitive_mode=sensitive_mode)
    tree_nn = load_tree_nn_results(
        query_output_dir / "stats" / f"{query_id}.tree_nn"
    )
    return summarizer.summarize_query_full(query_id, query_output_dir, tree_nn)
