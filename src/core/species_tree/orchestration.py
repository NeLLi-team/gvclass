"""Species-tree orchestration: per-query hook + opt-in combined step.

* Per-query (inside ``process_single_query``, after summary data exists and before
  the query dir is archived/deleted): the DEFAULT product. Gated on a registered
  panel matching ``taxonomy_majority``; build the dedicated per-group trees, select
  neighbors + the query's representative proteins, write a run-level sidecar, then
  build this query's own supermatrix/tree/placement (``k=NEIGHBORS_PER_QUERY_TREE``)
  under ``<out>/species_tree/<query>/`` and merge the placement columns into its
  ``summary_data`` (so the per-query placement owns the summary columns and
  persists in the per-query summary for ``--resume``).
* Combined (inside ``gvclass_flow`` after parallel processing): OPT-IN via
  ``--species-tree-combined``. Read this run's sidecars, build one
  supermatrix/tree/placement per panel (``k=NEIGHBORS_PER_COMBINED_TREE``), and
  write ``<out>/species_tree/`` as an additional artifact. It does NOT touch the
  summary columns — the per-query placement is the single source of truth there.

Both entry points are wrapped so a species-tree failure is logged and never blocks
the standard pipeline/summary.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO

from src.core.species_tree.config import (
    NEIGHBORS_PER_COMBINED_TREE,
    NEIGHBORS_PER_QUERY_TREE,
    PANELS,
    SpeciesTreePanel,
    neighbors_to_store,
    select_panel,
)
from src.core.species_tree.handoff import (
    read_sidecars,
    scratch_dir_for,
    write_sidecar,
)
from src.core.species_tree.neighbor_selection import (
    NeighborHit,
    build_query_group_trees,
    select_ncldv_neighbors_from_trees,
)
from src.core.species_tree.protein_gathering import (
    ReferenceProteinIndex,
    gather_query_representatives,
    gather_reference_representatives,
)
from src.core.species_tree.supermatrix import build_supermatrix, namespace_query
from src.core.species_tree.taxonomy_placement import (
    SpeciesTreeClassifier,
    write_species_tree_taxonomy,
)
from src.core.species_tree.tree_inference import infer_species_tree
from src.core.tree_analysis import TreeAnalyzer

logger = logging.getLogger(__name__)


def _ref_faa_dir(database_path: Path) -> Path:
    for candidate in (database_path / "database" / "faa", database_path / "faa"):
        if candidate.exists():
            return candidate
    return database_path / "database" / "faa"


def _fetch_sequence(faa: Path, protein_id: str) -> Optional[str]:
    if not faa.exists():
        return None
    for record in SeqIO.parse(str(faa), "fasta"):
        if record.id == protein_id:
            return str(record.seq)
    return None


# ---------------------------------------------------------------------------
# Per-query hook
# ---------------------------------------------------------------------------


def run_per_query_species_tree(
    query_id: str,
    query_output_dir: Path,
    database_path: Path,
    summary_data: Dict,
    output_base: Path,
    threads: int,
    trim_method: str = "witchi",
) -> None:
    """Per-query species-tree step (gated on the genome's domain panel).

    This is the default species-tree product: it builds the query's own tree and
    merges the placement into ``summary_data``.
    """
    try:
        panel = select_panel(summary_data.get("taxonomy_majority", ""))
        if panel is None:
            logger.debug("Species-tree: %s not a registered domain; no-op", query_id)
            return
        _run_per_query(
            panel,
            query_id,
            query_output_dir,
            database_path,
            summary_data,
            output_base,
            threads,
            trim_method,
        )
    except Exception as exc:  # never block the standard pipeline
        logger.warning("Species-tree per-query step failed for %s: %s", query_id, exc)


def _run_per_query(
    panel: SpeciesTreePanel,
    query_id: str,
    query_output_dir: Path,
    database_path: Path,
    summary_data: Dict,
    output_base: Path,
    threads: int,
    trim_method: str = "witchi",
) -> None:
    query_hits_dir = query_output_dir / "query_hits_faa"
    if not query_hits_dir.exists():
        return

    ref_faa_dir = _ref_faa_dir(database_path)
    analyzer = TreeAnalyzer(database_path / "labels.tsv")
    scratch = scratch_dir_for(output_base)
    work = scratch / query_id / "work"

    group_to_tree = build_query_group_trees(
        query_hits_dir,
        ref_faa_dir,
        panel.groups,
        query_id,
        work,
        threads=threads,
        keep_prefix=panel.domain_prefix,
    )
    if not group_to_tree:
        logger.info("No species-tree gene trees built for %s", query_id)
        return

    neighbors = select_ncldv_neighbors_from_trees(
        group_to_tree,
        analyzer,
        query_id,
        k=neighbors_to_store(),
        keep_prefix=panel.domain_prefix,
    )
    query_reps = gather_query_representatives(
        group_to_tree,
        query_id,
        query_hits_dir,
        min_markers=panel.min_markers,
        keep_prefix=panel.domain_prefix,
    )
    if query_reps is None:
        logger.info(
            "Query %s has < %d species-tree markers; excluded", query_id, panel.min_markers
        )
        return

    query_rep_seqs: Dict[str, str] = {}
    for group, protein_id in query_reps.items():
        seq = _fetch_sequence(query_hits_dir / f"{group}.faa", protein_id)
        if seq is not None:
            query_rep_seqs[group] = seq

    payload = {
        "query_id": query_id,
        "panel": panel.name,
        "query_reps": query_rep_seqs,
        "neighbors": [
            {
                "genome_id": hit.genome_id,
                "protein_id": hit.protein_id,
                "group": hit.group,
                "distance": hit.distance,
            }
            for hit in neighbors
        ],
    }
    write_sidecar(scratch, query_id, payload)

    # Build this query's own species tree (the default product) and merge the
    # placement into summary_data, which is written to the per-query summary just
    # after this hook returns. The per-query placement OWNS the summary columns;
    # the opt-in combined tree is an additional artifact that never overwrites
    # them. Output lands in the run-level <out>/species_tree/<query>/ (which
    # survives the per-query tar+rmtree, unlike the query dir) so it is the
    # discoverable default deliverable.
    _build_and_place(
        panel,
        {query_id: payload},
        database_path,
        analyzer,
        work / "query_supermatrix",
        output_base / "species_tree" / query_id,
        threads,
        summary_by_query={query_id: summary_data},
        trim_method=trim_method,
        k_neighbors=NEIGHBORS_PER_QUERY_TREE,
        basename=query_id,
    )


# ---------------------------------------------------------------------------
# Combined step
# ---------------------------------------------------------------------------


def run_combined_species_tree(
    output_base: Path,
    database_path: Path,
    results: List[Dict],
    threads: int,
    trim_method: str = "witchi",
) -> None:
    """Combined species-tree step over this run's sidecars (one tree per panel)."""
    try:
        _run_combined(output_base, database_path, results, threads, trim_method)
    except Exception as exc:  # never block the final summary
        logger.warning("Combined species-tree step failed: %s", exc)


def _run_combined(
    output_base: Path,
    database_path: Path,
    results: List[Dict],
    threads: int,
    trim_method: str = "witchi",
) -> None:
    scratch = scratch_dir_for(output_base)
    # Only place queries whose standard pipeline COMPLETED — a query can write its
    # sidecar+.done and then fail in _post_process_query; such a query is reported
    # failed and must not appear as a leaf in the combined tree / taxonomy table.
    completed = {res["query"] for res in results if res.get("status") == "complete"}
    sidecars = {
        query_id: sidecar
        for query_id, sidecar in read_sidecars(scratch).items()
        if query_id in completed
    }
    if not sidecars:
        logger.info("No completed species-tree sidecars for this run; skipping combined tree")
        return

    by_panel: Dict[str, Dict[str, Dict]] = defaultdict(dict)
    for query_id, sidecar in sidecars.items():
        panel_name = sidecar.get("panel")
        if panel_name not in PANELS:
            logger.warning(
                "Species-tree sidecar %s has unknown panel %r; skipping", query_id, panel_name
            )
            continue
        by_panel[panel_name][query_id] = sidecar

    analyzer = TreeAnalyzer(database_path / "labels.tsv")

    multi = len(by_panel) > 1
    for panel_name, panel_sidecars in by_panel.items():
        panel = PANELS.get(panel_name)
        if panel is None:
            continue
        base_out = output_base / "species_tree"
        # Single panel (the common case): combined.* at the species_tree/ root.
        # Multiple panels: namespace under species_tree/_combined/<panel>/ so the
        # combined dirs can never collide with a per-query subdir whose genome
        # stem happens to equal a panel name (NCLDV/PPV/MIRUS).
        out_dir = base_out / "_combined" / panel_name if multi else base_out
        # The combined tree is an additional artifact only: pass an empty
        # summary map so it never overwrites the per-query placement columns
        # (which the per-query hook already wrote to each query's summary, so
        # they also survive --resume without any rewrite here).
        _build_and_place(
            panel,
            panel_sidecars,
            database_path,
            analyzer,
            scratch / "combined_work" / panel_name,
            out_dir,
            threads,
            summary_by_query={},
            trim_method=trim_method,
            k_neighbors=NEIGHBORS_PER_COMBINED_TREE,
        )


# ---------------------------------------------------------------------------
# Shared assembly: sidecars -> supermatrix -> tree -> placement -> merge
# ---------------------------------------------------------------------------


def _truncate_neighbors(sidecars: Dict[str, Dict], k_neighbors: int) -> List[NeighborHit]:
    """Flatten sidecar neighbors, keeping each query's top-k per group by distance.

    Sidecars store up to ``neighbors_to_store()`` neighbors so both the per-query
    and combined trees can be served from the same data; this truncates to the
    requested ``k_neighbors`` per (query, group) deterministically (nearest first,
    genome-id tie-break) before the reference set is assembled.
    """
    hits: List[NeighborHit] = []
    for sidecar in sidecars.values():
        by_group: Dict[str, List[Dict]] = defaultdict(list)
        for hit in sidecar.get("neighbors", []):
            by_group[hit["group"]].append(hit)
        for group_hits in by_group.values():
            group_hits.sort(key=lambda h: (h["distance"], h["genome_id"]))
            for hit in group_hits[:k_neighbors]:
                hits.append(
                    NeighborHit(
                        genome_id=hit["genome_id"],
                        protein_id=hit["protein_id"],
                        group=hit["group"],
                        distance=hit["distance"],
                    )
                )
    return hits


def _build_and_place(
    panel: SpeciesTreePanel,
    sidecars: Dict[str, Dict],
    database_path: Path,
    analyzer: TreeAnalyzer,
    work_dir: Path,
    out_dir: Path,
    threads: int,
    summary_by_query: Dict[str, Dict],
    trim_method: str = "witchi",
    k_neighbors: int = NEIGHBORS_PER_COMBINED_TREE,
    basename: str = "combined",
) -> None:
    ref_faa_dir = _ref_faa_dir(database_path)
    all_neighbors = _truncate_neighbors(sidecars, k_neighbors)

    with ExitStack() as stack:
        ref_indexes: Dict[str, ReferenceProteinIndex] = {}
        for group in panel.groups:
            faa = ref_faa_dir / f"{group}.faa"
            if faa.exists():
                ref_indexes[group] = stack.enter_context(
                    ReferenceProteinIndex(faa, analyzer)
                )

        ref_reps, dropped = gather_reference_representatives(
            all_neighbors, ref_indexes, panel.groups, min_markers=panel.min_markers
        )
        logger.info(
            "Species-tree (%s): %d reference genomes (%d dropped < %d markers), %d queries",
            panel.name,
            len(ref_reps),
            len(dropped),
            panel.min_markers,
            len(sidecars),
        )

        taxa: Dict[str, Dict[str, str]] = {}
        for query_id, sidecar in sidecars.items():
            reps = {g: s for g, s in sidecar.get("query_reps", {}).items() if s}
            if len(reps) >= panel.min_markers:
                taxa[namespace_query(query_id)] = reps
        for genome_id, group_protein in ref_reps.items():
            taxa[genome_id] = {
                group: ref_indexes[group].sequence(protein_id)
                for group, protein_id in group_protein.items()
                if group in ref_indexes
            }

        result = build_supermatrix(
            taxa,
            panel.groups,
            work_dir,
            out_dir,
            basename=basename,
            threads=threads,
            min_taxa_per_group=panel.min_markers,
            trim_method=trim_method,
        )

    if result is None:
        logger.warning("Species-tree supermatrix empty for panel %s", panel.name)
        return

    tree = infer_species_tree(result.supermatrix_faa, out_dir / f"{basename}.treefile", threads)
    classifier = SpeciesTreeClassifier(database_path / "labels.tsv")
    rows = classifier.classify(tree)
    write_species_tree_taxonomy(rows, out_dir / "species_tree_taxonomy.tsv")

    merged = 0
    for row in rows:
        summary_data = summary_by_query.get(row.query)
        if summary_data is not None:
            summary_data.update(row.summary_dict())
            merged += 1
    logger.info(
        "Species-tree (%s): %d placements written to %s (%d merged into summary)",
        panel.name,
        len(rows),
        out_dir,
        merged,
    )
