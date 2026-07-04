"""NCLDV neighbor selection for the species-tree supermatrix (Section 2).

For each query, gather the nearest NCLDV reference *genomes* per GVOG8 group via
a **dedicated NCLDV-only similarity search + per-group gene tree**, then
dereplicate to a reference genome set that seeds the supermatrix.

Why a dedicated search rather than reading the production gene trees: the shared
``queryrefs_genetrees/<group>.treefile`` are built from a mixed top-100 BLAST, so
for a sparse group like ``GVOGm0461`` (2,322 NCLDV / 18,298 total references) far
fewer than ``k`` NCLDV leaves survive. Searching a prefix-filtered NCLDV-only
target with a raised per-query cap guarantees ``>= k`` reachable NCLDV neighbors
for every GVOG8 group. Crucially this path is **isolated** from the production
gene-tree path: it writes only into the species-tree scratch dir and never
touches ``queryrefs_genetrees`` or ``taxonomy_majority``, so the golden output is
unchanged whether or not ``--species-tree`` is on.

The pure tree-based selection primitive (:func:`select_ncldv_neighbors_from_trees`)
operates on any newick trees and is unit tested with synthetic trees; the
search/tree builders reuse the shared FAMSA/VeryFastTree helpers from
:mod:`src.core.alignment`.
"""

from __future__ import annotations

import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from Bio import SeqIO

from src.core.alignment import align_sequences_pyfamsa, run_veryfasttree
from src.core.blast import run_blastp
from src.core.tree_analysis import TreeAnalyzer

logger = logging.getLogger(__name__)

NCLDV_PREFIX = "NCLDV__"
PrefixFilter = Union[str, Sequence[str]]


def _prefix_tuple(keep_prefix: PrefixFilter) -> Tuple[str, ...]:
    if isinstance(keep_prefix, str):
        return (keep_prefix,)
    return tuple(keep_prefix)


def _prefix_label(keep_prefix: PrefixFilter) -> str:
    return ",".join(_prefix_tuple(keep_prefix))


def _prefix_slug(keep_prefix: PrefixFilter) -> str:
    parts = []
    for prefix in _prefix_tuple(keep_prefix):
        clean = prefix.strip("_").replace("-", "_").lower()
        parts.append(re.sub(r"[^a-z0-9_]+", "_", clean).strip("_"))
    return "_".join(part for part in parts if part) or "refs"


@dataclass(frozen=True)
class NeighborHit:
    """A single placed reference neighbor for one query, in one marker group.

    ``protein_id`` is the *specific* reference protein placed nearest the query in
    that group's tree (Section 3 preserves it as the genome's representative);
    ``genome_id`` is its resolved genome; ``distance`` is the patristic distance to
    the nearest query leaf.
    """

    genome_id: str
    protein_id: str
    group: str
    distance: float


# ---------------------------------------------------------------------------
# Pure tree-based selection primitives (unit tested with synthetic trees).
# ---------------------------------------------------------------------------


def select_ncldv_neighbors_from_trees(
    group_to_tree: Dict[str, Path],
    analyzer: TreeAnalyzer,
    query_id: str,
    k: int = 50,
    keep_prefix: PrefixFilter = NCLDV_PREFIX,
) -> List[NeighborHit]:
    """Top-``k`` reference neighbors per group, flattened to one list of hits.

    Args:
        group_to_tree: Mapping of GVOG8 group name -> its dedicated newick tree.
        analyzer: A :class:`TreeAnalyzer` (provides ``get_top_k_neighbors`` and the
            genome-id resolution against the taxonomy labels).
        query_id: Query genome identifier (query leaves start with it).
        k: Max neighbor genomes per group.
        keep_prefix: Reference-prefix gate applied to leaf names. A tuple allows
            auxiliary references, for example ``("NCLDV__", "EUK-pEVE__")``.

    Returns:
        A flat list of :class:`NeighborHit`; groups whose tree is missing/empty
        simply contribute nothing.
    """
    hits: List[NeighborHit] = []
    for group, tree_file in group_to_tree.items():
        for protein_id, genome_id, distance in analyzer.get_top_k_neighbors(
            tree_file, query_id, k=k, keep_prefix=keep_prefix
        ):
            hits.append(
                NeighborHit(
                    genome_id=genome_id,
                    protein_id=protein_id,
                    group=group,
                    distance=distance,
                )
            )
    return hits


# ---------------------------------------------------------------------------
# Dedicated NCLDV-only search + per-group gene-tree builders (the breadth
# mechanism). These reuse the shared alignment/tree helpers from Section 1.
# ---------------------------------------------------------------------------


def build_ncldv_reference_subset(
    group_faa: Path, out_faa: Path, prefix: PrefixFilter = NCLDV_PREFIX
) -> int:
    """Write the ``prefix``-filtered subset of a group reference faa.

    This is the dedicated panel-specific search *target*; building it once per
    group (and sharing it across queries in combined mode) is what makes ``>= k``
    eligible neighbors reachable for sparse groups.

    Returns the number of sequences written.
    """
    # Stream line-by-line rather than SeqIO.parse: the group faa can be tens of MB
    # (grp_COG0086 is ~49 MB / 63k records) and per-record SeqIO object creation
    # dominates the per-query wall clock. A header's id-token carries the domain
    # prefix, so ``startswith(">" + prefix)`` is exactly the id-token test. Written
    # atomically (tmp -> os.replace) so one subset file can be shared across the
    # parallel per-query workers without a torn read.
    tokens = tuple(">" + p for p in _prefix_tuple(prefix))
    count = 0
    tmp = out_faa.with_name(out_faa.name + ".tmp")
    with open(group_faa) as in_handle, open(tmp, "w") as out_handle:
        keep = False
        for line in in_handle:
            if line.startswith(">"):
                keep = line.startswith(tokens)
                if keep:
                    count += 1
            if keep:
                out_handle.write(line)
    os.replace(tmp, out_faa)
    return count


def gather_ncldv_neighbor_proteins(
    query_faa: Path,
    ncldv_ref_faa: Path,
    blast_out: Path,
    max_candidates: int = 300,
    threads: int = 4,
    min_identity: float = 30.0,
    overfetch: int = 1000,
) -> List[str]:
    """pyswrd-search query proteins vs the NCLDV-only references; return up to
    ``max_candidates`` candidate reference protein IDs for this group.

    The result is a deterministic **per-group** candidate set: unique subjects are
    ranked by ``(-score, subject_id)`` and truncated to ``max_candidates``. This
    bounds the dedicated gene tree to ~``max_candidates`` leaves regardless of how
    many query paralogs hit the group (a per-query cap would let ``N`` paralogs
    inflate the pool to ``max_candidates * N``), and makes equal-score ties
    resolve deterministically rather than depending on pyswrd's thread-sensitive
    emission order. The ranking is deliberately local (not the production
    ``parse_blastp``) so this route stays both deterministic and isolated from the
    gene-tree path.

    ``run_blastp`` is asked to emit ``overfetch`` hits per query (well above
    ``max_candidates``) so its own tie-sensitive per-query top-N cut sits far from
    the deterministic ``max_candidates`` cut applied here — making the kept set
    reproducible even when many references share a score.
    """
    run_blastp(
        str(query_faa),
        str(ncldv_ref_faa),
        str(blast_out),
        threads=threads,
        top_per_query=max(overfetch, max_candidates),
    )
    return _rank_unique_subjects(blast_out, max_candidates, min_identity)


def _rank_unique_subjects(
    m8_path: Path, max_candidates: int, min_identity: float
) -> List[str]:
    """Deterministically rank unique subject IDs from a pyswrd m8 file.

    Keeps the best (highest) score per subject, filters by ``min_identity``, sorts
    by ``(-score, subject_id)`` and truncates to ``max_candidates``.
    """
    if not m8_path.exists():
        return []

    best_score: Dict[str, float] = {}
    with open(m8_path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            try:
                subject = parts[1]
                identity = float(parts[2])
                score = float(parts[11])
            except (ValueError, IndexError):
                continue
            if identity < min_identity:
                continue
            if subject not in best_score or score > best_score[subject]:
                best_score[subject] = score

    ranked = sorted(best_score.items(), key=lambda item: (-item[1], item[0]))
    return [subject for subject, _ in ranked[:max_candidates]]


def build_species_tree_gene_tree(
    query_faa: Path,
    candidate_protein_ids: Iterable[str],
    ncldv_ref_faa: Path,
    out_tree: Path,
    threads: int = 4,
    min_sequences: int = 3,
) -> Optional[Path]:
    """Build a dedicated per-group gene tree from query proteins + candidate NCLDV
    reference proteins, reusing the shared FAMSA + VeryFastTree helpers.

    Writes the alignment alongside ``out_tree`` (``<stem>.aln.fasta``) and the
    newick to ``out_tree``. Returns ``out_tree``, or ``None`` if fewer than
    ``min_sequences`` total sequences are available (too few to place).
    """
    candidate_set = set(candidate_protein_ids)

    records = list(SeqIO.parse(str(query_faa), "fasta"))
    for record in SeqIO.parse(str(ncldv_ref_faa), "fasta"):
        # ``candidate_set`` holds id-tokens from the m8 (``record.id`` form).
        if record.id in candidate_set:
            records.append(record)

    if len(records) < min_sequences:
        logger.info(
            "Skipping species-tree gene tree %s: only %d sequences (< %d)",
            out_tree.name,
            len(records),
            min_sequences,
        )
        return None

    aligned = align_sequences_pyfamsa(records, threads)
    aln_path = out_tree.with_suffix(".aln.fasta")
    with open(aln_path, "w") as aln_handle:
        for seq_id, seq in aligned:
            # Replace stop codons (*) with gaps before tree inference, mirroring
            # the production gene-tree path: real NCLDV references carry * (e.g.
            # GVOGm0461 has 252/2322 such records) which VeryFastTree mishandles.
            aln_handle.write(f">{seq_id}\n{seq.replace('*', '-')}\n")

    # Multi-threaded for speed. VeryFastTree is nondeterministic at threads>=2, so
    # the exact selected reference SET (the top-k boundary) can vary run-to-run;
    # the user-facing placement taxonomy is robust to this (verified: identical
    # nearest-reference lineages across runs with different reference sets). Use
    # --threads 1 if a byte-identical reference set is required.
    newick = run_veryfasttree(aln_path, threads)
    out_tree.write_text(newick)
    return out_tree


# ---------------------------------------------------------------------------
# Per-query orchestration (composes the above; called by the Section 6 hook).
# ---------------------------------------------------------------------------


def build_query_group_trees(
    query_hits_dir: Path,
    ref_faa_dir: Path,
    groups: Sequence[str],
    query_id: str,
    work_dir: Path,
    max_candidates: int = 300,
    threads: int = 4,
    keep_prefix: PrefixFilter = NCLDV_PREFIX,
    subset_cache: Optional[Dict[Tuple[str, Tuple[str, ...]], Path]] = None,
) -> Dict[str, Path]:
    """Build the dedicated per-group gene trees for one query.

    For each group with query hits: build (or reuse, via ``subset_cache``) the
    prefix-filtered reference subset, gather candidate neighbor proteins by
    similarity, and build a dedicated gene tree. Groups with no query hits, no
    reference faa, no eligible references, no candidates, or too few sequences
    are skipped.

    These trees are the shared substrate for **both** neighbor selection (Section
    2) and query-representative selection (Section 3), so the Section 6 hook builds
    them once and feeds them to both.

    Returns:
        ``{group -> tree_path}`` for the groups that produced a tree.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    group_to_tree: Dict[str, Path] = {}

    for group in groups:
        query_faa = query_hits_dir / f"{group}.faa"
        if not query_faa.exists() or query_faa.stat().st_size == 0:
            continue

        group_ref = ref_faa_dir / f"{group}.faa"
        if not group_ref.exists():
            logger.warning("Reference faa missing for group %s: %s", group, group_ref)
            continue

        ncldv_subset = _resolve_ncldv_subset(
            group, group_ref, work_dir, keep_prefix, subset_cache
        )
        if ncldv_subset is None:
            continue

        candidates = gather_ncldv_neighbor_proteins(
            query_faa,
            ncldv_subset,
            work_dir / f"{query_id}.{group}.m8",
            max_candidates=max_candidates,
            threads=threads,
        )
        if not candidates:
            continue

        tree_path = build_species_tree_gene_tree(
            query_faa,
            candidates,
            ncldv_subset,
            work_dir / f"{query_id}.{group}.treefile",
            threads=threads,
        )
        if tree_path is not None:
            group_to_tree[group] = tree_path

    return group_to_tree


def _resolve_ncldv_subset(
    group: str,
    group_ref: Path,
    work_dir: Path,
    keep_prefix: PrefixFilter,
    subset_cache: Optional[Dict[Tuple[str, Tuple[str, ...]], Path]],
) -> Optional[Path]:
    """Return the prefix-filtered subset faa for ``group``, building it once and reusing
    a cached path across queries (combined mode) when ``subset_cache`` is given."""
    prefixes = _prefix_tuple(keep_prefix)
    cache_key = (group, prefixes)
    if subset_cache is not None and cache_key in subset_cache:
        return subset_cache[cache_key]

    ncldv_subset = work_dir / f"{group}.{_prefix_slug(prefixes)}.faa"
    if not ncldv_subset.exists():
        n_refs = build_ncldv_reference_subset(
            group_ref, ncldv_subset, prefix=keep_prefix
        )
        if n_refs == 0:
            logger.warning("No %s references in %s", _prefix_label(prefixes), group_ref)
            return None

    if subset_cache is not None:
        subset_cache[cache_key] = ncldv_subset
    return ncldv_subset
