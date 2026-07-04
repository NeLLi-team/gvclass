"""Representative-protein gathering for the species-tree supermatrix (Section 3).

Picks exactly one protein per (taxon, group) for the supermatrix, then drops taxa
present in fewer than ``min_markers`` of the GVOG8 groups:

* **Query representative** per (query, group): the query protein placed nearest any
  NCLDV reference in that group's dedicated gene tree (tie-break on equal patristic
  distance: length -> id). Reuses the Section 2 per-group trees.
* **Reference representative** per (genome, group): the *exact* protein placed as
  the neighbor (preserving the paralog Section 2 chose) when the genome was a
  neighbor in that group; otherwise the genome's *longest* protein in that group,
  gathered to raise the genome's marker count toward ``>= min_markers``.

The ``>= min_markers`` (default 3/8) filter is applied to both queries and
references *before* alignment, so sparse taxa never enter the supermatrix.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

from Bio import SeqIO
from ete3 import Tree

from src.core.species_tree.neighbor_selection import NCLDV_PREFIX, NeighborHit
from src.core.tree_analysis import TreeAnalyzer

logger = logging.getLogger(__name__)

PrefixFilter = Union[str, Sequence[str]]


def _prefix_tuple(keep_prefix: PrefixFilter) -> Tuple[str, ...]:
    if isinstance(keep_prefix, str):
        return (keep_prefix,)
    return tuple(keep_prefix)


class ReferenceProteinIndex:
    """Random-access index over one group's reference faa.

    Provides genome -> proteins resolution (via :meth:`TreeAnalyzer._extract_genome_id`,
    so it matches neighbor selection's genome ids exactly), longest-protein
    selection, and lazy sequence fetch. Backed by ``Bio.SeqIO.index`` (offset
    index; ~0.4 s and ~50 MB for a 63k-record group faa).
    """

    def __init__(self, group_faa: Path, analyzer: TreeAnalyzer):
        self.group_faa = Path(group_faa)
        self._index = SeqIO.index(str(self.group_faa), "fasta")
        try:
            self._genome_to_proteins: Dict[str, List[str]] = defaultdict(list)
            for protein_id in self._index.keys():
                genome_id = analyzer._extract_genome_id(protein_id)
                self._genome_to_proteins[genome_id].append(protein_id)
        except Exception:
            # Don't leak the file handle if genome resolution fails mid-build.
            self._index.close()
            raise

    def __enter__(self) -> "ReferenceProteinIndex":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:
        # GC backstop; callers should still close() explicitly (or use `with` /
        # contextlib.ExitStack) so handles are released deterministically.
        try:
            self.close()
        except Exception:
            pass

    def proteins_for(self, genome_id: str) -> List[str]:
        return list(self._genome_to_proteins.get(genome_id, []))

    def sequence(self, protein_id: str) -> str:
        return str(self._index[protein_id].seq)

    def longest_protein(self, genome_id: str) -> Optional[str]:
        """The genome's longest protein in this group; tie -> lexicographic id."""
        proteins = self._genome_to_proteins.get(genome_id)
        if not proteins:
            return None
        return min(proteins, key=lambda pid: (-len(self._index[pid].seq), pid))

    def __contains__(self, protein_id: str) -> bool:
        return protein_id in self._index

    def close(self) -> None:
        index = getattr(self, "_index", None)
        if index is not None:
            index.close()
            self._index = None


def query_protein_lengths(query_faa: Path) -> Dict[str, int]:
    """Map ``protein_id -> length`` for a query group faa (query-rep tie-break)."""
    lengths: Dict[str, int] = {}
    if not query_faa.exists():
        return lengths
    for record in SeqIO.parse(str(query_faa), "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths


def select_query_representative_from_tree(
    tree_file: Path,
    query_id: str,
    keep_prefix: PrefixFilter = NCLDV_PREFIX,
    query_lengths: Optional[Dict[str, int]] = None,
) -> Optional[Tuple[str, float]]:
    """Return ``(query_protein_id, min_distance)`` for the query protein placed
    nearest any allowed reference leaf in the tree, or ``None``.

    Primary key is the minimum patristic distance to any reference leaf; equal
    distances (rare for real branch lengths) break on length (desc) then id (asc),
    so selection is deterministic.
    """
    if not tree_file.exists() or tree_file.with_suffix(".SKIPPED").exists():
        return None
    try:
        tree = Tree(str(tree_file))
    except Exception as exc:  # pragma: no cover - malformed newick guard
        logger.error("Error loading tree %s: %s", tree_file, exc)
        return None

    query_nodes = [
        leaf for leaf in tree.iter_leaves() if leaf.name.startswith(query_id)
    ]
    prefixes = _prefix_tuple(keep_prefix)
    ref_nodes = [
        leaf
        for leaf in tree.iter_leaves()
        if leaf.name.startswith(prefixes) and not leaf.name.startswith(query_id)
    ]
    if not query_nodes or not ref_nodes:
        return None

    lengths = query_lengths or {}

    best_key: Optional[Tuple[float, int, str]] = None
    best_protein: Optional[str] = None
    best_distance = float("inf")
    for query_node in query_nodes:
        min_dist = float("inf")
        for ref_node in ref_nodes:
            try:
                dist = tree.get_distance(query_node, ref_node)
            except Exception:
                continue
            if dist < min_dist:
                min_dist = dist
        if min_dist == float("inf"):
            continue
        key = (
            min_dist,
            -lengths.get(query_node.name, 0),
            query_node.name,
        )
        if best_key is None or key < best_key:
            best_key = key
            best_protein = query_node.name
            best_distance = min_dist

    if best_protein is None:
        return None
    return best_protein, best_distance


def gather_query_representatives(
    group_to_tree: Dict[str, Path],
    query_id: str,
    query_hits_dir: Optional[Path] = None,
    min_markers: int = 3,
    keep_prefix: PrefixFilter = NCLDV_PREFIX,
) -> Optional[Dict[str, str]]:
    """Per-group query representative proteins, or ``None`` if ``< min_markers``.

    Returns ``{group -> query_protein_id}`` over the groups whose dedicated tree
    placed the query, provided the query reaches ``min_markers`` of them.
    """
    reps: Dict[str, str] = {}
    for group, tree_file in group_to_tree.items():
        lengths = (
            query_protein_lengths(query_hits_dir / f"{group}.faa")
            if query_hits_dir is not None
            else None
        )
        chosen = select_query_representative_from_tree(
            tree_file,
            query_id,
            keep_prefix=keep_prefix,
            query_lengths=lengths,
        )
        if chosen is not None:
            reps[group] = chosen[0]

    if len(reps) < min_markers:
        logger.info(
            "Query %s has only %d/%d species-tree markers (< %d) — excluded",
            query_id,
            len(reps),
            len(group_to_tree),
            min_markers,
        )
        return None
    return reps


def gather_reference_representatives(
    neighbor_hits: Sequence[NeighborHit],
    ref_indexes: Dict[str, ReferenceProteinIndex],
    groups: Sequence[str],
    min_markers: int = 3,
) -> Tuple[Dict[str, Dict[str, str]], List[str]]:
    """Per-(genome, group) reference representative proteins + dropped genomes.

    For each genome that was a neighbor in any group: keep the placed neighbor
    protein for those groups, then fill *additional* groups with the genome's
    longest protein (when present) to raise its marker count. Genomes reaching
    ``>= min_markers`` groups are kept; the rest are returned as dropped.

    Returns ``({genome -> {group -> protein_id}}, dropped_genomes)``.

    When the same ``(genome, group)`` is placed by several queries (combined
    mode), the **nearest** placement wins (tie-break on protein id), so the
    representative is deterministic regardless of the order queries were
    processed in. Hits for groups outside ``groups`` are ignored.
    """
    group_set = set(groups)

    # Nearest placed paralog per (genome, group): (distance, protein_id) order is
    # independent of hit ordering across queries.
    best_hit: Dict[Tuple[str, str], Tuple[float, str]] = {}
    for hit in neighbor_hits:
        if hit.group not in group_set:
            continue
        key = (hit.genome_id, hit.group)
        candidate = (hit.distance, hit.protein_id)
        if key not in best_hit or candidate < best_hit[key]:
            best_hit[key] = candidate

    placed: Dict[str, Dict[str, str]] = defaultdict(dict)
    for (genome_id, group), (_distance, protein_id) in best_hit.items():
        placed[genome_id][group] = protein_id

    kept: Dict[str, Dict[str, str]] = {}
    dropped: List[str] = []
    for genome_id, group_protein in placed.items():
        reps = dict(group_protein)
        for group in groups:
            if group in reps:
                continue
            index = ref_indexes.get(group)
            if index is None:
                continue
            longest = index.longest_protein(genome_id)
            if longest is not None:
                reps[group] = longest
        if len(reps) >= min_markers:
            kept[genome_id] = reps
        else:
            dropped.append(genome_id)

    if dropped:
        logger.info(
            "Dropped %d reference genomes with < %d species-tree markers",
            len(dropped),
            min_markers,
        )
    return kept, dropped
