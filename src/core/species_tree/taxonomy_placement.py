"""Tree-placement taxonomy for the species tree (Section 5).

Assigns a root-invariant placement taxonomy to every query leaf of the
concatenated species tree:

* The tree is rooted on a deterministically chosen **reference** leaf. Because
  query-only clades exclude references, a reference root makes every maximal
  query-only clade a standard rooted monophyletic subtree, so clade detection is
  root-independent in result. Patristic distances are mathematically
  root-invariant (bit-level results can differ by <=1e-5 across rootings from
  floating-point summation order, so distance tie-breaks round to 9 dp to keep
  the chosen nearest reference rooting-invariant).
* Patristic distances use a precompute: one preorder pass for root->node depth
  sums plus cached leaf ancestor sets, then ``dist(a,b) = depth[a] + depth[b]
  - 2*depth[lca]`` (O(path length) per pair; parity-tested vs ``tree.get_distance``).
* Each query gets exactly one row; references never produce rows; ids are
  de-namespaced for output. Queries sharing a query-only clade (>= 2 members)
  share ``species_tree_clade_id`` and the clade's nearest-reference taxonomy.
"""

from __future__ import annotations

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from ete3 import Tree, TreeNode

from src.core.species_tree.supermatrix import denamespace, is_query_leaf
from src.core.tree_analysis import TreeAnalyzer

logger = logging.getLogger(__name__)

_RANK_PREFIXES = ["d_", "p_", "c_", "o_", "f_", "g_", "s_"]

SPECIES_TREE_TAXONOMY_COLUMNS = [
    "query",
    "species_tree_nn_genome",
    "species_tree_nn_taxonomy",
    "species_tree_nn_distance",
    "species_tree_clade_id",
    "species_tree_clade_members",
    "species_tree_clade_mean_dist",
    "species_tree_clade_ref_dist",
]


def format_lineage(labels_str: str) -> str:
    """Render a ``DOMAIN|phylum|...|species`` labels string as gvclass-style
    ``d_DOMAIN;p_phylum;...`` (matching ``taxonomy_majority`` formatting)."""
    parts = labels_str.split("|")
    fields = [
        f"{_RANK_PREFIXES[i]}{part}"
        for i, part in enumerate(parts)
        if i < len(_RANK_PREFIXES) and part
    ]
    return ";".join(fields)


@dataclass
class PlacementRow:
    query: str
    nn_genome: str
    nn_taxonomy: str
    nn_distance: float
    clade_id: str = ""
    clade_members: str = ""
    clade_mean_dist: Optional[float] = None
    clade_ref_dist: Optional[float] = None

    def summary_dict(self) -> Dict[str, object]:
        """The four placement columns merged into a query's ``summary_data``."""
        return {
            "species_tree_nn_genome": self.nn_genome,
            "species_tree_nn_taxonomy": self.nn_taxonomy,
            "species_tree_nn_distance": self.nn_distance,
            "species_tree_clade_id": self.clade_id,
        }


class _Patristic:
    """Patristic distances over a fixed (rooted) tree via depth + LCA precompute."""

    def __init__(self, tree: Tree):
        self._depth: Dict[TreeNode, float] = {}
        for node in tree.traverse("preorder"):
            self._depth[node] = (
                0.0 if node.up is None else self._depth[node.up] + node.dist
            )
        self._ancestor_cache: Dict[TreeNode, Set[TreeNode]] = {}

    def _ancestors(self, leaf: TreeNode) -> Set[TreeNode]:
        cached = self._ancestor_cache.get(leaf)
        if cached is None:
            cached = set()
            node: Optional[TreeNode] = leaf
            while node is not None:
                cached.add(node)
                node = node.up
            self._ancestor_cache[leaf] = cached
        return cached

    def dist(self, a: TreeNode, b: TreeNode) -> float:
        if a is b:
            return 0.0
        ancestors_a = self._ancestors(a)
        node: Optional[TreeNode] = b
        while node is not None and node not in ancestors_a:
            node = node.up
        lca_depth = self._depth[node] if node is not None else 0.0
        return self._depth[a] + self._depth[b] - 2.0 * lca_depth


class SpeciesTreeClassifier:
    """Assign placement taxonomy to query leaves of a species tree."""

    def __init__(self, labels_file: Path):
        self.analyzer = TreeAnalyzer(labels_file)

    def _taxonomy(self, genome_id: str) -> str:
        labels_str = self.analyzer.labels_dict.get(genome_id)
        if labels_str is None:
            labels_str = self.analyzer.labels_dict.get(
                self.analyzer._extract_genome_id(genome_id)
            )
        if labels_str is None:
            logger.warning(
                "No taxonomy label for species-tree reference %s; using 'unknown'",
                genome_id,
            )
            return "unknown"
        return format_lineage(labels_str)

    def classify(self, tree_file: Path) -> List[PlacementRow]:
        """Return one :class:`PlacementRow` per query leaf, sorted by query.

        Every query yields exactly one row **when the tree has at least one
        reference leaf**. A tree with no reference leaf (or no query leaf) admits
        no placement and returns ``[]`` (logged); combined mode then writes no
        ``species_tree_taxonomy.tsv``. Reference leaves are domain-gated upstream
        (NCLDV for the default route; PPV/MIRUS via the Section 7 config), so each
        lineage is taken from ``labels.tsv`` as-is and a missing label degrades to
        ``"unknown"`` rather than failing the run.
        """
        try:
            tree = Tree(str(tree_file))
        except Exception as exc:  # pragma: no cover - malformed newick guard
            logger.error("Error loading species tree %s: %s", tree_file, exc)
            return []

        query_leaves = [leaf for leaf in tree.get_leaves() if is_query_leaf(leaf.name)]
        ref_leaves = [
            leaf for leaf in tree.get_leaves() if not is_query_leaf(leaf.name)
        ]
        if not query_leaves or not ref_leaves:
            logger.warning(
                "Species tree %s has %d query / %d reference leaves; no placement",
                tree_file,
                len(query_leaves),
                len(ref_leaves),
            )
            return []

        # Root on a deterministic reference leaf so query-only clades are rooted
        # monophyletic subtrees (root-independent in result).
        outgroup = min(ref_leaves, key=lambda leaf: leaf.name)
        try:
            tree.set_outgroup(outgroup)
        except Exception as exc:
            # Fail closed: detecting clades on an arbitrary rooting would break the
            # root-invariance guarantee, so emit no placement rather than a wrong
            # one. VeryFastTree output is unrooted, so this is a pathological guard.
            logger.warning(
                "set_outgroup failed on %s (%s); skipping species-tree placement for %s",
                outgroup.name,
                exc,
                tree_file,
            )
            return []

        # Re-fetch node objects after rooting and rebuild distance precompute.
        query_leaves = [leaf for leaf in tree.get_leaves() if is_query_leaf(leaf.name)]
        ref_leaves = [leaf for leaf in tree.get_leaves() if not is_query_leaf(leaf.name)]
        patristic = _Patristic(tree)

        clade_root_of = self._maximal_query_clades(tree)
        clades: Dict[TreeNode, List[TreeNode]] = defaultdict(list)
        for query in query_leaves:
            clades[clade_root_of[query]].append(query)

        rows: List[PlacementRow] = []
        for members in clades.values():
            if len(members) >= 2:
                rows.extend(self._clade_rows(members, ref_leaves, patristic))
            else:
                rows.append(self._solo_row(members[0], ref_leaves, patristic))

        return sorted(rows, key=lambda row: row.query)

    @staticmethod
    def _maximal_query_clades(tree: Tree) -> Dict[TreeNode, TreeNode]:
        """Map each query leaf to its maximal all-query subtree root."""
        all_query: Dict[TreeNode, bool] = {}
        for node in tree.traverse("postorder"):
            if node.is_leaf():
                all_query[node] = is_query_leaf(node.name)
            else:
                all_query[node] = all(all_query[child] for child in node.children)

        clade_root_of: Dict[TreeNode, TreeNode] = {}
        for node in tree.traverse("preorder"):
            if all_query.get(node) and (
                node.up is None or not all_query.get(node.up)
            ):
                for leaf in node.get_leaves():
                    clade_root_of[leaf] = node
        return clade_root_of

    def _nearest_reference(
        self, leaf: TreeNode, ref_leaves: List[TreeNode], patristic: _Patristic
    ) -> Tuple[TreeNode, float]:
        """Reference leaf at minimum patristic distance.

        The tie-break key rounds distance to 9 dp (below the 6-dp output precision,
        above the ~1e-5 re-rooting drift) so the chosen reference — and thus
        nn_genome/nn_taxonomy — is invariant to the input tree's rooting even when
        two references are near-equidistant. The raw distance is returned.
        """
        best_ref = None
        best_key: Optional[Tuple[float, str]] = None
        best_dist = 0.0
        for ref in ref_leaves:
            dist = patristic.dist(leaf, ref)
            key = (round(dist, 9), ref.name)
            if best_key is None or key < best_key:
                best_key = key
                best_ref = ref
                best_dist = dist
        return best_ref, best_dist

    def _solo_row(
        self, query: TreeNode, ref_leaves: List[TreeNode], patristic: _Patristic
    ) -> PlacementRow:
        nn_ref, nn_dist = self._nearest_reference(query, ref_leaves, patristic)
        return PlacementRow(
            query=denamespace(query.name),
            nn_genome=nn_ref.name,
            nn_taxonomy=self._taxonomy(nn_ref.name),
            nn_distance=nn_dist,
        )

    def _clade_rows(
        self,
        members: List[TreeNode],
        ref_leaves: List[TreeNode],
        patristic: _Patristic,
    ) -> List[PlacementRow]:
        member_names = sorted(denamespace(member.name) for member in members)
        clade_id = member_names[0]
        clade_members = ";".join(member_names)
        clade_mean = _mean_pairwise(members, patristic)

        # Clade -> nearest external reference: min over member x reference. The
        # tie-break key rounds distance to 9 dp so the choice is rooting-invariant.
        best_ref = None
        best_key: Optional[Tuple[float, str]] = None
        clade_ref_dist = 0.0
        for member in members:
            for ref in ref_leaves:
                dist = patristic.dist(member, ref)
                key = (round(dist, 9), ref.name)
                if best_key is None or key < best_key:
                    best_key = key
                    best_ref = ref
                    clade_ref_dist = dist
        nn_genome = best_ref.name
        nn_taxonomy = self._taxonomy(nn_genome)

        rows = []
        for member in members:
            rows.append(
                PlacementRow(
                    query=denamespace(member.name),
                    nn_genome=nn_genome,
                    nn_taxonomy=nn_taxonomy,
                    nn_distance=patristic.dist(member, best_ref),
                    clade_id=clade_id,
                    clade_members=clade_members,
                    clade_mean_dist=clade_mean,
                    clade_ref_dist=clade_ref_dist,
                )
            )
        return rows


def _mean_pairwise(members: List[TreeNode], patristic: _Patristic) -> float:
    total = 0.0
    pairs = 0
    for i in range(len(members)):
        for j in range(i + 1, len(members)):
            total += patristic.dist(members[i], members[j])
            pairs += 1
    return total / pairs if pairs else 0.0


def _fmt(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.6f}"


def write_species_tree_taxonomy(rows: List[PlacementRow], out_path: Path) -> Path:
    """Write the ``species_tree_taxonomy.tsv`` table."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(SPECIES_TREE_TAXONOMY_COLUMNS)
        for row in rows:
            writer.writerow(
                [
                    row.query,
                    row.nn_genome,
                    row.nn_taxonomy,
                    _fmt(row.nn_distance),
                    row.clade_id,
                    row.clade_members,
                    _fmt(row.clade_mean_dist),
                    _fmt(row.clade_ref_dist),
                ]
            )
    return out_path
