"""Taxonomy counting for the GVClass summary.

Stateless helpers that turn tree nearest-neighbour placements into taxonomy
counters.

The label index these functions need (``labels``, a
:class:`~src.core.labels.LabelIndex`) is passed in as an explicit first argument,
so nothing here reads instance state.
"""

import logging
from collections import Counter, defaultdict
from collections.abc import Mapping

from src.config.marker_sets import (
    CAPS_MARKER_GROUP,
    PLV_PREFIX,
)
from src.core.labels import (
    LabelIndex,
    Placement,
    common_lineage,
    domain_prefix,
    nearest_placements,
)
from src.core.marker_extraction import (
    count_unique_proteins_for_prefixes,
)
from src.core.summarize.consensus import TAX_LEVELS, Lineage
from src.core.summarize.tax_format import _build_taxonomy_key

logger = logging.getLogger(__name__)

TAX_LEVEL_MAPPING = {
    "domain": 0,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
    "species": 6,
}

# Bellas&Sommaruga caps group names (the family rank of caps reference lineages).
# Used to recognise a caps-group nearest neighbour in the MCP trees and ignore generic
# (non-caps) refs whose family rank is some other taxon.
CAPS_GROUPS = frozenset(CAPS_MARKER_GROUP.values()) - {"unclassified"}
CAPSID_PHYLA = frozenset({"Nucleocytoviricota", "Mirusviricota"})
MCP_MARKERS = ("mcp_ncldv", "mcp_mirus", "mcp_plv")


def collect_taxonomy_counts(
    labels: LabelIndex,
    tree_nn_results: dict[str, dict[str, dict[str, float]]],
) -> tuple[dict[str, Counter], dict[str, Counter[Lineage]], list[float]]:
    """Aggregate taxonomy counts from NN tree hits.

    One physical query protein can hit several overlapping marker models (e.g.
    a genome-packaging ATPase matching the NCLDV ATPase markers AND the PPV A32
    marker PLV_PC_054). Each protein counts ONCE, for its nearest placement, so
    a shared core gene cannot cast multiple or cross-domain votes. The nearest
    placement is the set of labelled neighbours at the minimal patristic
    distance over all of the protein's markers. When several neighbours tie,
    the protein votes with the lineage they share and abstains at the ranks
    below it. The marker credited is the smallest marker name among them.

    Returns a tuple of:
    * ``tax_counters``: flat per-level Counter, used by contamination scoring
      and the per-level detailed columns. It counts the shared ranks only.
    * ``per_marker_lineage_counters``: ``{marker: Counter[lineage]}``, the
      seven-rank lineage of each protein credited to the marker. Consensus
      taxonomy uses it to constrain lower ranks to the selected parent lineage.
    * ``distances``: flat list of the nearest patristic distances, for the mean.
    """
    tax_counters: dict[str, Counter] = {level: Counter() for level in TAX_LEVELS}
    per_marker_lineage_counters: dict[str, Counter[Lineage]] = defaultdict(Counter)
    distances: list[float] = []
    nearest_by_protein = _nearest_placements_by_protein(labels, tree_nn_results)
    for nearest in nearest_by_protein.values():
        lineage = common_lineage([placement.lineage for placement in nearest])
        credited = nearest[0]
        _add_lineage_counts(credited.genome_id, lineage, tax_counters)
        if any(lineage):
            per_marker_lineage_counters[credited.marker][lineage] += 1
        distances.append(credited.distance)
    return tax_counters, per_marker_lineage_counters, distances


def _nearest_placements_by_protein(
    labels: LabelIndex,
    tree_nn_results: Mapping[str, Mapping[str, Mapping[str, float]]],
) -> dict[str, list[Placement]]:
    """Map each query protein to its labelled neighbours at the minimal distance.

    The minimum runs over every marker tree that placed the protein. Neighbours
    without a reference label are skipped. Each list is sorted, so no caller
    depends on the insertion order of the tied neighbours.
    """
    placements: dict[str, list[Placement]] = defaultdict(list)
    for marker, query_neighbors in tree_nn_results.items():
        for protein, neighbors in query_neighbors.items():
            placements[protein].extend(labels.placements(marker, neighbors))
    nearest_by_protein: dict[str, list[Placement]] = {}
    for protein, found in placements.items():
        nearest = nearest_placements(found)
        if nearest:
            nearest_by_protein[protein] = nearest
    return nearest_by_protein


def _add_lineage_counts(
    genome_id: str, lineage: Lineage, tax_counters: dict[str, Counter]
) -> None:
    """Count each non-empty rank of ``lineage`` under its domain-prefixed key."""
    prefix = domain_prefix(genome_id)
    for level, tax_value in zip(TAX_LEVELS, lineage, strict=True):
        if tax_value:
            tax_key = _build_taxonomy_key(level, prefix, tax_value)
            tax_counters[level][tax_key] += 1


def plv_count_tree_aware(
    labels: LabelIndex,
    marker_hits: dict[str, set[str]] | None,
    tree_nn_results: dict[str, dict[str, dict[str, float]]] | None,
) -> int:
    """Count A32/PLV-marker proteins that genuinely place with PPV references.

    ``PLV_PC_054`` (the sole PLV marker) is the A32 packaging ATPase shared
    by NCLDV and PPV, so a plain NCLDV A32 hits the HMM too and the raw HMM
    count over-flags every giant virus. With the marker tree available we
    keep only proteins whose nearest references are all PPV (a tie between
    PPV and another domain does not count); we fall back to the raw HMM count
    only when no tree is available (counts-only resume).
    """
    if not tree_nn_results:
        return count_unique_proteins_for_prefixes(marker_hits, [PLV_PREFIX])
    plv_proteins: set[str] = set()
    for marker, proteins in (marker_hits or {}).items():
        if not marker.startswith(PLV_PREFIX):
            continue
        nn = tree_nn_results.get(marker, {})
        for protein in proteins:
            nearest = labels.nearest_lineage(marker, nn.get(protein, {}))
            if nearest is not None and nearest[0][0] == "PPV":
                plv_proteins.add(protein)
    return len(plv_proteins)


def protein_nearest_domain(
    labels: LabelIndex,
    tree_nn_results: dict[str, dict[str, dict[str, float]]] | None,
) -> dict[str, str]:
    """Map each query protein to the domain of its nearest tree neighbours.

    The domain is NCLDV, PPV, MIRUS, BAC, EUK or similar. One physical protein
    is placed in several group trees; the minimum patristic distance over all
    of them decides its true lineage. Tied nearest neighbours of different
    domains leave the protein out of the map, as for an unplaced protein. Used
    to bucket cross-group panel hits — a VP or Mirus marker on a protein that
    actually places with NCLDV is not counted toward that panel.
    """
    nearest_by_protein = _nearest_placements_by_protein(labels, tree_nn_results or {})
    domains: dict[str, str] = {}
    for protein, nearest in nearest_by_protein.items():
        prefixes = {domain_prefix(placement.genome_id) for placement in nearest}
        if len(prefixes) == 1:
            domains[protein] = prefixes.pop()
    return domains


def capsid_group_counter(
    labels: LabelIndex, tree_nn_results: Mapping[str, Mapping[str, Mapping[str, float]]]
) -> str:
    """Unified capsid-type tally across all MCP panels.

    One vote per physical protein. The protein's nearest labelled neighbours over
    every MCP group tree (``MCP_MARKERS``) decide: tied neighbours, in one tree or
    several, count with the lineage they share. That lineage maps to a
    capsid-type label: a Bellas&Sommaruga caps family when the family rank is in
    CAPS_GROUPS, otherwise the phylum (Nucleocytoviricota / Mirusviricota). A
    protein whose nearest lineage has no such label casts no vote; a farther
    placement never replaces it.
    Returns "label:count" pairs comma-joined, sorted by count desc then label
    (e.g. "Nucleocytoviricota:4,Gossevirus:1").
    """
    mcp_results = {marker: tree_nn_results.get(marker, {}) for marker in MCP_MARKERS}
    counts: Counter = Counter()
    for nearest in _nearest_placements_by_protein(labels, mcp_results).values():
        label = _capsid_label(
            common_lineage([placement.lineage for placement in nearest])
        )
        if label is not None:
            counts[label] += 1
    return ",".join(
        f"{lbl}:{n}"
        for lbl, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    )


def _capsid_label(lineage: Lineage) -> str | None:
    """Return the caps family of ``lineage``, else its capsid phylum, else None."""
    family = lineage[TAX_LEVEL_MAPPING["family"]]
    if family in CAPS_GROUPS:
        return family
    phylum = lineage[TAX_LEVEL_MAPPING["phylum"]]
    return phylum if phylum in CAPSID_PHYLA else None
