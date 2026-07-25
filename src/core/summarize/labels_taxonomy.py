"""
Reference labels and taxonomy counting for the GVClass summary.

Stateless counterparts of the ``FullSummarizer`` methods that read the
reference ``labels.tsv`` and turn tree nearest-neighbour placements into
taxonomy counters; the class keeps thin delegating wrappers for every name
here.

The label index these functions need (``labels_dict``) is passed in as an
explicit first argument, so nothing here reads instance state.
"""

from collections import Counter, defaultdict
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from src.config.marker_sets import (
    CAPS_MARKER_GROUP,
    PLV_PREFIX,
)
from src.core.marker_extraction import (
    count_unique_proteins_for_prefixes,
)
from src.core.summarize.consensus import TAX_LEVELS
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


def _dominant(counter: Counter) -> str:
    """Single dominant key, or "" on an empty counter or a top-2 tie."""
    top = counter.most_common(2)
    if not top or (len(top) > 1 and top[0][1] == top[1][1]):
        return ""
    return top[0][0]


def load_labels(labels_file: Path) -> Dict[str, List[str]]:
    """Load taxonomy labels from file."""
    labels = {}
    try:
        with open(labels_file, "r") as f:
            for line in f:
                parsed = _parse_label_line(line)
                if parsed is None:
                    continue
                genome_id, tax_parts = parsed
                labels[genome_id] = tax_parts
    except Exception as e:
        logger.error(f"Error loading labels: {e}")
    return labels


def _parse_label_line(line: str) -> Optional[Tuple[str, List[str]]]:
    if line.startswith("#"):
        return None

    parts = line.strip().split("\t")
    if len(parts) < 2:
        return None

    genome_id = parts[0]
    tax_parts = parts[1].split("|") if "|" in parts[1] else [parts[1]]
    while len(tax_parts) < 7:
        tax_parts.append("")
    return genome_id, tax_parts


def _extract_genome_id(labels_dict: Dict[str, List[str]], protein_id: str) -> str:
    """
    Extract genome ID from a protein ID for labels lookup.

    Protein IDs can have formats like:
    - PPV__IMGVR_UViG_3300044959|000235_9 -> PPV__IMGVR_UViG_3300044959|000235
    - NCLDV__GCA_000123456|contig_1_42 -> NCLDV__GCA_000123456|contig_1

    The protein suffix is typically _N where N is a number at the end.

    Args:
        protein_id: Full protein identifier from tree

    Returns:
        Genome/contig ID suitable for labels lookup
    """
    # First try: check if the full ID (minus trailing _digits) is in labels
    # This handles: PPV__IMGVR_UViG_3300044959|000235_9 -> PPV__IMGVR_UViG_3300044959|000235
    stripped = re.sub(r"_\d+$", "", protein_id)
    if stripped in labels_dict:
        return stripped

    # Second try: just the part before | (for older format entries)
    # This handles: genome_id|protein_id -> genome_id
    if "|" in protein_id:
        base_id = protein_id.split("|")[0]
        if base_id in labels_dict:
            return base_id

    # Third try: strip protein suffix from full ID even if not in labels
    return stripped


def _collect_taxonomy_counts(
    labels_dict: Dict[str, List[str]],
    tree_nn_results: Dict[str, Dict[str, Dict[str, float]]],
) -> Tuple[
    Dict[str, Counter],
    Dict[str, Dict[str, Counter]],
    Dict[str, Counter[Tuple[str, ...]]],
    List[float],
]:
    """Aggregate taxonomy counts from NN tree hits.

    Returns a tuple of:
    * ``tax_counters`` — flat per-level Counter (legacy shape, used by
      contamination scoring and the per-level detailed columns).
    * ``per_marker_counters`` — ``Dict[level, Dict[marker, Counter]]``
      that retains the per-marker breakdown. The per-marker majority
      rule in :meth:`_build_consensus_taxonomies` uses this to avoid
      the prior bug where a single marker with many paralog hits could
      dominate the majority vote by sheer count.
    * ``per_marker_lineage_counters`` — ``Dict[marker, Counter[lineage]]``
      retaining the full 7-rank lineage for each hit. Consensus taxonomy
      uses this to constrain lower ranks to the already selected parent
      lineage.
    * ``distances`` — flat list of patristic distances for the mean.
    """
    tax_counters = {level: Counter() for level in TAX_LEVELS}
    per_marker_counters: Dict[str, Dict[str, Counter]] = {
        level: defaultdict(Counter) for level in TAX_LEVELS
    }
    per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]] = defaultdict(
        Counter
    )
    distances: List[float] = []

    # One physical query protein can hit several overlapping marker models
    # (e.g. a genome-packaging ATPase matching the NCLDV ATPase markers AND
    # the PPV A32 marker PLV_PC_054). Count each protein ONCE, for its single
    # nearest (minimum patristic distance) placement, so a shared core gene
    # cannot cast multiple or cross-domain votes. Deterministic tie-break on
    # (distance, marker, neighbor).
    placements_by_protein: Dict[str, List[Tuple[float, str, str]]] = defaultdict(list)
    for marker, query_neighbors in tree_nn_results.items():
        for query_protein, neighbors in query_neighbors.items():
            for neighbor, distance in neighbors.items():
                placements_by_protein[query_protein].append(
                    (distance, marker, neighbor)
                )

    for placements in placements_by_protein.values():
        resolved = []
        for distance, marker, neighbor in placements:
            genome_id = _extract_genome_id(labels_dict, neighbor)
            tax_info = labels_dict.get(genome_id)
            if tax_info:
                resolved.append((distance, marker, neighbor, genome_id, tax_info))
        if not resolved:
            continue
        distance, marker, neighbor, genome_id, tax_info = min(
            resolved, key=lambda r: (r[0], r[1], r[2])
        )
        _add_taxonomy_counts_for_neighbor(genome_id, tax_info, tax_counters)
        _add_taxonomy_counts_for_marker(
            genome_id, tax_info, per_marker_counters, marker
        )
        lineage = _build_taxonomy_lineage(genome_id, tax_info)
        if any(lineage):
            per_marker_lineage_counters[marker][lineage] += 1
        distances.append(distance)

    return tax_counters, per_marker_counters, per_marker_lineage_counters, distances


def _build_taxonomy_lineage(genome_id: str, tax_info: List[str]) -> Tuple[str, ...]:
    """Return a normalized 7-rank lineage tuple for one labeled genome."""
    domain_prefix = genome_id.split("__", 1)[0] if "__" in genome_id else ""
    lineage: List[str] = []
    for level in TAX_LEVELS:
        idx = TAX_LEVEL_MAPPING[level]
        tax_value = tax_info[idx].strip() if idx < len(tax_info) else ""
        if level == "domain":
            tax_value = domain_prefix or tax_value
        lineage.append(tax_value)
    return tuple(lineage)


def _add_taxonomy_counts_for_neighbor(
    genome_id: str, tax_info: List[str], tax_counters: Dict[str, Counter]
) -> None:
    domain_prefix = genome_id.split("__")[0] if "__" in genome_id else ""
    for level, idx in TAX_LEVEL_MAPPING.items():
        if idx >= len(tax_info):
            continue
        tax_value = tax_info[idx].strip()
        if not tax_value:
            continue
        tax_key = _build_taxonomy_key(level, domain_prefix, tax_value)
        tax_counters[level][tax_key] += 1


def _add_taxonomy_counts_for_marker(
    genome_id: str,
    tax_info: List[str],
    per_marker_counters: Dict[str, Dict[str, Counter]],
    marker: str,
) -> None:
    domain_prefix = genome_id.split("__")[0] if "__" in genome_id else ""
    for level, idx in TAX_LEVEL_MAPPING.items():
        if idx >= len(tax_info):
            continue
        tax_value = tax_info[idx].strip()
        if not tax_value:
            continue
        tax_key = _build_taxonomy_key(level, domain_prefix, tax_value)
        per_marker_counters[level][marker][tax_key] += 1


def _plv_count_tree_aware(
    labels_dict: Dict[str, List[str]],
    marker_hits: Optional[Dict[str, Set[str]]],
    tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]],
) -> int:
    """Count A32/PLV-marker proteins ONLY when they genuinely place with PPV
    references in the marker tree.

    ``PLV_PC_054`` (the sole PLV marker) is the A32 packaging ATPase shared
    by NCLDV and PPV, so a plain NCLDV A32 hits the HMM too and the raw HMM
    count over-flags every giant virus. With the marker tree available we
    keep only proteins whose nearest reference is PPV; we fall back to the
    raw HMM count only when no tree is available (counts-only resume).
    """
    if not tree_nn_results:
        return count_unique_proteins_for_prefixes(marker_hits, [PLV_PREFIX])
    plv_proteins: Set[str] = set()
    for marker, proteins in (marker_hits or {}).items():
        if not marker.startswith(PLV_PREFIX):
            continue
        nn = tree_nn_results.get(marker, {})
        for protein in proteins:
            neighbors = nn.get(protein)
            if not neighbors:
                continue
            nearest, _dist = min(neighbors.items(), key=lambda kv: kv[1])
            tax = labels_dict.get(_extract_genome_id(labels_dict, nearest))
            if tax and tax[0] == "PPV":
                plv_proteins.add(protein)
    return len(plv_proteins)


def _protein_nearest_domain(
    labels_dict: Dict[str, List[str]],
    tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]],
) -> Dict[str, str]:
    """Map each query protein to the domain (NCLDV/PPV/MIRUS/BAC/EUK/...) of
    its single nearest tree neighbour across all marker trees.

    One physical protein is placed in several group trees; the minimum
    patristic distance (deterministic tie-break on marker then neighbour)
    decides its true lineage. Used to bucket cross-group panel hits — a VP or
    Mirus marker on a protein that actually places with NCLDV is not counted
    toward that panel.
    """
    placements: Dict[str, List[Tuple[float, str, str]]] = defaultdict(list)
    for marker, query_neighbors in (tree_nn_results or {}).items():
        for protein, neighbors in query_neighbors.items():
            for neighbor, distance in neighbors.items():
                placements[protein].append((distance, marker, neighbor))
    domains: Dict[str, str] = {}
    for protein, placement_list in placements.items():
        for _distance, _marker, neighbor in sorted(placement_list):
            genome_id = _extract_genome_id(labels_dict, neighbor)
            if genome_id in labels_dict:
                domains[protein] = (
                    genome_id.split("__", 1)[0] if "__" in genome_id else ""
                )
                break
    return domains


def _caps_group_from_tree(labels_dict: Dict[str, List[str]], tree_nn_results) -> str:
    """Dominant Bellas&Sommaruga caps group from the caps MCP tree-NN placement.

    A caps reference lineage carries its Bellas group at the family rank; the
    CAPS_GROUPS whitelist gates out everything else (caps group names occur only on
    caps references — under PPV or NCLDV domains — never on a non-caps genome).
    Each query protein votes ONCE: its single nearest neighbour is read per MCP group
    tree, and a cross-reactive MCP that places in BOTH mcp_plv and mcp_ncldv is
    resolved to its CLOSEST caps-group placement so it cannot double-count (mirrors
    contamination_scoring._per_protein_dominant_lineage). Dominant group across
    proteins; a tie, or no caps-group neighbour, resolves to "".
    """
    fam_idx = TAX_LEVEL_MAPPING["family"]
    per_protein: Dict[str, Tuple[Optional[str], float]] = {}
    for marker in ("mcp_plv", "mcp_ncldv"):
        for protein, neighbors in tree_nn_results.get(marker, {}).items():
            if not neighbors:
                continue
            nearest, dist = min(neighbors.items(), key=lambda kv: kv[1])
            tax = labels_dict.get(_extract_genome_id(labels_dict, nearest))
            if not tax or len(tax) <= fam_idx:
                continue
            family = tax[fam_idx].strip()
            if family not in CAPS_GROUPS:
                continue
            prev = per_protein.get(protein)
            if prev is None or dist < prev[1]:
                per_protein[protein] = (family, dist)
            elif dist == prev[1] and family != prev[0]:
                # equal-distance placements in different caps groups -> ambiguous, drop
                per_protein[protein] = (None, dist)
    groups: Counter = Counter(
        fam for fam, _dist in per_protein.values() if fam is not None
    )
    return _dominant(groups)


def _capsid_group_counter(labels_dict: Dict[str, List[str]], tree_nn_results) -> str:
    """Unified capsid-type tally across all MCP panels.

    For each query protein placed in an MCP group tree (mcp_ncldv / mcp_mirus
    / mcp_plv), read its single nearest neighbour and map it to a capsid-type
    label: a Bellas&Sommaruga caps family when the neighbour's family rank is
    in CAPS_GROUPS, otherwise the neighbour's phylum (Nucleocytoviricota /
    Mirusviricota). One vote per physical protein — a protein cross-placed in
    several MCP trees is resolved to its closest placement; an equal-distance
    cross-type tie is dropped. Returns "label:count" pairs comma-joined,
    sorted by count desc then label (e.g. "Nucleocytoviricota:4,Gossevirus:1").
    """
    fam_idx = TAX_LEVEL_MAPPING["family"]
    phy_idx = TAX_LEVEL_MAPPING["phylum"]
    capsid_phyla = {"Nucleocytoviricota", "Mirusviricota"}
    per_protein: Dict[str, Tuple[Optional[str], float]] = {}
    for marker in ("mcp_ncldv", "mcp_mirus", "mcp_plv"):
        for protein, neighbors in tree_nn_results.get(marker, {}).items():
            if not neighbors:
                continue
            nearest, dist = min(neighbors.items(), key=lambda kv: kv[1])
            tax = labels_dict.get(_extract_genome_id(labels_dict, nearest))
            if not tax:
                continue
            label: Optional[str] = None
            if len(tax) > fam_idx and tax[fam_idx].strip() in CAPS_GROUPS:
                label = tax[fam_idx].strip()
            elif len(tax) > phy_idx and tax[phy_idx].strip() in capsid_phyla:
                label = tax[phy_idx].strip()
            if label is None:
                continue
            prev = per_protein.get(protein)
            if prev is None or dist < prev[1]:
                per_protein[protein] = (label, dist)
            elif dist == prev[1] and label != prev[0]:
                # equal-distance placements of different capsid types -> drop
                per_protein[protein] = (None, dist)
    counts: Counter = Counter(
        lbl for lbl, _dist in per_protein.values() if lbl is not None
    )
    return ",".join(
        f"{lbl}:{n}"
        for lbl, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    )
