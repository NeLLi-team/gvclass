"""
Marker panel metrics for the GVClass summary.

Stateless counterparts of the ``FullSummarizer`` methods that turn marker
counts, hits and per-protein scores into panel completeness and duplication
values; the class keeps thin delegating wrappers for every name here.

The values that need the reference labels (``_protein_nearest_domain``,
``_plv_count_tree_aware``, ``_capsid_group_counter``) are computed by the class
and passed in, so nothing here touches ``labels_dict``.
"""

import logging
from typing import Any, Dict, Optional, Set, Tuple

from src.config.marker_sets import (
    MIRUS_CATEGORY_MODELS,
    BUSCO_MODELS,
    PHAGE_MODELS,
    GVOG4M_MODELS,
    GVOG8M_MODELS,
    UNI56_MODELS,
    NCLDV_MCP_MODELS,
    MRYA_MODELS,
    VP_CATEGORY_PREFIXES,
    PLV_PREFIX,
)
from src.core.marker_extraction import (
    consolidate_grouped_panel,
    count_unique_proteins_for_markers,
)

logger = logging.getLogger(__name__)


def _panel_has_lineage_protein(
    proteins: Set[str],
    protein_domain: Optional[Dict[str, str]],
    expected_domain: str,
) -> bool:
    """Whether a category has supporting evidence after tree bucketing.

    Without ``protein_domain`` (no tree placements) any hit counts (raw
    behaviour). With it, a protein still supports the category unless it places
    with a *different* known domain; proteins placing with ``expected_domain``
    or left unplaced are kept, so genuine but un-placed members are not
    dropped.
    """
    if not proteins:
        return False
    if protein_domain is None:
        return True
    for protein in proteins:
        domain = protein_domain.get(protein)
        if domain is None or domain == expected_domain:
            return True
    return False


def calculate_vp_metrics(
    marker_counts: Dict[str, int],
    marker_hits: Optional[Dict[str, Set[str]]],
    protein_domain: Optional[Dict[str, str]],
    plv_tree_count: int,
) -> Tuple[str, int, int, float]:
    """Calculate VP (Virophage / Preplasmiviricota) completeness metrics.

    VP has 4 core categories: MCP, Penton, ATPase, Protease.
    Uses prefix matching (e.g., VP_MCP_* matches MCP category).

    When ``protein_domain`` (each query protein's nearest tree-neighbour
    domain) is provided, a category is credited only if at least one hitting
    protein does not place with a non-PPV lineage — i.e. an NCLDV core gene
    that cross-hits a VP marker (e.g. the shared packaging ATPase) and places
    with NCLDV refs is excluded. Without tree placements the raw category
    presence is used.

    ``plv_tree_count`` is the tree-aware PLV protein count from
    ``FullSummarizer._plv_count_tree_aware``; it is only consulted on the
    marker-hits path, where the raw count would over-flag giant viruses.

    Returns:
        Tuple of (completeness_str, vp_mcp_count, plv_count, vp_df)
        - completeness_str: "n/N" format (N = number of core categories)
        - vp_mcp_count: count of unique proteins hitting VP_MCP markers
        - plv_count: count of unique proteins hitting PLV markers
        - vp_df: duplication factor (total VP hits / 4)
    """
    if marker_hits:
        category_proteins: Dict[str, Set[str]] = {
            category: set() for category in VP_CATEGORY_PREFIXES
        }
        for marker, proteins in marker_hits.items():
            for category, prefix in VP_CATEGORY_PREFIXES.items():
                if marker.startswith(prefix):
                    category_proteins[category].update(proteins)
        categories_present = {
            category
            for category, proteins in category_proteins.items()
            if _panel_has_lineage_protein(proteins, protein_domain, "PPV")
        }
        vp_mcp_count = len(category_proteins["MCP"])
        total_vp_hits = sum(len(proteins) for proteins in category_proteins.values())
        plv_count = plv_tree_count
        completeness = len(categories_present)
        vp_df = total_vp_hits / 4.0 if total_vp_hits > 0 else 0.0
        comp = f"{completeness}/{len(VP_CATEGORY_PREFIXES)}"
        return comp, vp_mcp_count, plv_count, vp_df

    categories_present = set()
    vp_mcp_count = 0
    total_vp_hits = 0
    plv_count = 0
    for marker, count in marker_counts.items():
        if count <= 0:
            continue

        # Check VP categories by prefix
        for category, prefix in VP_CATEGORY_PREFIXES.items():
            if marker.startswith(prefix):
                categories_present.add(category)
                total_vp_hits += count
                if category == "MCP":
                    vp_mcp_count += count

        # Check PLV markers
        if marker.startswith(PLV_PREFIX):
            plv_count += count

    completeness = len(categories_present)
    vp_df = total_vp_hits / 4.0 if total_vp_hits > 0 else 0.0

    comp = f"{completeness}/{len(VP_CATEGORY_PREFIXES)}"
    return comp, vp_mcp_count, plv_count, vp_df


def calculate_mirus_completeness(
    marker_counts: Dict[str, int],
    marker_hits: Optional[Dict[str, Set[str]]],
    protein_domain: Optional[Dict[str, str]],
) -> Tuple[str, float]:
    """Calculate Mirusviricota completeness metrics.

    Mirus has 4 core categories:
    - MCP: Mirus_MCP, Mirus_JellyRoll
    - ATPase: Mirus_Terminase_ATPase, Mirus_Terminase_merged
    - Portal: Mirus_Portal
    - Triplex: Mirus_Triplex1, Mirus_Triplex2

    When ``protein_domain`` (each query protein's nearest tree-neighbour
    domain) is provided, a category is credited only if a hitting protein does
    not place with a non-MIRUS lineage, so an NCLDV protein that cross-hits a
    Mirus marker (e.g. the terminase ATPase) and places with NCLDV refs is
    excluded. Without tree placements the raw category presence is used.

    Returns:
        Tuple of (completeness_str, mirus_df)
        - completeness_str: "n/N" format (N = number of core categories)
        - mirus_df: duplication factor (total Mirus hits / 4)
    """
    if marker_hits:
        category_proteins: Dict[str, Set[str]] = {
            category: set() for category in MIRUS_CATEGORY_MODELS
        }
        for category, models in MIRUS_CATEGORY_MODELS.items():
            for model in models:
                category_proteins[category].update(marker_hits.get(model, set()))
        categories_present = {
            category
            for category, proteins in category_proteins.items()
            if _panel_has_lineage_protein(proteins, protein_domain, "MIRUS")
        }
        total_mirus_hits = sum(len(proteins) for proteins in category_proteins.values())
        completeness = len(categories_present)
        mirus_df = total_mirus_hits / 4.0 if total_mirus_hits > 0 else 0.0
        return f"{completeness}/{len(MIRUS_CATEGORY_MODELS)}", mirus_df

    categories_present = set()
    total_mirus_hits = 0
    for marker, count in marker_counts.items():
        if count <= 0:
            continue

        # Check each Mirus category
        for category, models in MIRUS_CATEGORY_MODELS.items():
            if marker in models:
                categories_present.add(category)
                total_mirus_hits += count

    completeness = len(categories_present)
    mirus_df = total_mirus_hits / 4.0 if total_mirus_hits > 0 else 0.0

    return f"{completeness}/{len(MIRUS_CATEGORY_MODELS)}", mirus_df


def _add_marker_metrics(
    marker_counts: Dict[str, int],
    marker_hits: Optional[Dict[str, Set[str]]],
    marker_scores: Optional[Dict[str, Dict[str, float]]],
    protein_domain: Optional[Dict[str, str]],
    plv_tree_count: int,
    capsid_group: str,
) -> Dict[str, Any]:
    """Build the marker-based metrics for the query.

    Each core marker panel is reported as ``{panel}_completeness`` ("n/N":
    distinct models present out of the panel size) and ``{panel}_dup`` (the
    duplication factor = total hits / distinct models present), mirroring the
    vp_/mirus_completeness convention. The cellular single-copy markers are
    split into the BUSCO and COG (UNI56) panels.

    The MRYA panel is consolidated against its marker groups: several MRYA
    HMMs (``VLTF3``/``ATPase``/``gamadvirusMCP``/``VLTF2``) share gene-family
    groups with GVOG/NCLDV markers, so a generic NCLDV core protein hits them
    as a secondary match. Using per-protein bit scores, a grouped MRYA model
    is credited only when it is the protein's best hit within its group;
    otherwise the protein is a higher-scoring GVOG/NCLDV gene and is not
    counted toward MRYA. Falls back to the raw panel when no scores are
    available (e.g. a counts-only resume).
    """
    result: Dict[str, Any] = {}

    def _panel(models):
        present = sum(1 for m in models if marker_counts.get(m, 0) > 0)
        total = sum(marker_counts.get(m, 0) for m in models)
        return present, (total / present if present > 0 else 0.0)

    gvog4_present, result["gvog4_dup"] = _panel(GVOG4M_MODELS)
    result["gvog4_completeness"] = f"{gvog4_present}/{len(GVOG4M_MODELS)}"

    gvog8_present, result["gvog8_dup"] = _panel(GVOG8M_MODELS)
    result["gvog8_completeness"] = f"{gvog8_present}/{len(GVOG8M_MODELS)}"

    busco_present, result["busco_dup"] = _panel(BUSCO_MODELS)
    result["busco_completeness"] = f"{busco_present}/{len(BUSCO_MODELS)}"

    if marker_scores:
        cog_counts = consolidate_grouped_panel(marker_scores, UNI56_MODELS)
        cog_present = len(cog_counts)
        cog_total = sum(cog_counts.values())
        result["cog_dup"] = cog_total / cog_present if cog_present else 0.0
    else:
        cog_present, result["cog_dup"] = _panel(UNI56_MODELS)
    result["cog_completeness"] = f"{cog_present}/{len(UNI56_MODELS)}"

    if marker_scores:
        mrya_counts = consolidate_grouped_panel(marker_scores, MRYA_MODELS)
        mrya_present = len(mrya_counts)
        mrya_total = sum(mrya_counts.values())
        result["mrya_dup"] = mrya_total / mrya_present if mrya_present else 0.0
    else:
        if marker_counts:
            logger.warning(
                "No per-protein marker scores available (missing/empty "
                "models.out.filtered); MRYA panel falls back to raw counts, "
                "which overlapping GVOG/NCLDV markers can inflate."
            )
        mrya_present, result["mrya_dup"] = _panel(MRYA_MODELS)
    result["mrya_completeness"] = f"{mrya_present}/{len(MRYA_MODELS)}"

    phage_present, result["phage_dup"] = _panel(PHAGE_MODELS)
    result["phage_completeness"] = f"{phage_present}/{len(PHAGE_MODELS)}"

    # Standalone NCLDV major-capsid-protein count (capsid typing is reported
    # in capsid_group). Protein-aware when marker hits are available.
    if marker_hits:
        result["ncldv_mcp_total"] = count_unique_proteins_for_markers(
            marker_hits, NCLDV_MCP_MODELS
        )
    else:
        result["ncldv_mcp_total"] = sum(
            marker_counts.get(m, 0) for m in NCLDV_MCP_MODELS
        )

    # Internal contamination-model features: consumed by the trained model's
    # feature vector (contamination_scoring.CONTAMINATION_MODEL_FEATURES), not
    # emitted as summary columns. cellular = COG (UNI56) + BUSCO; counts are
    # protein-aware when marker hits are available.
    cellular_models = UNI56_MODELS + BUSCO_MODELS
    result["cellular_unique"] = sum(
        1 for m in cellular_models if marker_counts.get(m, 0) > 0
    )
    result["phage_unique"] = phage_present
    if marker_hits:
        result["cellular_total"] = count_unique_proteins_for_markers(
            marker_hits, cellular_models
        )
        result["phage_total"] = count_unique_proteins_for_markers(
            marker_hits, PHAGE_MODELS
        )
    else:
        result["cellular_total"] = sum(marker_counts.get(m, 0) for m in cellular_models)
        result["phage_total"] = sum(marker_counts.get(m, 0) for m in PHAGE_MODELS)

    vp_completeness, vp_mcp, plv_count, vp_df = calculate_vp_metrics(
        marker_counts,
        marker_hits,
        protein_domain=protein_domain,
        plv_tree_count=plv_tree_count,
    )
    result["vp_completeness"] = vp_completeness
    result["vp_mcp"] = vp_mcp
    result["plv"] = plv_count
    result["vp_df"] = vp_df

    mirus_completeness, mirus_df = calculate_mirus_completeness(
        marker_counts, marker_hits, protein_domain=protein_domain
    )
    result["mirus_completeness"] = mirus_completeness
    result["mirus_df"] = mirus_df

    # Unified capsid-type tally across the MCP panels (NCLDV/Mirus/PPV caps
    # groups), e.g. "Nucleocytoviricota:4,Gossevirus:1". Needs the MCP tree-NN
    # placements; blank on a counts-only resume that skipped the MCP trees.
    result["capsid_group"] = capsid_group
    return result


def _set_default_marker_metrics(result: Dict[str, Any]) -> None:
    """Apply default marker metric values when counts are unavailable."""
    # Standalone count + virophage detail columns, plus the internal
    # contamination-model features (cellular_*/phage_unique/phage_total) that
    # are read by the trained model but not emitted as summary columns.
    for metric in [
        "ncldv_mcp_total",
        "vp_mcp",
        "plv",
        "cellular_unique",
        "cellular_total",
        "phage_unique",
        "phage_total",
    ]:
        result[metric] = 0

    for metric in [
        "gvog4_dup",
        "gvog8_dup",
        "busco_dup",
        "cog_dup",
        "mrya_dup",
        "phage_dup",
        "vp_df",
        "mirus_df",
    ]:
        result[metric] = 0.0

    result["gvog4_completeness"] = f"0/{len(GVOG4M_MODELS)}"
    result["gvog8_completeness"] = f"0/{len(GVOG8M_MODELS)}"
    result["busco_completeness"] = f"0/{len(BUSCO_MODELS)}"
    result["cog_completeness"] = f"0/{len(UNI56_MODELS)}"
    result["mrya_completeness"] = f"0/{len(MRYA_MODELS)}"
    result["phage_completeness"] = f"0/{len(PHAGE_MODELS)}"
    result["vp_completeness"] = f"0/{len(VP_CATEGORY_PREFIXES)}"
    result["mirus_completeness"] = f"0/{len(MIRUS_CATEGORY_MODELS)}"
    result["capsid_group"] = ""
