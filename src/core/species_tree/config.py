"""Per-domain species-tree panel registry (Section 6 seam; Section 7 extends it).

A :class:`SpeciesTreePanel` parameterizes the whole route by
``(domain_prefix, marker_panel, min_markers, gate_token)`` so adding PPV / MIRUS
in Section 7 is a config addition, not new pipeline code. The default route
registers only NCLDV (GVOG8).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

from src.config.marker_sets import (
    GVOG8M_MODELS,
    MIRUS_CATEGORY_MODELS,
    MODEL_TO_GROUP,
)

# --- Neighbor breadth (adjustable knobs) ------------------------------------
# Top-k nearest reference *genomes* kept per marker gene tree. The per-query
# species tree (the default product) uses the larger set for resolution; the
# opt-in combined one-leaf-per-dataset tree uses the smaller set for a lighter
# supermatrix. Both are read at run time, so editing them here changes behaviour
# with no other code change.
NEIGHBORS_PER_QUERY_TREE = 30
NEIGHBORS_PER_COMBINED_TREE = 20


def neighbors_to_store() -> int:
    """Neighbors a sidecar must retain so either tree's k is satisfiable.

    The per-query hook selects (and stores) this many neighbors per group; the
    per-query and combined assembly steps each truncate to their own k, so the
    two counts stay independent even though the combined step reuses the
    per-query sidecars.
    """
    return max(NEIGHBORS_PER_QUERY_TREE, NEIGHBORS_PER_COMBINED_TREE)


@dataclass(frozen=True)
class SpeciesTreePanel:
    """One domain's species-tree configuration.

    Attributes:
        name: Human-readable panel name (e.g. ``"NCLDV"``).
        domain_prefix: Reference-leaf domain gate (e.g. ``"NCLDV__"``) — only
            references with this prefix are eligible neighbors.
        gate_token: ``taxonomy_majority`` prefix that routes a query to this panel
            (e.g. ``"d_NCLDV"``).
        groups: Ordered marker-group names (the supermatrix partitions).
        min_markers: Minimum groups a taxon must have to enter the supermatrix.
    """

    name: str
    domain_prefix: str
    gate_token: str
    groups: Tuple[str, ...]
    min_markers: int = 3


def groups_for_models(models: Sequence[str]) -> Tuple[str, ...]:
    """Map marker models to their group/faa basenames, order-preserving + unique.

    Uses ``MODEL_TO_GROUP`` (combinable groups) and falls back to the model name
    for ungrouped singletons, matching how ``extract_marker_hits`` and the
    reference faa files are named.
    """
    ordered = []
    for model in models:
        group = MODEL_TO_GROUP.get(model, model)
        if group not in ordered:
            ordered.append(group)
    return tuple(ordered)


NCLDV_PANEL = SpeciesTreePanel(
    name="NCLDV",
    domain_prefix="NCLDV__",
    gate_token="d_NCLDV",
    groups=groups_for_models(GVOG8M_MODELS),
    min_markers=3,
)

# PPV (Preplasmiviricota) has no GVOG-style panel. Its core structural markers are
# the conserved virophage/PLV capsid + packaging groups carrying the most PPV
# references (mcp_vp ~7.1k, mcp_plv ~5.6k, atpase_vp ~4.4k, penton_vp ~4.0k PPV__
# headers). PPV genomes are small, so the marker floor is 2/4.
PPV_PANEL = SpeciesTreePanel(
    name="PPV",
    domain_prefix="PPV__",
    gate_token="d_PPV",
    groups=("mcp_vp", "mcp_plv", "penton_vp", "atpase_vp"),
    min_markers=2,
)

# MIRUS (Mirusviricota) core = the four marker categories, grouped (MCP, terminase
# ATPase, portal, triplex), resolving to six group faa.
_MIRUS_MODELS = [model for models in MIRUS_CATEGORY_MODELS.values() for model in models]
MIRUS_PANEL = SpeciesTreePanel(
    name="MIRUS",
    domain_prefix="MIRUS__",
    gate_token="d_MIRUS",
    groups=groups_for_models(_MIRUS_MODELS),
    min_markers=3,
)

# Registry of all domain panels routed by select_panel().
PANELS: Dict[str, SpeciesTreePanel] = {
    panel.name: panel for panel in (NCLDV_PANEL, PPV_PANEL, MIRUS_PANEL)
}


def select_panel(taxonomy_majority: str) -> Optional[SpeciesTreePanel]:
    """Return the panel whose ``gate_token`` prefixes ``taxonomy_majority``.

    A genome is routed to at most one domain panel based on the domain gvclass
    already called for it; returns ``None`` (a logged no-op) for any genome not
    matching a registered panel.
    """
    if not taxonomy_majority:
        return None
    for panel in PANELS.values():
        if taxonomy_majority.startswith(panel.gate_token):
            return panel
    return None
