"""Tests for species-tree Section 7: per-domain panel config (NCLDV/PPV/MIRUS).

The orchestration is domain-agnostic (parameterized by SpeciesTreePanel), so PPV
and MIRUS are pure config additions. These tests pin panel registration, domain
routing, and marker-group -> reference-faa resolution.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.core.species_tree.config import (
    MIRUS_PANEL,
    NCLDV_PANEL,
    NEIGHBORS_PER_COMBINED_TREE,
    NEIGHBORS_PER_QUERY_TREE,
    PANELS,
    PEVE_PREFIX,
    PPV_PANEL,
    groups_for_models,
    neighbors_to_store,
    select_panel,
)
from src.utils.resource_store import ResourceStore

REPO = Path(__file__).resolve().parents[1]
RESOURCE_DIR = REPO / "resources"
STORE = ResourceStore(RESOURCE_DIR)


def test_all_three_panels_registered() -> None:
    assert set(PANELS) == {"NCLDV", "PPV", "MIRUS"}


def test_neighbor_knobs_and_store_count() -> None:
    # Adjustable defaults: per-query tree (the default product) is broader than
    # the opt-in combined tree.
    assert NEIGHBORS_PER_QUERY_TREE == 30
    assert NEIGHBORS_PER_COMBINED_TREE == 20
    # Sidecars must retain enough neighbors to satisfy whichever k is larger so
    # both trees can be served from the same stored data.
    assert neighbors_to_store() == max(
        NEIGHBORS_PER_QUERY_TREE, NEIGHBORS_PER_COMBINED_TREE
    )


def test_select_panel_routes_by_domain_token() -> None:
    assert select_panel("d_NCLDV;p_Nucleocytoviricota").name == "NCLDV"
    assert select_panel("d_PPV;p_Preplasmiviricota").name == "PPV"
    assert select_panel("d_MIRUS;p_Mirusviricota").name == "MIRUS"
    # Non-registered domains and empty/low-support are a no-op.
    assert select_panel("d_BAC;p_Pseudomonadota") is None
    assert select_panel("d_EUK;p_Chordata") is None
    assert select_panel("d_EUK-pEVE;p_Discosea-pEVE") is None
    assert select_panel("") is None
    assert select_panel("low_support") is None


def test_peve_is_auxiliary_reference_prefix_not_routing_panel() -> None:
    assert set(PANELS) == {"NCLDV", "PPV", "MIRUS"}
    for panel in PANELS.values():
        assert panel.reference_prefixes == (panel.domain_prefix, PEVE_PREFIX)


def test_gvog8_resolves_to_eight_distinct_groups() -> None:
    assert NCLDV_PANEL.groups == (
        "grp_GVOGm0013",
        "grp_COG0085",
        "grp_COG0086",
        "polb",
        "tfiib_cyclin",
        "GVOGm0461",
        "atpase_ncldv_mrya",
        "vltf3",
    )


def test_mirus_resolves_to_six_groups() -> None:
    assert len(MIRUS_PANEL.groups) == 6
    assert "mcp_mirus" in MIRUS_PANEL.groups
    assert "terminase" in MIRUS_PANEL.groups  # ATPase + merged collapse to one
    assert MIRUS_PANEL.min_markers == 3


def test_ppv_core_groups() -> None:
    assert PPV_PANEL.groups == ("mcp_vp", "mcp_plv", "penton_vp", "atpase_vp")
    assert PPV_PANEL.min_markers == 2


def test_ppv_and_mirus_thresholds_are_half_panel() -> None:
    assert PPV_PANEL.min_markers * 2 == len(PPV_PANEL.groups)
    assert MIRUS_PANEL.min_markers * 2 == len(MIRUS_PANEL.groups)


def test_groups_for_models_dedupes_preserving_order() -> None:
    # Two models mapping to the same group collapse once; order preserved.
    assert groups_for_models(["GVOGm0022", "GVOGm0023", "GVOGm0022"]) == (
        "grp_COG0085",
        "grp_COG0086",
    )


def _count_prefixed_headers(faa: Path, prefix: str) -> int:
    token = ">" + prefix
    count = 0
    with open(faa) as handle:
        for line in handle:
            if line.startswith(token):
                count += 1
    return count


@pytest.mark.requires_db
@pytest.mark.skipif(
    not ((RESOURCE_DIR / "database" / "faa").exists() or STORE.has_parquet_faa()),
    reason="reference DB not installed",
)
@pytest.mark.parametrize("panel", [NCLDV_PANEL, PPV_PANEL, MIRUS_PANEL])
def test_panel_groups_resolve_to_faa_with_domain_refs(panel) -> None:
    for group in panel.groups:
        faa = STORE.marker_faa_path(group)
        assert faa.exists(), f"{panel.name} group {group} has no reference faa"
        n_refs = _count_prefixed_headers(faa, panel.domain_prefix)
        assert (
            n_refs >= panel.min_markers * 50
        ), f"{panel.name}/{group}: only {n_refs} {panel.domain_prefix} refs"
