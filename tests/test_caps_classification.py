"""capscan caps marker groups and their derivation from hits and tree placements."""

import pytest

from src.config.marker_sets import CAPS_MARKER_GROUP


def test_caps_marker_group_map():
    assert CAPS_MARKER_GROUP["plv_mcp_caps_PgVV_Aquinto"] == "PgVV"
    assert CAPS_MARKER_GROUP["plv_mcp_caps_GKS2_Alpenseevirus_NCV_like"] == "Alpenseevirus"
    assert all(isinstance(v, str) and v for v in CAPS_MARKER_GROUP.values())


from src.core.summarize_full import derive_capscan_group


def test_derive_capscan_group_best_per_protein():
    mh = {
        "plv_mcp_caps_PgVV_Aquinto": {"p1", "p2", "p3"},
        "plv_mcp_caps_Trimcap_A1_Aquinto": {"p4"},
        "plv_mcp_caps_Trimcap_A2_Aquinto": {"p4"},  # p4: 2 Trimcap markers -> one group vote
    }
    assert derive_capscan_group(mh, {}) == "PgVV"  # 3 PgVV proteins vs 1 Trimcap


def test_derive_capscan_group_profile_multiplicity_handled():
    # one protein cross-hits 3 Trimcap profiles + 1 PgVV -> per-protein dominant = Trimcap
    mh = {
        "plv_mcp_caps_Trimcap_A1_Aquinto": {"x"},
        "plv_mcp_caps_Trimcap_A2_Aquinto": {"x"},
        "plv_mcp_caps_Trimcap_A3_Aquinto": {"x"},
        "plv_mcp_caps_PgVV_Aquinto": {"x"},
    }
    assert derive_capscan_group(mh, {}) == "Trimcap_cluster_1"


def test_derive_capscan_group_tie_and_empty():
    assert derive_capscan_group({}, {}) == ""
    assert derive_capscan_group(None, {}) == ""
    mh = {"plv_mcp_caps_PgVV_Aquinto": {"p1"}, "plv_mcp_caps_SP_Aquinto": {"p2"}}
    assert derive_capscan_group(mh, {}) == ""  # 1 PgVV vs 1 SP -> tie


def test_derive_capscan_group_fallback_to_counts():
    assert derive_capscan_group(None, {"plv_mcp_caps_PgVV_Aquinto": 5}) == "PgVV"


# --- tree-NN capscan_group (the phase-C primary path) ------------------------
from src.core.summarize import labels_taxonomy

_LBL = {
    "PPV__pgvv": ["PPV", "Preplasmiviricota", "Aquintoviricetes", "Aquintoviricetes_order", "PgVV", "PgVV_genus", "PgVV_sp"],
    "PPV__sp": ["PPV", "Preplasmiviricota", "Aquintoviricetes", "Aquintoviricetes_order", "SP", "SP_genus", "SP_sp"],
    "PPV__generic": ["PPV", "Preplasmiviricota", "Polintoviricetes", "PLV_unclassified", "PLV_unclassified", "PLV_unclassified", "acc"],
    "NCLDV__zeph": ["NCLDV", "Nucleocytoviricota", "Zephyrvirus_NCV_like", "Zephyrvirus_order", "Zephyrvirus", "Zephyrvirus_genus", "Zephyrvirus_sp"],
    "NCLDV__realncldv": ["NCLDV", "Nucleocytoviricota", "Megaviricetes", "Imitervirales", "IM_01", "g1", "sp1"],
}


def test_caps_from_tree_single_caps_neighbor():
    assert labels_taxonomy._caps_group_from_tree(_LBL, {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.1}}}) == "PgVV"


def test_caps_from_tree_generic_and_real_ncldv_neighbor_empty():
    assert labels_taxonomy._caps_group_from_tree(_LBL, {"mcp_plv": {"q1": {"PPV__generic|p": 0.1}}}) == ""  # PLV_unclassified
    assert labels_taxonomy._caps_group_from_tree(_LBL, {"mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.1}}}) == ""  # IM_01


def test_caps_from_tree_tie_and_dominant():
    assert labels_taxonomy._caps_group_from_tree(_LBL, {"mcp_plv": {"a": {"PPV__pgvv|p": 0.1}, "b": {"PPV__sp|p": 0.1}}}) == ""
    assert labels_taxonomy._caps_group_from_tree(
        _LBL, {"mcp_plv": {"a": {"PPV__pgvv|p": 0.1}, "b": {"PPV__pgvv|p": 0.1}, "c": {"PPV__sp|p": 0.1}}}
    ) == "PgVV"


def test_caps_from_tree_nearest_within_marker_wins():
    assert labels_taxonomy._caps_group_from_tree(
        _LBL, {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.05, "PPV__sp|p2": 0.9}}}
    ) == "PgVV"


def test_caps_from_tree_no_double_count_across_markers():
    # SAME protein q1 in BOTH trees; must vote ONCE, resolving to the closest placement
    # (PgVV at 0.2, not Zephyrvirus at 0.5). Pre-fix this double-counted -> tie -> "".
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.2}}, "mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.5}}}
    assert labels_taxonomy._caps_group_from_tree(_LBL, tnn) == "PgVV"


def test_caps_from_tree_equal_distance_cross_marker_is_ambiguous():
    # SAME protein at EQUAL distance to PgVV (mcp_plv) and Zephyrvirus (mcp_ncldv):
    # ambiguous -> dropped -> "" (not silently resolved by marker iteration order).
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.3}}, "mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.3}}}
    assert labels_taxonomy._caps_group_from_tree(_LBL, tnn) == ""


def test_caps_from_tree_empty():
    assert labels_taxonomy._caps_group_from_tree(_LBL, {}) == ""


# --- unified capsid_group counter --------------------------------------------
_LBL["MIRUS__mir"] = ["MIRUS", "Mirusviricota", "c", "o", "f", "g", "s"]


def test_capsid_group_phylum_fallback_ncldv():
    # generic NCLDV (family not a caps group) -> phylum label
    assert labels_taxonomy._capsid_group_counter(
        _LBL, {"mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.1}}}
    ) == "Nucleocytoviricota:1"


def test_capsid_group_mirus_phylum():
    assert labels_taxonomy._capsid_group_counter(
        _LBL, {"mcp_mirus": {"q1": {"MIRUS__mir|p": 0.1}}}
    ) == "Mirusviricota:1"


def test_capsid_group_caps_family_preferred_over_phylum():
    # NCLDV__zeph carries caps family Zephyrvirus AND phylum Nucleocytoviricota;
    # the caps family wins.
    assert labels_taxonomy._capsid_group_counter(
        _LBL, {"mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.1}}}
    ) == "Zephyrvirus:1"


def test_capsid_group_generic_ppv_skipped():
    # PLV_unclassified family + Preplasmiviricota phylum (not a capsid phylum) -> no label
    assert labels_taxonomy._capsid_group_counter(
        _LBL, {"mcp_plv": {"q1": {"PPV__generic|p": 0.1}}}
    ) == ""


def test_capsid_group_counts_sorted_by_count_then_label():
    tnn = {
        "mcp_ncldv": {"a": {"NCLDV__realncldv|p": 0.1}, "b": {"NCLDV__realncldv|p": 0.1}},
        "mcp_plv": {"c": {"PPV__pgvv|p": 0.1}},
    }
    assert labels_taxonomy._capsid_group_counter(_LBL, tnn) == "Nucleocytoviricota:2,PgVV:1"


def test_capsid_group_cross_panel_closest_wins():
    # same protein in two MCP trees -> one vote for the closest placement
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.2}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.5}}}
    assert labels_taxonomy._capsid_group_counter(_LBL, tnn) == "PgVV:1"


def test_capsid_group_equal_distance_cross_type_dropped():
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.3}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.3}}}
    assert labels_taxonomy._capsid_group_counter(_LBL, tnn) == ""


def test_capsid_group_empty():
    assert labels_taxonomy._capsid_group_counter(_LBL, {}) == ""


# --- tree-aware plv count (shared A32 must not flag NCLDV) --------------------
def test_plv_count_tree_aware_ppv_placement_counted():
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"PPV__pgvv|p": 0.1}}}
    assert labels_taxonomy._plv_count_tree_aware(_LBL, mh, tnn) == 1


def test_plv_count_tree_aware_ncldv_placement_not_counted():
    # the A32 case: hits the PLV HMM but places with an NCLDV reference -> not PLV
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"NCLDV__realncldv|p": 0.1}}}
    assert labels_taxonomy._plv_count_tree_aware(_LBL, mh, tnn) == 0


def test_plv_count_tree_aware_fallback_to_hmm_without_tree():
    mh = {"PLV_PC_054": {"q1", "q2"}}
    assert labels_taxonomy._plv_count_tree_aware(_LBL, mh, None) == 2


# --- tied nearest neighbours: lineage in common, any insertion order ----------
_TIED_LBL = {**_LBL, "PPV__pgvv_twin": list(_LBL["PPV__pgvv"])}


def _tied(neighbors: list[tuple[str, float]], reverse: bool) -> dict[str, float]:
    return dict(neighbors[::-1] if reverse else neighbors)


@pytest.mark.parametrize("reverse", [False, True])
def test_caps_from_tree_tie_between_caps_groups_names_no_group(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy._caps_group_from_tree(_TIED_LBL, tnn) == ""


@pytest.mark.parametrize("reverse", [False, True])
def test_caps_from_tree_tie_between_caps_and_generic_names_no_group(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__generic|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy._caps_group_from_tree(_TIED_LBL, tnn) == ""


@pytest.mark.parametrize("reverse", [False, True])
def test_caps_from_tree_tie_inside_one_caps_group_keeps_the_group(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__pgvv_twin|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy._caps_group_from_tree(_TIED_LBL, tnn) == "PgVV"


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_between_caps_groups_falls_back_to_their_phylum(reverse):
    # PgVV and SP share Preplasmiviricota, which is not a capsid phylum.
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy._capsid_group_counter(_TIED_LBL, tnn) == ""


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_between_caps_and_generic_ncldv_uses_the_phylum(reverse):
    neighbors = _tied([("NCLDV__zeph|p", 0.1), ("NCLDV__realncldv|p", 0.1)], reverse)
    tnn = {"mcp_ncldv": {"q1": neighbors}}
    assert (
        labels_taxonomy._capsid_group_counter(_TIED_LBL, tnn) == "Nucleocytoviricota:1"
    )


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_inside_one_caps_group_keeps_the_group(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__pgvv_twin|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy._capsid_group_counter(_TIED_LBL, tnn) == "PgVV:1"


@pytest.mark.parametrize("reverse", [False, True])
def test_plv_count_tree_aware_tie_between_ppv_and_ncldv_is_not_counted(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("NCLDV__realncldv|p", 0.1)], reverse)
    tnn = {"PLV_PC_054": {"q1": neighbors}}
    mh = {"PLV_PC_054": {"q1"}}
    assert labels_taxonomy._plv_count_tree_aware(_TIED_LBL, mh, tnn) == 0


@pytest.mark.parametrize("reverse", [False, True])
def test_plv_count_tree_aware_tie_between_ppv_references_is_counted(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {"PLV_PC_054": {"q1": neighbors}}
    mh = {"PLV_PC_054": {"q1"}}
    assert labels_taxonomy._plv_count_tree_aware(_TIED_LBL, mh, tnn) == 1
