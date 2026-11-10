"""capscan caps marker groups and their derivation from hits and tree placements."""

import pytest

from src.config.marker_sets import CAPS_MARKER_GROUP
from src.core.labels import LabelIndex


def test_caps_marker_group_map():
    assert CAPS_MARKER_GROUP["plv_mcp_caps_PgVV_Aquinto"] == "PgVV"
    assert CAPS_MARKER_GROUP["plv_mcp_caps_GKS2_Alpenseevirus_NCV_like"] == "Alpenseevirus"
    assert all(isinstance(v, str) and v for v in CAPS_MARKER_GROUP.values())


# --- tree-NN capscan_group (the phase-C primary path) ------------------------
from src.core.summarize import labels_taxonomy

_TAXONOMY = {
    "PPV__pgvv": [
        "PPV",
        "Preplasmiviricota",
        "Aquintoviricetes",
        "Aquintoviricetes_order",
        "PgVV",
        "PgVV_genus",
        "PgVV_sp",
    ],
    "PPV__sp": [
        "PPV",
        "Preplasmiviricota",
        "Aquintoviricetes",
        "Aquintoviricetes_order",
        "SP",
        "SP_genus",
        "SP_sp",
    ],
    "PPV__generic": [
        "PPV",
        "Preplasmiviricota",
        "Polintoviricetes",
        "PLV_unclassified",
        "PLV_unclassified",
        "PLV_unclassified",
        "acc",
    ],
    "NCLDV__zeph": [
        "NCLDV",
        "Nucleocytoviricota",
        "Zephyrvirus_NCV_like",
        "Zephyrvirus_order",
        "Zephyrvirus",
        "Zephyrvirus_genus",
        "Zephyrvirus_sp",
    ],
    "NCLDV__realncldv": [
        "NCLDV",
        "Nucleocytoviricota",
        "Megaviricetes",
        "Imitervirales",
        "IM_01",
        "g1",
        "sp1",
    ],
}
_TAXONOMY["MIRUS__mir"] = ["MIRUS", "Mirusviricota", "c", "o", "f", "g", "s"]
_LBL = LabelIndex.from_taxonomy(_TAXONOMY)


# --- unified capsid_group counter --------------------------------------------
def test_capsid_group_phylum_fallback_ncldv():
    # generic NCLDV (family not a caps group) -> phylum label
    assert (
        labels_taxonomy.capsid_group_counter(
            _LBL, {"mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.1}}}
        )
        == "Nucleocytoviricota:1"
    )


def test_capsid_group_mirus_phylum():
    assert (
        labels_taxonomy.capsid_group_counter(
            _LBL, {"mcp_mirus": {"q1": {"MIRUS__mir|p": 0.1}}}
        )
        == "Mirusviricota:1"
    )


def test_capsid_group_caps_family_preferred_over_phylum():
    # NCLDV__zeph carries caps family Zephyrvirus AND phylum Nucleocytoviricota;
    # the caps family wins.
    assert (
        labels_taxonomy.capsid_group_counter(
            _LBL, {"mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.1}}}
        )
        == "Zephyrvirus:1"
    )


def test_capsid_group_generic_ppv_skipped():
    # PLV_unclassified family + Preplasmiviricota phylum (not a capsid phylum) -> no label
    assert (
        labels_taxonomy.capsid_group_counter(
            _LBL, {"mcp_plv": {"q1": {"PPV__generic|p": 0.1}}}
        )
        == ""
    )


def test_capsid_group_counts_sorted_by_count_then_label():
    tnn = {
        "mcp_ncldv": {"a": {"NCLDV__realncldv|p": 0.1}, "b": {"NCLDV__realncldv|p": 0.1}},
        "mcp_plv": {"c": {"PPV__pgvv|p": 0.1}},
    }
    assert (
        labels_taxonomy.capsid_group_counter(_LBL, tnn) == "Nucleocytoviricota:2,PgVV:1"
    )


def test_capsid_group_cross_panel_closest_wins():
    # same protein in two MCP trees -> one vote for the closest placement
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.2}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.5}}}
    assert labels_taxonomy.capsid_group_counter(_LBL, tnn) == "PgVV:1"


def test_capsid_group_equal_distance_cross_type_dropped():
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.3}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.3}}}
    assert labels_taxonomy.capsid_group_counter(_LBL, tnn) == ""


def test_capsid_group_empty():
    assert labels_taxonomy.capsid_group_counter(_LBL, {}) == ""


# --- tree-aware plv count (shared A32 must not flag NCLDV) --------------------
def test_plv_count_tree_aware_ppv_placement_counted():
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"PPV__pgvv|p": 0.1}}}
    assert labels_taxonomy.plv_count_tree_aware(_LBL, mh, tnn) == 1


def test_plv_count_tree_aware_ncldv_placement_not_counted():
    # the A32 case: hits the PLV HMM but places with an NCLDV reference -> not PLV
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"NCLDV__realncldv|p": 0.1}}}
    assert labels_taxonomy.plv_count_tree_aware(_LBL, mh, tnn) == 0


def test_plv_count_tree_aware_fallback_to_hmm_without_tree():
    mh = {"PLV_PC_054": {"q1", "q2"}}
    assert labels_taxonomy.plv_count_tree_aware(_LBL, mh, None) == 2


# --- tied nearest neighbours: lineage in common, any insertion order ----------
_TIED_LBL = LabelIndex.from_taxonomy(
    {**_TAXONOMY, "PPV__pgvv_twin": list(_TAXONOMY["PPV__pgvv"])}
)


def _tied(neighbors: list[tuple[str, float]], reverse: bool) -> dict[str, float]:
    return dict(neighbors[::-1] if reverse else neighbors)


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_between_caps_groups_falls_back_to_their_phylum(reverse):
    # PgVV and SP share Preplasmiviricota, which is not a capsid phylum.
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy.capsid_group_counter(_TIED_LBL, tnn) == ""


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_between_caps_and_generic_ncldv_uses_the_phylum(reverse):
    neighbors = _tied([("NCLDV__zeph|p", 0.1), ("NCLDV__realncldv|p", 0.1)], reverse)
    tnn = {"mcp_ncldv": {"q1": neighbors}}
    assert (
        labels_taxonomy.capsid_group_counter(_TIED_LBL, tnn) == "Nucleocytoviricota:1"
    )


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_tie_inside_one_caps_group_keeps_the_group(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__pgvv_twin|p", 0.1)], reverse)
    tnn = {"mcp_plv": {"q1": neighbors}}
    assert labels_taxonomy.capsid_group_counter(_TIED_LBL, tnn) == "PgVV:1"


@pytest.mark.parametrize("reverse", [False, True])
def test_capsid_group_unreportable_nearest_placement_casts_no_vote(reverse):
    # PgVV and SP tie at 0.1 in the PPV tree and share no capsid label. The NCLDV
    # placement at 0.2 is farther, so it must not take over the vote.
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {
        "mcp_plv": {"q1": neighbors},
        "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.2}},
    }
    assert labels_taxonomy.capsid_group_counter(_TIED_LBL, tnn) == ""


def test_capsid_group_tie_across_trees_inside_one_caps_group_keeps_the_group():
    tnn = {
        "mcp_plv": {"q1": {"PPV__pgvv|p": 0.1}},
        "mcp_ncldv": {"q1": {"PPV__pgvv_twin|p": 0.1}},
    }
    assert labels_taxonomy.capsid_group_counter(_TIED_LBL, tnn) == "PgVV:1"


@pytest.mark.parametrize("reverse", [False, True])
def test_plv_count_tree_aware_tie_between_ppv_and_ncldv_is_not_counted(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("NCLDV__realncldv|p", 0.1)], reverse)
    tnn = {"PLV_PC_054": {"q1": neighbors}}
    mh = {"PLV_PC_054": {"q1"}}
    assert labels_taxonomy.plv_count_tree_aware(_TIED_LBL, mh, tnn) == 0


@pytest.mark.parametrize("reverse", [False, True])
def test_plv_count_tree_aware_tie_between_ppv_references_is_counted(reverse):
    neighbors = _tied([("PPV__pgvv|p", 0.1), ("PPV__sp|p", 0.1)], reverse)
    tnn = {"PLV_PC_054": {"q1": neighbors}}
    mh = {"PLV_PC_054": {"q1"}}
    assert labels_taxonomy.plv_count_tree_aware(_TIED_LBL, mh, tnn) == 1
