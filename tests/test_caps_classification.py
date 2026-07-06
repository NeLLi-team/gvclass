"""capscan caps markers now follow the standard per-model GA cutoff logic (no special floor):
GA enforced in normal mode, dropped in sensitive mode, exactly like every other model."""
from src.core.hmm_search import _should_keep_filtered_hit
from src.config.marker_sets import CAPS_MARKER_GROUP


def test_caps_markers_follow_standard_ga():
    m = "plv_mcp_caps_PgVV_Aquinto"
    # sensitive mode (no GA cutoffs): caps kept regardless of score, like any other model
    assert _should_keep_filtered_hit(m, 50.0, 50.0, False, {}) is True
    # normal mode with the stamped GA=75: enforced the standard way (both seq and dom)
    assert _should_keep_filtered_hit(m, 80.0, 80.0, True, {m: (75.0, 75.0)}) is True
    assert _should_keep_filtered_hit(m, 74.9, 80.0, True, {m: (75.0, 75.0)}) is False
    assert _should_keep_filtered_hit(m, 80.0, 74.9, True, {m: (75.0, 75.0)}) is False


def test_non_caps_markers_follow_standard_ga():
    # no GA cutoffs -> kept regardless of score (sensitive-mode behaviour)
    assert _should_keep_filtered_hit("VP_MCP_A", 10.0, 10.0, False, {}) is True
    # GA cutoffs applied in normal mode
    assert _should_keep_filtered_hit("OG1352", 50.0, 50.0, True, {"OG1352": (60.0, 60.0)}) is False
    assert _should_keep_filtered_hit("OG1352", 65.0, 65.0, True, {"OG1352": (60.0, 60.0)}) is True


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
from src.core.summarize_full import FullSummarizer

_LBL = {
    "PPV__pgvv": ["PPV", "Preplasmiviricota", "Aquintoviricetes", "Aquintoviricetes_order", "PgVV", "PgVV_genus", "PgVV_sp"],
    "PPV__sp": ["PPV", "Preplasmiviricota", "Aquintoviricetes", "Aquintoviricetes_order", "SP", "SP_genus", "SP_sp"],
    "PPV__generic": ["PPV", "Preplasmiviricota", "Polintoviricetes", "PLV_unclassified", "PLV_unclassified", "PLV_unclassified", "acc"],
    "NCLDV__zeph": ["NCLDV", "Nucleocytoviricota", "Zephyrvirus_NCV_like", "Zephyrvirus_order", "Zephyrvirus", "Zephyrvirus_genus", "Zephyrvirus_sp"],
    "NCLDV__realncldv": ["NCLDV", "Nucleocytoviricota", "Megaviricetes", "Imitervirales", "IM_01", "g1", "sp1"],
}


def _summarizer():
    s = FullSummarizer.__new__(FullSummarizer)  # bypass heavy __init__
    s.labels_dict = _LBL
    return s


def test_caps_from_tree_single_caps_neighbor():
    assert _summarizer()._caps_group_from_tree({"mcp_plv": {"q1": {"PPV__pgvv|p": 0.1}}}) == "PgVV"


def test_caps_from_tree_generic_and_real_ncldv_neighbor_empty():
    s = _summarizer()
    assert s._caps_group_from_tree({"mcp_plv": {"q1": {"PPV__generic|p": 0.1}}}) == ""  # PLV_unclassified
    assert s._caps_group_from_tree({"mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.1}}}) == ""  # IM_01


def test_caps_from_tree_tie_and_dominant():
    s = _summarizer()
    assert s._caps_group_from_tree({"mcp_plv": {"a": {"PPV__pgvv|p": 0.1}, "b": {"PPV__sp|p": 0.1}}}) == ""
    assert s._caps_group_from_tree(
        {"mcp_plv": {"a": {"PPV__pgvv|p": 0.1}, "b": {"PPV__pgvv|p": 0.1}, "c": {"PPV__sp|p": 0.1}}}
    ) == "PgVV"


def test_caps_from_tree_nearest_within_marker_wins():
    assert _summarizer()._caps_group_from_tree(
        {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.05, "PPV__sp|p2": 0.9}}}
    ) == "PgVV"


def test_caps_from_tree_no_double_count_across_markers():
    # SAME protein q1 in BOTH trees; must vote ONCE, resolving to the closest placement
    # (PgVV at 0.2, not Zephyrvirus at 0.5). Pre-fix this double-counted -> tie -> "".
    s = _summarizer()
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.2}}, "mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.5}}}
    assert s._caps_group_from_tree(tnn) == "PgVV"


def test_caps_from_tree_equal_distance_cross_marker_is_ambiguous():
    # SAME protein at EQUAL distance to PgVV (mcp_plv) and Zephyrvirus (mcp_ncldv):
    # ambiguous -> dropped -> "" (not silently resolved by marker iteration order).
    s = _summarizer()
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.3}}, "mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.3}}}
    assert s._caps_group_from_tree(tnn) == ""


def test_caps_from_tree_empty():
    assert _summarizer()._caps_group_from_tree({}) == ""


# --- unified capsid_group counter --------------------------------------------
_LBL["MIRUS__mir"] = ["MIRUS", "Mirusviricota", "c", "o", "f", "g", "s"]


def test_capsid_group_phylum_fallback_ncldv():
    # generic NCLDV (family not a caps group) -> phylum label
    assert _summarizer()._capsid_group_counter(
        {"mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.1}}}
    ) == "Nucleocytoviricota:1"


def test_capsid_group_mirus_phylum():
    assert _summarizer()._capsid_group_counter(
        {"mcp_mirus": {"q1": {"MIRUS__mir|p": 0.1}}}
    ) == "Mirusviricota:1"


def test_capsid_group_caps_family_preferred_over_phylum():
    # NCLDV__zeph carries caps family Zephyrvirus AND phylum Nucleocytoviricota;
    # the caps family wins.
    assert _summarizer()._capsid_group_counter(
        {"mcp_ncldv": {"q1": {"NCLDV__zeph|p": 0.1}}}
    ) == "Zephyrvirus:1"


def test_capsid_group_generic_ppv_skipped():
    # PLV_unclassified family + Preplasmiviricota phylum (not a capsid phylum) -> no label
    assert _summarizer()._capsid_group_counter(
        {"mcp_plv": {"q1": {"PPV__generic|p": 0.1}}}
    ) == ""


def test_capsid_group_counts_sorted_by_count_then_label():
    tnn = {
        "mcp_ncldv": {"a": {"NCLDV__realncldv|p": 0.1}, "b": {"NCLDV__realncldv|p": 0.1}},
        "mcp_plv": {"c": {"PPV__pgvv|p": 0.1}},
    }
    assert _summarizer()._capsid_group_counter(tnn) == "Nucleocytoviricota:2,PgVV:1"


def test_capsid_group_cross_panel_closest_wins():
    # same protein in two MCP trees -> one vote for the closest placement
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.2}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.5}}}
    assert _summarizer()._capsid_group_counter(tnn) == "PgVV:1"


def test_capsid_group_equal_distance_cross_type_dropped():
    tnn = {"mcp_plv": {"q1": {"PPV__pgvv|p": 0.3}}, "mcp_ncldv": {"q1": {"NCLDV__realncldv|p": 0.3}}}
    assert _summarizer()._capsid_group_counter(tnn) == ""


def test_capsid_group_empty():
    assert _summarizer()._capsid_group_counter({}) == ""


# --- tree-aware plv count (shared A32 must not flag NCLDV) --------------------
def test_plv_count_tree_aware_ppv_placement_counted():
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"PPV__pgvv|p": 0.1}}}
    assert _summarizer()._plv_count_tree_aware(mh, tnn) == 1


def test_plv_count_tree_aware_ncldv_placement_not_counted():
    # the A32 case: hits the PLV HMM but places with an NCLDV reference -> not PLV
    mh = {"PLV_PC_054": {"q1"}}
    tnn = {"PLV_PC_054": {"q1": {"NCLDV__realncldv|p": 0.1}}}
    assert _summarizer()._plv_count_tree_aware(mh, tnn) == 0


def test_plv_count_tree_aware_fallback_to_hmm_without_tree():
    mh = {"PLV_PC_054": {"q1", "q2"}}
    assert _summarizer()._plv_count_tree_aware(mh, None) == 2
