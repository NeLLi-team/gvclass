"""Taxonomy consensus: per-marker votes, ties, support floors and confidence flags."""

from __future__ import annotations

from collections import Counter

import pytest

from src.core.summarize import consensus, order_completeness, tax_format

Lineage = tuple[str, ...]


def _ncldv(order: str, family: str = "", genus: str = "", species: str = "") -> Lineage:
    return (
        "NCLDV",
        "Nucleocytoviricota",
        "Megaviricetes",
        order,
        family,
        genus,
        species,
    )


IMITERVIRALES = _ncldv("Imitervirales", "Mimiviridae", "g1", "s1")
ALGAVIRALES = _ncldv("Algavirales", "Phycodnaviridae", "g2", "s2")
BACTERIA = (
    "BAC",
    "Pseudomonadota",
    "Gammaproteobacteria",
    "Enterobacterales",
    "Enterobacteriaceae",
    "g3",
    "s3",
)

NCLDV_PREFIX = "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes"
IMITERVIRALES_CALL = f"{NCLDV_PREFIX};o_Imitervirales;f_Mimiviridae;g_g1;s_s1"
ALGAVIRALES_CALL = f"{NCLDV_PREFIX};o_Algavirales;f_Phycodnaviridae;g_g2;s_s2"
ORDER_UNASSIGNED_CALL = f"{NCLDV_PREFIX};o_;f_;g_;s_"
UNASSIGNED_CALL = "d_;p_;c_;o_;f_;g_;s_"


def _markers(prefix: str, count: int, lineage: Lineage) -> dict[str, Counter]:
    """Return ``count`` markers that each place once on ``lineage``."""
    return {f"{prefix}{index}": Counter({lineage: 1}) for index in range(count)}


def _consensus(
    markers: dict[str, Counter], *, mode_fast: bool = False
) -> tuple[str, str]:
    return consensus._build_consensus_taxonomies(markers, mode_fast=mode_fast)


def test_equal_order_votes_leave_the_order_blank_and_flag_a_tie() -> None:
    markers = {
        **_markers("imi", 3, IMITERVIRALES),
        **_markers("alg", 3, ALGAVIRALES),
    }

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "tie:order"


def test_viral_cellular_domain_tie_keeps_viral_lineage_and_flags_tie() -> None:
    markers = {
        **_markers("ncldv", 3, _ncldv("Asfuvirales", "Asfarviridae", "g4", "s4")),
        **_markers("euk", 3, ("EUK", "EukPhylum", "EukClass", "EukOrder")),
    }

    majority, confidence = _consensus(markers)

    assert majority == (
        "d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;"
        "o_Asfuvirales;f_Asfarviridae;g_g4;s_s4"
    )
    assert confidence == "tie:domain,plurality:domain"


@pytest.mark.parametrize("viral_domain", ["MIRUS", "PHAGE", "VP", "PLV", "PPV", "MRYA"])
def test_other_viral_domains_win_cellular_domain_ties(viral_domain: str) -> None:
    markers = {
        **_markers("viral", 3, (viral_domain,)),
        **_markers("bac", 3, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == f"d_{viral_domain};p_;c_;o_;f_;g_;s_"
    assert confidence == "tie:domain,plurality:domain"


@pytest.mark.parametrize("cellular_domain", ["ARC", "EUK-pEVE", "MITO", "PLASTID"])
def test_viral_domain_wins_other_cellular_domain_ties(cellular_domain: str) -> None:
    markers = {
        **_markers("ncldv", 3, IMITERVIRALES),
        **_markers("cellular", 3, (cellular_domain,)),
    }

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "tie:domain,plurality:domain"


def test_one_viral_domain_wins_a_three_way_tie_with_two_cellular_domains() -> None:
    markers = {
        **_markers("ncldv", 3, IMITERVIRALES),
        **_markers("bac", 3, BACTERIA),
        **_markers("euk", 3, ("EUK",)),
    }

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "tie:domain,plurality:domain"


@pytest.mark.parametrize("cellular_markers", [0, 3])
def test_two_viral_domains_tied_leave_every_rank_blank(cellular_markers: int) -> None:
    markers = {
        **_markers("ncldv", 3, IMITERVIRALES),
        **_markers("mirus", 3, ("MIRUS",)),
        **_markers("bac", cellular_markers, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == UNASSIGNED_CALL
    assert confidence == "tie:domain"


def test_two_cellular_domains_tied_leave_every_rank_blank() -> None:
    markers = {
        **_markers("bac", 3, BACTERIA),
        **_markers("euk", 3, ("EUK",)),
        **_markers("ncldv", 2, IMITERVIRALES),
    }

    majority, confidence = _consensus(markers)

    assert majority == UNASSIGNED_CALL
    assert confidence == "tie:domain"


def test_viral_domain_tie_does_not_resolve_a_lower_order_tie() -> None:
    markers = {
        **_markers("imi", 3, IMITERVIRALES),
        **_markers("alg", 3, ALGAVIRALES),
        **_markers("bac", 6, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "tie:domain,plurality:domain"


def test_viral_domain_tie_still_requires_normal_order_support() -> None:
    markers = {
        **_markers("ncldv", 2, IMITERVIRALES),
        **_markers("bac", 2, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "low_support,tie:domain,plurality:domain"


def test_viral_domain_tie_still_requires_normal_domain_support() -> None:
    markers = {
        **_markers("ncldv", 1, IMITERVIRALES),
        **_markers("bac", 1, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == UNASSIGNED_CALL
    assert confidence == "low_support,tie:domain"


def test_viral_cellular_domain_vote_without_a_tie_is_unchanged() -> None:
    markers = {
        **_markers("ncldv", 3, IMITERVIRALES),
        **_markers("bac", 2, BACTERIA),
    }

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "high"


def test_three_way_phylum_tie_is_a_tie_not_low_support() -> None:
    markers = {
        "m1": Counter({("NCLDV", "Phylum_A", "Class_Z"): 1}),
        "m2": Counter({("NCLDV", "Phylum_B", "Class_Z"): 1}),
        "m3": Counter({("NCLDV", "Phylum_C", "Class_Z"): 1}),
    }

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == "d_NCLDV;p_;c_;o_;f_;g_;s_"
    assert confidence == "tie:phylum"


def test_tie_at_species_keeps_the_ranks_above() -> None:
    markers = {
        **_markers("a", 2, _ncldv("Imitervirales", "Mimiviridae", "g1", "s_a")),
        **_markers("b", 2, _ncldv("Imitervirales", "Mimiviridae", "g1", "s_b")),
    }

    majority, confidence = _consensus(markers)

    assert majority == f"{NCLDV_PREFIX};o_Imitervirales;f_Mimiviridae;g_g1;s_"
    assert confidence == "tie:species"


def test_three_of_twelve_markers_are_called_but_flagged_plurality() -> None:
    others = {
        f"other{index}": Counter({_ncldv(f"Order_{index}", "F", "g", "s"): 1})
        for index in range(9)
    }
    markers = {**_markers("imi", 3, IMITERVIRALES), **others}

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "plurality:order"


def test_exactly_half_of_the_markers_is_a_plurality() -> None:
    others = {
        f"other{index}": Counter({_ncldv(f"Order_{index}", "F", "g", "s"): 1})
        for index in range(3)
    }
    markers = {**_markers("imi", 3, IMITERVIRALES), **others}

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "plurality:order"


def test_majority_of_markers_stays_high() -> None:
    markers = {
        **_markers("imi", 5, IMITERVIRALES),
        **_markers("alg", 1, ALGAVIRALES),
    }

    majority, confidence = _consensus(markers)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "high"


def test_marker_with_internal_tie_casts_no_vote_but_counts_in_the_total() -> None:
    both_orders = Counter({ALGAVIRALES: 1, _ncldv("Pithovirales", "Pithoviridae"): 1})
    markers = {
        **_markers("imi", 3, IMITERVIRALES),
        **{f"tied{index}": Counter(both_orders) for index in range(4)},
    }

    majority, confidence = _consensus(markers)

    # Algavirales would win with four alphabetical votes if tied markers voted.
    assert majority == IMITERVIRALES_CALL
    assert confidence == "plurality:order"


def test_marker_votes_for_the_taxon_with_most_placements() -> None:
    markers = {
        **_markers("imi", 2, IMITERVIRALES),
        "paralogs": Counter({IMITERVIRALES: 1, ALGAVIRALES: 3}),
        **_markers("alg", 2, ALGAVIRALES),
    }

    majority, confidence = _consensus(markers)

    assert majority == ALGAVIRALES_CALL
    assert confidence == "high"


def test_markers_that_all_tie_internally_give_low_support() -> None:
    both_orders = Counter({IMITERVIRALES: 1, ALGAVIRALES: 1})
    markers = {f"tied{index}": Counter(both_orders) for index in range(3)}

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "low_support"


def test_one_paralog_rich_marker_does_not_outvote_distinct_markers() -> None:
    markers = {
        **_markers("alg", 3, ALGAVIRALES),
        "paralog": Counter({IMITERVIRALES: 10}),
    }

    majority, confidence = _consensus(markers)

    assert majority == ALGAVIRALES_CALL
    assert confidence == "high"


@pytest.mark.parametrize("reverse", [False, True])
def test_result_does_not_depend_on_marker_or_lineage_order(reverse: bool) -> None:
    ordered = [
        ("m1", [IMITERVIRALES, ALGAVIRALES]),
        ("m2", [IMITERVIRALES]),
        ("m3", [ALGAVIRALES]),
        ("m4", [IMITERVIRALES]),
        ("m5", [ALGAVIRALES]),
        ("m6", [ALGAVIRALES, IMITERVIRALES]),
    ]
    if reverse:
        ordered = [(name, lineages[::-1]) for name, lineages in ordered[::-1]]
    markers = {name: Counter(dict.fromkeys(lineages, 1)) for name, lineages in ordered}

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "tie:order"


def test_below_the_order_floor_gives_low_support() -> None:
    markers = _markers("imi", 2, IMITERVIRALES)

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "low_support"


def test_no_markers_give_no_support() -> None:
    majority, confidence = _consensus({})

    assert majority == UNASSIGNED_CALL
    assert confidence == "no_support"


def test_ranks_without_labels_under_the_parent_are_blank_without_a_flag() -> None:
    markers = _markers("m", 3, _ncldv(""))

    majority, confidence = _consensus(markers)

    assert majority == ORDER_UNASSIGNED_CALL
    assert confidence == "high"


def test_fast_mode_relaxes_the_order_floor_and_flags_it() -> None:
    markers = _markers("imi", 2, IMITERVIRALES)

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "reduced_fastmode"


def test_fast_mode_is_not_flagged_when_support_meets_the_standard_floor() -> None:
    markers = _markers("imi", 3, IMITERVIRALES)

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "high"


def test_plurality_and_reduced_fastmode_are_listed_in_fixed_order() -> None:
    others = {
        f"other{index}": Counter({_ncldv(f"Order_{index}", "F", "g", "s"): 1})
        for index in range(3)
    }
    markers = {**_markers("imi", 2, IMITERVIRALES), **others}

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == IMITERVIRALES_CALL
    assert confidence == "plurality:order,reduced_fastmode"


def test_confidence_flags_are_joined_in_fixed_order() -> None:
    flag_ranks = {
        "reduced_fastmode": "order",
        "plurality": "genus",
        "mixed_cellular": "domain",
        "tie": "family",
        "low_support": "phylum",
        "no_support": "domain",
    }

    assert tax_format._aggregate_confidence_flags(flag_ranks) == (
        "no_support,low_support,tie:family,mixed_cellular,plurality:genus,"
        "reduced_fastmode"
    )


def test_plurality_names_the_shallowest_rank_that_has_one() -> None:
    # Family: 4 of 10 markers vote F1. Genus under F1: 2 of 4 markers vote g1.
    # Both ranks are a plurality; the flag names the family.
    def lineage(family: str, genus: str) -> Lineage:
        return _ncldv("Imitervirales", family, genus, "s1")

    markers = {
        "m1": Counter({lineage("F1", "g1"): 1}),
        "m2": Counter({lineage("F1", "g1"): 1}),
        "m3": Counter({lineage("F1", "g1"): 1, lineage("F1", "g2"): 1}),
        "m4": Counter({lineage("F1", "g3"): 1, lineage("F1", "g4"): 1}),
        **{
            f"other{index}": Counter({lineage(f"F_other{index}", "g"): 1})
            for index in range(6)
        },
    }

    majority, confidence = _consensus(markers)

    assert majority == f"{NCLDV_PREFIX};o_Imitervirales;f_F1;g_g1;s_s1"
    assert confidence == "plurality:family"


def test_plurality_at_the_genus_names_the_genus() -> None:
    def lineage(genus: str, species: str) -> Lineage:
        return _ncldv("Imitervirales", "Mimiviridae", genus, species)

    markers = {
        **_markers("g1_", 2, lineage("g1", "s1")),
        "g2": Counter({lineage("g2", "s2"): 1}),
        "g3": Counter({lineage("g3", "s3"): 1}),
        "g4": Counter({lineage("g4", "s4"): 1}),
    }

    majority, confidence = _consensus(markers)

    assert majority == f"{NCLDV_PREFIX};o_Imitervirales;f_Mimiviridae;g_g1;s_s1"
    assert confidence == "plurality:genus"


def test_confidence_without_flags_is_high() -> None:
    assert tax_format._aggregate_confidence_flags([]) == "high"


def test_lower_ranks_stay_in_the_majority_domain() -> None:
    mirus_class_04 = (
        "MIRUS",
        "Mirusviricota",
        "Class_04",
        "Demutovirales",
        "Demuto_02",
        "MIRUS_unclassified",
    )
    mirus_unclassified = (
        "MIRUS",
        "Mirusviricota",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "MIRUS_unclassified",
        "TARA_PSW_NCLDV_00076",
    )
    eukaryote = (
        "EUK",
        "EUK_unclassified",
        "EUK_unclassified",
        "EUK_unclassified",
        "EUK_unclassified",
        "EUK_unclassified",
        "Stramenopiles sp. TOSAG23-3",
    )
    markers = {
        **_markers("g838_", 3, (*mirus_class_04, "MIRUS__Delmont2025__MIRUS_G_0838")),
        **_markers("g1322_", 1, (*mirus_class_04, "MIRUS__Delmont2025__MIRUS_G_1322")),
        **_markers("tara_", 2, mirus_unclassified),
        **_markers("euk_", 3, eukaryote),
    }

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == (
        "d_MIRUS;p_Mirusviricota;c_Class_04;o_Demutovirales;"
        "f_Demuto_02;g_MIRUS_unclassified;s_MIRUS_G_0838"
    )
    assert confidence == "high"


def test_child_ranks_stay_in_the_selected_parent() -> None:
    def lineage(phylum: str, class_: str) -> Lineage:
        return ("NCLDV", phylum, class_, "Order", "Family", "Genus", "Species")

    markers = {
        **_markers("a1_", 3, lineage("Phylum_A", "Class_A1")),
        **_markers("a2_", 2, lineage("Phylum_A", "Class_A2")),
        **_markers("b_", 4, lineage("Phylum_B", "Class_B")),
    }

    majority, confidence = _consensus(markers, mode_fast=True)

    assert majority == (
        "d_NCLDV;p_Phylum_A;c_Class_A1;o_Order;f_Family;g_Genus;s_Species"
    )
    assert confidence == "high"


def test_quality_group_follows_the_accepted_taxonomy_not_paralog_counts() -> None:
    asfuvirales = _ncldv("Asfuvirales", "Asfarviridae")
    markers = {
        **_markers("asfu", 3, asfuvirales),
        "paralog": Counter({IMITERVIRALES: 10}),
    }

    majority, _ = _consensus(markers)
    accepted_lineage = tax_format.parse_taxonomy_majority(majority)

    assert accepted_lineage["order"] == "Asfuvirales"
    assert order_completeness._select_quality_group(accepted_lineage) == "Asfuvirales"


def test_tax_level_counts_normalize_domain_unclassified_placeholders() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV__NCLDV_unclassified": 6, "NCLDV__Megaviricetes": 4}),
        level="class",
    )

    assert "NCLDV__unclassified:6(60.00%)" in formatted
    assert "NCLDV__NCLDV_unclassified" not in formatted


def test_tax_level_counts_use_domain_prefix_for_grouped_minor_taxa() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV__S392": 18, "PLV__EP00001": 1, "NCLDV__singleton": 1}),
        level="species",
    )

    assert "PLV__other:1(5.00%)" in formatted
    assert "NCLDV__other:1(5.00%)" in formatted
    assert "PLV_other" not in formatted
    assert "NCLDV_other" not in formatted


def test_domain_counts_are_not_grouped_as_other() -> None:
    formatted = tax_format.format_tax_level_counts(
        Counter({"NCLDV": 19, "PLV": 1}),
        level="domain",
    )

    assert "NCLDV:19(95.00%)" in formatted
    assert "PLV:1(5.00%)" in formatted
    assert "PLV_other" not in formatted
