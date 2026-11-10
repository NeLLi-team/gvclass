"""A cellular domain needs 80% of the cellular-domain markers, not a plurality."""

from __future__ import annotations

from collections import Counter

from src.core.summarize import consensus

LINEAGES = {
    "EUK": (
        "EUK",
        "Ascomycota",
        "Sordariomycetes",
        "Hypocreales",
        "Nectriaceae",
        "g1",
        "s1",
    ),
    "BAC": (
        "BAC",
        "Pseudomonadota",
        "Gammaproteobacteria",
        "Enterobacterales",
        "Enterobacteriaceae",
        "g2",
        "s2",
    ),
    "ARC": (
        "ARC",
        "Thermoproteota",
        "Nitrososphaeria",
        "Nitrososphaerales",
        "Nitrososphaeraceae",
        "g3",
        "s3",
    ),
    "NCLDV": (
        "NCLDV",
        "Nucleocytoviricota",
        "Megaviricetes",
        "Imitervirales",
        "IM_01",
        "g4",
        "s4",
    ),
}


def _votes(**markers: int) -> dict[str, Counter]:
    """Return one marker per requested vote, each placed on its domain lineage."""
    return {
        f"{domain}_{index}": Counter({LINEAGES[domain]: 1})
        for domain, count in markers.items()
        for index in range(count)
    }


def _call(**markers: int) -> tuple[str, str]:
    return consensus._build_consensus_taxonomies(_votes(**markers), mode_fast=False)


def _domain_label(majority: str) -> str:
    return majority.partition(";")[0]


def test_cellular_plurality_without_supermajority_leaves_domain_unassigned() -> None:
    majority, confidence = _call(EUK=6, BAC=5)

    assert majority == "d_;p_;c_;o_;f_;g_;s_"
    assert confidence == "mixed_cellular"


def test_cellular_supermajority_keeps_domain_and_lower_ranks() -> None:
    majority, confidence = _call(EUK=9, BAC=1)

    assert _domain_label(majority) == "d_EUK"
    assert "p_Ascomycota" in majority
    assert confidence == "high"


def test_threshold_is_inclusive() -> None:
    majority, _ = _call(EUK=4, BAC=1)

    assert _domain_label(majority) == "d_EUK"


def test_viral_winner_is_not_subject_to_the_share_rule() -> None:
    majority, confidence = _call(NCLDV=6, EUK=4)

    assert _domain_label(majority) == "d_NCLDV"
    assert "mixed_cellular" not in confidence


def test_viral_markers_do_not_dilute_the_cellular_share() -> None:
    # EUK wins the plurality over NCLDV; among cellular markers it holds 10 of 11.
    majority, _ = _call(EUK=10, NCLDV=8, BAC=1)

    assert _domain_label(majority) == "d_EUK"


def test_bacterial_and_archaeal_pluralities_keep_the_call() -> None:
    majority, confidence = _call(BAC=7, ARC=3)

    assert _domain_label(majority) == "d_BAC"
    assert "mixed_cellular" not in confidence


def test_eukaryotic_votes_are_diluted_by_bacterial_and_archaeal_markers() -> None:
    majority, _ = _call(EUK=8, BAC=1, ARC=1)

    assert _domain_label(majority) == "d_EUK"


def test_equal_eukaryotic_and_bacterial_votes_are_a_tie_not_mixed_cellular() -> None:
    majority, confidence = _call(EUK=5, BAC=5)

    assert majority == "d_;p_;c_;o_;f_;g_;s_"
    assert confidence == "tie"
