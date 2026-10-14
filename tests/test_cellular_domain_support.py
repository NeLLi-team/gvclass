"""A cellular domain needs 80% of the cellular-domain markers, not a plurality."""

from __future__ import annotations

from collections import Counter, defaultdict

import pytest

from src.core.summarize import consensus, tax_format

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


def _votes(**markers: int):
    """Return flat, per-marker and lineage counters with one lineage vote per marker."""
    flat: dict[str, Counter] = {level: Counter() for level in consensus.TAX_LEVELS}
    per_marker: dict[str, dict[str, Counter]] = {
        level: defaultdict(Counter) for level in consensus.TAX_LEVELS
    }
    lineage_counters: dict[str, Counter[tuple[str, ...]]] = defaultdict(Counter)
    for domain, count in markers.items():
        lineage = LINEAGES[domain]
        for index in range(count):
            marker = f"{domain}_{index}"
            for level, value in zip(consensus.TAX_LEVELS, lineage, strict=True):
                key = tax_format._build_taxonomy_key(level, domain, value)
                flat[level][key] += 1
                per_marker[level][marker][key] += 1
            lineage_counters[marker][lineage] += 1
    return flat, per_marker, lineage_counters


def _domain_label(majority: str) -> str:
    return majority.partition(";")[0]


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_cellular_plurality_without_supermajority_leaves_domain_unassigned(
    path,
) -> None:
    flat, per_marker, lineages = _votes(EUK=6, BAC=5)

    majority, strict, confidence = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert majority == "d_;p_;c_;o_;f_;g_;s_"
    assert strict == "d_;p_;c_;o_;f_;g_;s_"
    assert confidence == "mixed_cellular"


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_cellular_supermajority_keeps_domain_and_lower_ranks(path) -> None:
    flat, per_marker, lineages = _votes(EUK=9, BAC=1)

    majority, _, confidence = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_EUK"
    assert "p_Ascomycota" in majority
    assert confidence == "high"


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_threshold_is_inclusive(path) -> None:
    flat, per_marker, lineages = _votes(EUK=4, BAC=1)

    majority, _, _ = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_EUK"


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_viral_winner_is_not_subject_to_the_share_rule(path) -> None:
    flat, per_marker, lineages = _votes(NCLDV=6, EUK=4)

    majority, _, confidence = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_NCLDV"
    assert "mixed_cellular" not in confidence


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_viral_markers_do_not_dilute_the_cellular_share(path) -> None:
    # EUK wins the plurality over NCLDV; among cellular markers it holds 10 of 11.
    flat, per_marker, lineages = _votes(EUK=10, NCLDV=8, BAC=1)

    majority, _, _ = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_EUK"


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_bacterial_and_archaeal_pluralities_keep_the_call(path) -> None:
    flat, per_marker, lineages = _votes(BAC=7, ARC=3)

    majority, _, confidence = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_BAC"
    assert "mixed_cellular" not in confidence


@pytest.mark.parametrize("path", ["lineage", "flat"])
def test_eukaryotic_votes_are_diluted_by_bacterial_and_archaeal_markers(path) -> None:
    flat, per_marker, lineages = _votes(EUK=8, BAC=1, ARC=1)

    majority, _, _ = consensus._build_consensus_taxonomies(
        flat, per_marker, False, lineages if path == "lineage" else None
    )

    assert _domain_label(majority) == "d_EUK"


def test_mixed_cellular_is_reported_before_reduced_fastmode() -> None:
    assert (
        tax_format._aggregate_confidence_flags(["reduced_fastmode", "mixed_cellular"])
        == "mixed_cellular,reduced_fastmode"
    )
