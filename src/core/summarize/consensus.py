"""Taxonomy consensus for the GVClass summary.

Turns the per-marker neighbor lineages into the ``taxonomy_majority`` string and
its ``taxonomy_confidence`` flags.

At each rank, every marker votes for its unique top taxon among its placements
under the lineage already selected. The unique top taxon across markers wins.
At domain, exactly one viral domain tied with cellular domains wins while the
tie flag stays. Other ties, too few supporting markers, or a eukaryotic call
without enough cellular support stop the lineage at that rank.
"""

from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import TypeAlias

from src.core.summarize.tax_format import (
    _aggregate_confidence_flags,
    _format_tax_label,
    _lineage_matches_prefix,
)

Lineage: TypeAlias = tuple[str, ...]

TAX_LEVELS = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

# Each level has a minimum number of distinct supporting markers required
# before a consensus assignment is accepted. Below that floor the level
# is emitted as an empty label (``d_``, ``p_``, ...) and
# ``taxonomy_confidence`` is set to ``low_support``. Order-level markers
# are skipped under ``mode_fast=True``, so the order threshold drops to
# 2 in that regime and ``taxonomy_confidence`` becomes
# ``reduced_fastmode`` when the relaxed threshold is what let the
# assignment through.
_MIN_MARKERS_FOR_LEVEL: dict[str, int] = {
    "domain": 2,
    "phylum": 2,
    "class": 2,
    "order": 3,
    "family": 2,
    "genus": 2,
    "species": 2,
}
_MIN_MARKERS_FOR_LEVEL_FAST: dict[str, int] = {
    **_MIN_MARKERS_FOR_LEVEL,
    "order": 2,
}

# A eukaryotic domain call must hold at least this share of the markers that
# vote for a cellular domain. Without it, a plurality of two markers decides:
# on 1,470 prokaryotic gut genomes distant homologs and contaminated
# eukaryotic references gave a eukaryotic call to 59 of them (median share
# 0.58), while 12 of 12 eukaryotic single cells held at least 0.95.
# Bacterial, archaeal and viral calls keep the plurality rule.
MIN_CELLULAR_DOMAIN_SHARE = 0.8

# Domain label -> the family it counts toward in the share test. Putative
# endogenous viral element references count with their eukaryotic source.
_CELLULAR_DOMAIN_FAMILY: dict[str, str] = {
    "BAC": "BAC",
    "ARC": "ARC",
    "EUK": "EUK",
    "EUK-pEVE": "EUK",
}
# Families whose call needs the share. Their false calls are the ones seen.
_SHARE_REQUIRED_FAMILIES = frozenset({"EUK"})


@dataclass(frozen=True, slots=True)
class _RankVote:
    """Across-marker vote at one rank under the selected parent lineage.

    ``total`` counts the markers with a placement at the rank, including markers
    that cast no vote because their own top two taxa tie. ``winner`` is None
    when no marker voted or when tied leaders cannot be resolved. A viral
    domain can win a tie with cellular domains while ``tied`` stays True.
    """

    winner: str | None
    supporting: int
    total: int
    tied: bool
    votes: Counter[str]


def _build_consensus_taxonomies(
    per_marker_lineage_counters: Mapping[str, Counter[Lineage]],
    *,
    mode_fast: bool,
) -> tuple[str, str]:
    """Return ``(taxonomy_majority, taxonomy_confidence)`` from marker lineages.

    Ranks are decided top-down. Each rank considers only placements under the
    lineage selected at the ranks above. The first rank that is not accepted
    leaves itself and all lower ranks blank. A flag that more than one rank sets
    names the shallowest of them. The result does not depend on the order of the
    markers or of their lineages.
    """
    labels: list[str] = []
    flag_ranks: dict[str, str] = {}
    selected_prefix: Lineage = ()
    for rank_index, level in enumerate(TAX_LEVELS):
        vote = _rank_vote(per_marker_lineage_counters, rank_index, selected_prefix)
        winner, rank_flags = _judge_rank(vote, level, mode_fast=mode_fast)
        for flag in rank_flags:
            flag_ranks.setdefault(flag, level)
        if winner is None:
            break
        labels.append(_format_tax_label(level, winner))
        selected_prefix += (winner,)
    labels.extend(f"{level[0]}_" for level in TAX_LEVELS[len(labels) :])
    return ";".join(labels), _aggregate_confidence_flags(flag_ranks)


def _min_markers_for_level(level: str, *, mode_fast: bool) -> int:
    threshold_map = _MIN_MARKERS_FOR_LEVEL_FAST if mode_fast else _MIN_MARKERS_FOR_LEVEL
    return threshold_map[level]


def _judge_rank(
    vote: _RankVote, level: str, *, mode_fast: bool
) -> tuple[str | None, list[str]]:
    """Return the accepted taxon (None when the rank stops) and the rank's flags."""
    if vote.tied and vote.winner is None:
        return None, ["tie"]
    if vote.total == 0:
        # No placement at this rank under the selected parent: the references
        # carry no label here. At the domain rank this means no evidence at all.
        return None, ["no_support"] if level == "domain" else []

    flags = ["tie"] if vote.tied else []
    threshold = _min_markers_for_level(level, mode_fast=mode_fast)
    if vote.winner is None or vote.supporting < threshold:
        return None, ["low_support", *flags]
    if level == "domain" and not _cellular_domain_supported(vote.votes, vote.winner):
        return None, ["mixed_cellular"]

    if 2 * vote.supporting <= vote.total:
        flags.append("plurality")
    # The call passed the threshold above. Support below the extended-mode
    # floor means the relaxed fast-mode threshold let it through. A fast-mode
    # run with enough order markers is not downgraded.
    if vote.supporting < _min_markers_for_level(level, mode_fast=False):
        flags.append("reduced_fastmode")
    return vote.winner, flags


def _rank_vote(
    per_marker_lineage_counters: Mapping[str, Counter[Lineage]],
    rank_index: int,
    selected_prefix: Lineage,
) -> _RankVote:
    """Collect one vote per marker at a rank and pick the across-marker winner."""
    votes: Counter[str] = Counter()
    total = 0
    for lineage_counter in per_marker_lineage_counters.values():
        rank_counter = _marker_rank_counter(
            lineage_counter, rank_index, selected_prefix
        )
        if not rank_counter:
            continue
        total += 1
        marker_top, _ = _unique_top(rank_counter)
        if marker_top is not None:
            votes[marker_top] += 1
    winner, tied = _unique_top(votes)
    if tied and rank_index == 0:
        winner = _viral_domain_tie_winner(votes)
    supporting = votes[winner] if winner is not None else 0
    return _RankVote(winner, supporting, total, tied, votes)


def _marker_rank_counter(
    lineage_counter: Counter[Lineage],
    rank_index: int,
    selected_prefix: Lineage,
) -> Counter[str]:
    """Count one marker's placements per taxon at a rank under the parent."""
    rank_counter: Counter[str] = Counter()
    for lineage, count in lineage_counter.items():
        if not _lineage_matches_prefix(lineage, selected_prefix):
            continue
        if rank_index >= len(lineage):
            continue
        taxon = lineage[rank_index].strip()
        if taxon:
            rank_counter[taxon] += count
    return rank_counter


def _unique_top(counter: Counter[str]) -> tuple[str | None, bool]:
    """Return ``(taxon, tied)`` for the entry with the highest count.

    ``taxon`` is None when the counter is empty or when the two highest counts
    are equal (``tied`` is then True). Only counts decide, so the result does
    not depend on insertion order.
    """
    top_two = counter.most_common(2)
    if not top_two:
        return None, False
    if len(top_two) == 2 and top_two[0][1] == top_two[1][1]:
        return None, True
    return top_two[0][0], False


def _viral_domain_tie_winner(votes: Counter[str]) -> str | None:
    """Return the sole viral domain tied for the lead with cellular domains."""
    from src.core.contamination_scoring import (  # noqa: PLC0415  # labels imports consensus rank types.
        CELLULAR_PREFIXES,
        GIANT_VIRUS_PREFIXES,
        PHAGE_PREFIXES,
    )

    top_count = max(votes.values())
    tied_domains = {domain for domain, count in votes.items() if count == top_count}
    viral_domains = tied_domains & (GIANT_VIRUS_PREFIXES | PHAGE_PREFIXES)
    cellular_domains = CELLULAR_PREFIXES | _CELLULAR_DOMAIN_FAMILY.keys()
    if len(viral_domains) == 1 and tied_domains - viral_domains <= cellular_domains:
        return viral_domains.pop()
    return None


def _cellular_domain_supported(marker_votes: Counter[str], winner: str) -> bool:
    """Return whether a winning eukaryotic domain holds enough cellular votes.

    The share is taken among markers that vote for a cellular domain, so viral
    markers do not dilute it. Any other winner is always supported.
    """
    family = _CELLULAR_DOMAIN_FAMILY.get(winner)
    if family not in _SHARE_REQUIRED_FAMILIES:
        return True
    cellular_votes = 0
    family_votes = 0
    for domain, count in marker_votes.items():
        domain_family = _CELLULAR_DOMAIN_FAMILY.get(domain)
        if domain_family is None:
            continue
        cellular_votes += count
        if domain_family == family:
            family_votes += count
    return family_votes >= MIN_CELLULAR_DOMAIN_SHARE * cellular_votes
