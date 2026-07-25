"""
Taxonomy consensus for the GVClass summary.

Stateless counterparts of the ``FullSummarizer`` methods that turn per-marker
taxonomy counters and neighbor lineages into consensus taxonomy strings; the
class keeps thin delegating wrappers for every name here.
"""

from collections import Counter
from typing import Dict, List, Optional, Tuple

from src.core.summarize.tax_format import (
    _aggregate_confidence_flags,
    _build_taxonomy_key,
    _format_tax_label,
    _lineage_matches_prefix,
    _stable_top,
    _tax_key_matches_domain,
    get_tax_consensus,
)

TAX_LEVELS = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

# --- Per-marker majority -------------------------------------------------
#
# Each level has a minimum number of distinct supporting markers required
# before a consensus assignment is accepted. Below that floor the level
# is emitted as an empty label (``d_``, ``p_``, ...) and
# ``taxonomy_confidence`` is set to ``low_support``. Order-level markers
# are skipped under ``mode_fast=True``, so the order threshold drops to
# 2 in that regime and ``taxonomy_confidence`` becomes
# ``reduced_fastmode`` when the relaxed threshold is what let the
# assignment through.
_MIN_MARKERS_FOR_LEVEL: Dict[str, int] = {
    "domain": 2,
    "phylum": 2,
    "class": 2,
    "order": 3,
    "family": 2,
    "genus": 2,
    "species": 2,
}
_MIN_MARKERS_FOR_LEVEL_FAST: Dict[str, int] = {
    **_MIN_MARKERS_FOR_LEVEL,
    "order": 2,
}


def _min_markers_for_level(level: str, mode_fast: bool) -> int:
    threshold_map = _MIN_MARKERS_FOR_LEVEL_FAST if mode_fast else _MIN_MARKERS_FOR_LEVEL
    return threshold_map.get(level, 2)


def _per_marker_majority(
    per_marker_level_counters: Dict[str, Counter],
) -> Tuple[Optional[str], int, int]:
    """Tally one winning taxon per marker and return the across-marker vote.

    For each marker, the per-marker winner is chosen with a stable
    lexicographic tiebreaker on taxon name. The across-marker winner is
    chosen the same way, so two pipelines processing the same data
    cannot disagree on the consensus call purely because they iterated
    markers in different filesystem orders.

    Returns ``(winning_taxon, supporting_markers, total_markers)``.
    """
    marker_votes: Counter = Counter()
    total_markers = 0
    for counter in per_marker_level_counters.values():
        if not counter:
            continue
        total_markers += 1
        winner_pair = _stable_top(counter)
        if winner_pair is None:
            continue
        marker_votes[winner_pair[0]] += 1

    if not marker_votes:
        return None, 0, total_markers

    winning_taxon, supporting_markers = _stable_top(marker_votes)
    return winning_taxon, supporting_markers, total_markers


def _filter_counter_to_domain(counter: Counter, domain: Optional[str]) -> Counter:
    if not domain:
        return counter
    return Counter(
        {
            tax_key: count
            for tax_key, count in counter.items()
            if _tax_key_matches_domain(tax_key, domain)
        }
    )


def _filter_per_marker_counters_to_domain(
    per_marker_level_counters: Dict[str, Counter],
    domain: Optional[str],
) -> Dict[str, Counter]:
    if not domain:
        return per_marker_level_counters
    filtered: Dict[str, Counter] = {}
    for marker, counter in per_marker_level_counters.items():
        marker_counter = _filter_counter_to_domain(counter, domain)
        if marker_counter:
            filtered[marker] = marker_counter
    return filtered


def _build_consensus_taxonomies(
    tax_counters: Dict[str, Counter],
    per_marker_counters: Dict[str, Dict[str, Counter]],
    mode_fast: bool,
    per_marker_lineage_counters: Optional[Dict[str, Counter[Tuple[str, ...]]]] = None,
) -> Tuple[str, str, str]:
    """Build majority and strict taxonomy strings and a confidence label.

    The majority call now uses a per-marker vote: one winning taxon per
    marker, then the across-marker winner. A minimum number of distinct
    supporting markers is required (see ``_min_markers_for_level``);
    otherwise the level is emitted unassigned. Strict consensus
    continues to require 100% agreement within the flat counter for
    backward compatibility with pre-1.4.3 callers.
    """
    if per_marker_lineage_counters is not None:
        return _build_consensus_taxonomies_from_lineages(
            per_marker_lineage_counters, mode_fast
        )

    majority_parts = []
    strict_parts = []
    confidence_flags: List[str] = []
    saw_taxonomy_evidence = False
    selected_domain: Optional[str] = None

    for level in TAX_LEVELS:
        flat_counter = tax_counters[level]
        per_marker_level = per_marker_counters.get(level, {})
        if level != "domain":
            flat_counter = _filter_counter_to_domain(flat_counter, selected_domain)
            per_marker_level = _filter_per_marker_counters_to_domain(
                per_marker_level, selected_domain
            )
        threshold = _min_markers_for_level(level, mode_fast)
        winning, supporting, total = _per_marker_majority(per_marker_level)
        if flat_counter or total > 0:
            saw_taxonomy_evidence = True

        if winning is not None and supporting >= threshold:
            majority = _format_tax_label(level, winning)
            if level == "domain":
                selected_domain = winning
            # Only flag reduced_fastmode when the relaxed threshold is
            # what actually allowed the call — i.e. the support falls
            # below the standard-mode threshold but clears the fast-mode
            # one. A fast-mode run that already has 3+ supporting
            # order-level markers should not be downgraded.
            standard_threshold = _MIN_MARKERS_FOR_LEVEL.get(level, 2)
            if (
                mode_fast
                and threshold < standard_threshold
                and supporting < standard_threshold
            ):
                confidence_flags.append("reduced_fastmode")
        else:
            majority = f"{level[0]}_"
            if total > 0:
                confidence_flags.append("low_support")

        # Strict keeps the flat-counter 100% rule, but must never be more
        # specific than majority: if the per-marker majority withheld this
        # level (insufficient distinct-marker support), strict is withheld
        # too, so the "conservative" column is always at least as
        # conservative as the majority column.
        _, strict = get_tax_consensus(flat_counter, level)
        if majority == f"{level[0]}_":
            strict = f"{level[0]}_"
        majority_parts.append(majority)
        strict_parts.append(strict)

    if not saw_taxonomy_evidence:
        confidence_flags.append("no_support")

    majority_str = ";".join(majority_parts)
    strict_str = ";".join(strict_parts)
    confidence = _aggregate_confidence_flags(confidence_flags)
    return majority_str, strict_str, confidence


def _lineage_rank_counter(
    per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
    rank_index: int,
    selected_prefix: Tuple[str, ...],
) -> Counter:
    rank_counter: Counter = Counter()
    for marker_counter in per_marker_lineage_counters.values():
        for lineage, count in marker_counter.items():
            if not _lineage_matches_prefix(lineage, selected_prefix):
                continue
            if rank_index >= len(lineage):
                continue
            tax_value = lineage[rank_index].strip()
            if tax_value:
                rank_counter[tax_value] += count
    return rank_counter


def _per_marker_lineage_majority(
    per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
    rank_index: int,
    selected_prefix: Tuple[str, ...],
) -> Tuple[Optional[str], int, int]:
    marker_votes: Counter = Counter()
    total_markers = 0

    for marker_counter in per_marker_lineage_counters.values():
        rank_counter: Counter = Counter()
        for lineage, count in marker_counter.items():
            if not _lineage_matches_prefix(lineage, selected_prefix):
                continue
            if rank_index >= len(lineage):
                continue
            tax_value = lineage[rank_index].strip()
            if tax_value:
                rank_counter[tax_value] += count
        if not rank_counter:
            continue
        total_markers += 1
        winner_pair = _stable_top(rank_counter)
        if winner_pair is None:
            continue
        marker_votes[winner_pair[0]] += 1

    if not marker_votes:
        return None, 0, total_markers

    winning_taxon, supporting_markers = _stable_top(marker_votes)
    return winning_taxon, supporting_markers, total_markers


def _format_lineage_tax_label(
    level: str, selected_prefix: Tuple[str, ...], tax_value: str
) -> str:
    if not tax_value:
        return f"{level[0]}_"
    if level == "domain":
        return _format_tax_label(level, tax_value)
    domain = selected_prefix[0] if selected_prefix else ""
    return _format_tax_label(level, _build_taxonomy_key(level, domain, tax_value))


def _build_consensus_taxonomies_from_lineages(
    per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
    mode_fast: bool,
) -> Tuple[str, str, str]:
    """Build top-down consensus taxonomy from full neighbor lineages."""
    majority_parts = []
    strict_parts = []
    confidence_flags: List[str] = []
    saw_taxonomy_evidence = False
    selected_prefix: Tuple[str, ...] = ()
    strict_prefix: Optional[Tuple[str, ...]] = ()
    stop_lower_ranks = False

    for rank_index, level in enumerate(TAX_LEVELS):
        blank = f"{level[0]}_"
        if stop_lower_ranks:
            majority_parts.append(blank)
            strict_parts.append(blank)
            continue

        flat_counter = _lineage_rank_counter(
            per_marker_lineage_counters, rank_index, selected_prefix
        )
        threshold = _min_markers_for_level(level, mode_fast)
        winning, supporting, total = _per_marker_lineage_majority(
            per_marker_lineage_counters, rank_index, selected_prefix
        )
        if flat_counter or total > 0:
            saw_taxonomy_evidence = True

        if winning is not None and supporting >= threshold:
            majority = _format_lineage_tax_label(level, selected_prefix, winning)
            strict = blank
            if strict_prefix is not None:
                strict_counter = _lineage_rank_counter(
                    per_marker_lineage_counters,
                    rank_index,
                    strict_prefix,
                )
                if len(strict_counter) == 1:
                    strict_value = next(iter(strict_counter))
                    if strict_value == winning:
                        strict = _format_lineage_tax_label(
                            level, strict_prefix, strict_value
                        )
                        strict_prefix = strict_prefix + (strict_value,)
                    else:
                        strict_prefix = None
                else:
                    strict_prefix = None
            selected_prefix = selected_prefix + (winning,)

            standard_threshold = _MIN_MARKERS_FOR_LEVEL.get(level, 2)
            if (
                mode_fast
                and threshold < standard_threshold
                and supporting < standard_threshold
            ):
                confidence_flags.append("reduced_fastmode")
        else:
            majority = blank
            strict = blank
            stop_lower_ranks = True
            strict_prefix = None
            if total > 0:
                confidence_flags.append("low_support")

        majority_parts.append(majority)
        strict_parts.append(strict)

    if not saw_taxonomy_evidence:
        confidence_flags.append("no_support")

    majority_str = ";".join(majority_parts)
    strict_str = ";".join(strict_parts)
    confidence = _aggregate_confidence_flags(confidence_flags)
    return majority_str, strict_str, confidence
