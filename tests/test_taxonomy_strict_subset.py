"""M4: taxonomy_strict must never be more specific than taxonomy_majority.

When the per-marker majority withholds a level for insufficient distinct-marker
support, strict (legacy 100%-agreement on the flat counter) must also be blank
at that level, so the "conservative" column is never less conservative than the
majority column.
"""

from __future__ import annotations

from collections import Counter

from src.core.summarize_full import FullSummarizer


def _summarizer():
    # Bypass the heavy __init__ (model loading); the method under test only
    # uses class attributes and pure helper methods.
    return FullSummarizer.__new__(FullSummarizer)


def test_strict_not_more_specific_than_majority_below_floor():
    summ = _summarizer()
    levels = FullSummarizer.TAX_LEVELS

    tax_counters = {lvl: Counter() for lvl in levels}
    # One taxon unanimously present in the flat order-level counter ...
    tax_counters["order"] = Counter({"NCLDV__Imitervirales": 5})

    per_marker = {lvl: {} for lvl in levels}
    # ... but only ONE distinct marker supports it (order floor is 3), so the
    # majority rule withholds the order call.
    per_marker["order"] = {"markerA": Counter({"NCLDV__Imitervirales": 5})}

    majority_str, strict_str, _conf = summ._build_consensus_taxonomies(
        tax_counters, per_marker, mode_fast=False
    )
    majority = dict(zip(levels, majority_str.split(";")))
    strict = dict(zip(levels, strict_str.split(";")))

    assert majority["order"] == "o_"  # withheld for low support
    assert strict["order"] == "o_"  # strict must not be more specific


def test_strict_keeps_call_when_majority_assigns():
    summ = _summarizer()
    levels = FullSummarizer.TAX_LEVELS

    tax_counters = {lvl: Counter() for lvl in levels}
    tax_counters["order"] = Counter({"NCLDV__Imitervirales": 4})
    per_marker = {lvl: {} for lvl in levels}
    # Four distinct markers all agree -> clears the order floor (3).
    per_marker["order"] = {
        f"marker{i}": Counter({"NCLDV__Imitervirales": 1}) for i in range(4)
    }

    majority_str, strict_str, _conf = summ._build_consensus_taxonomies(
        tax_counters, per_marker, mode_fast=False
    )
    majority = dict(zip(levels, majority_str.split(";")))
    strict = dict(zip(levels, strict_str.split(";")))

    assert majority["order"] == "o_Imitervirales"
    assert strict["order"] == "o_Imitervirales"  # 100% agreement preserved
