"""
Taxonomy label, key and count formatting for the GVClass summary.

Stateless helpers that turn taxonomy counters and lineage keys into output
strings.
"""

from collections import Counter
from typing import List, Optional, Tuple


def format_tax_level_counts(tax_counter: Counter, level: str | None = None) -> str:
    """Format taxonomy counts with percentages, grouping low-frequency taxa."""
    if not tax_counter:
        return ""

    normalized_counter = Counter()
    for tax, count in tax_counter.items():
        normalized_counter[_format_tax_count_key(tax, level)] += count

    total = sum(normalized_counter.values())
    if total <= 0:
        return ""

    grouped_counter = Counter(normalized_counter)
    low_threshold = 5.0
    low_counts_by_group = Counter()

    for tax, count in normalized_counter.items():
        percentage = (count / total) * 100
        group_key = _low_frequency_tax_group_key(tax, level)
        if percentage <= low_threshold and group_key:
            low_counts_by_group[group_key] += count
            grouped_counter.pop(tax, None)

    for group_key, count in low_counts_by_group.items():
        grouped_counter[group_key] += count

    sorted_counts = sorted(grouped_counter.items(), key=lambda x: x[1], reverse=True)

    formatted = []
    for tax, count in sorted_counts:
        percentage = (count / total) * 100
        formatted.append(f"{tax}:{count}({percentage:.2f}%)")

    return ",".join(formatted)


def _format_tax_count_key(tax_key: str, level: str | None = None) -> str:
    if level == "domain":
        return tax_key.split("__", 1)[0]

    if "__" not in tax_key:
        return tax_key

    domain, tax_name = tax_key.split("__", 1)
    if tax_name == f"{domain}_unclassified":
        tax_name = "unclassified"
    elif tax_name == f"{domain}_other":
        tax_name = "other"
    return f"{domain}__{tax_name}"


def _low_frequency_tax_group_key(tax_key: str, level: str | None = None) -> str:
    if level == "domain":
        return ""
    if "__" in tax_key:
        domain = tax_key.split("__", 1)[0]
        return f"{domain}__other"
    if "_" in tax_key:
        domain = tax_key.split("_", 1)[0]
        return f"{domain}__other"
    return f"{tax_key}_other"


def get_tax_consensus(tax_counter: Counter, level: str) -> Tuple[str, str]:
    """Get majority and strict consensus for a taxonomic level."""
    if not tax_counter:
        return f"{level[0]}_", f"{level[0]}_"

    total = sum(tax_counter.values())
    most_common = tax_counter.most_common(1)[0]
    tax, count = most_common

    # For domain level, use the full domain name (e.g., "NCLDV" -> "d_NCLDV")
    # For other levels, use the full taxonomy with prefix (e.g., "NCLDV__Mesomimiviridae" -> "f_Mesomimiviridae")
    if level == "domain":
        # Domain is just the prefix part (e.g., "NCLDV", "BAC", "EUK")
        tax_name = tax
    else:
        # For other levels, extract the name after the domain prefix
        if "__" in tax:
            parts = tax.split("__")
            tax_name = parts[-1] if len(parts) > 1 and parts[-1] else tax
        else:
            tax_name = tax

    # Strict consensus: 100% agreement
    strict = f"{level[0]}_{tax_name}" if count == total else f"{level[0]}_"

    # Majority consensus: >50% agreement
    majority = f"{level[0]}_{tax_name}" if count > total * 0.5 else f"{level[0]}_"

    return majority, strict


def _build_taxonomy_key(level: str, domain_prefix: str, tax_value: str) -> str:
    if level == "domain":
        return domain_prefix or tax_value
    if domain_prefix:
        return f"{domain_prefix}__{tax_value}"
    return tax_value


def _stable_top(counter: Counter) -> Optional[Tuple[str, int]]:
    """Return the counter entry with the highest count, breaking ties on
    lexicographic key so the result does not depend on upstream insertion
    order. ``Counter.most_common(1)`` alone is insertion-stable, which is
    non-deterministic when the feeder iterates an unordered glob
    (``tree_dir.glob("*.treefile")``)."""
    if not counter:
        return None
    # Sort descending by (count, key) with key as tiebreaker. Negative
    # count keeps the primary sort descending; key asc ensures stable
    # alphabetic tiebreak across runs.
    return min(counter.items(), key=lambda kv: (-kv[1], kv[0]))


def _tax_key_matches_domain(tax_key: str, domain: Optional[str]) -> bool:
    if not domain:
        return True
    return tax_key.startswith(f"{domain}__")


def _lineage_matches_prefix(
    lineage: Tuple[str, ...], selected_prefix: Tuple[str, ...]
) -> bool:
    return lineage[: len(selected_prefix)] == selected_prefix


def _format_tax_label(level: str, tax_key: str) -> str:
    if level == "domain":
        return f"{level[0]}_{tax_key}"
    if "__" in tax_key:
        parts = tax_key.split("__")
        tax_name = parts[-1] if len(parts) > 1 and parts[-1] else tax_key
    else:
        tax_name = tax_key
    return f"{level[0]}_{tax_name}"


def _aggregate_confidence_flags(flags: List[str]) -> str:
    if not flags:
        return "high"
    # Preserve a stable priority so the string always reports the
    # strongest caveat first.
    priority = {"no_support": 0, "low_support": 1, "reduced_fastmode": 2}
    unique = sorted(set(flags), key=lambda flag: priority.get(flag, 99))
    return ",".join(unique)
