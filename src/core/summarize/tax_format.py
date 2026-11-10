"""
Taxonomy label, key and count formatting for the GVClass summary.

Stateless helpers that turn taxonomy counters and lineage keys into output
strings.
"""

from collections import Counter
from collections.abc import Iterable

# Order of the flags in ``taxonomy_confidence``.
_CONFIDENCE_FLAG_ORDER = (
    "no_support",
    "low_support",
    "tie",
    "mixed_cellular",
    "plurality",
    "reduced_fastmode",
)

_TAXONOMY_PREFIX_TO_LEVEL = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}


def parse_taxonomy_majority(taxonomy_majority: str) -> dict[str, str]:
    """Return unprefixed ranks from an accepted consensus lineage."""
    lineage = {level: "" for level in _TAXONOMY_PREFIX_TO_LEVEL.values()}
    for token in str(taxonomy_majority or "").split(";"):
        prefix, separator, value = token.strip().partition("_")
        level = _TAXONOMY_PREFIX_TO_LEVEL.get(prefix)
        if separator and level:
            lineage[level] = value
    return lineage


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


def _build_taxonomy_key(level: str, domain_prefix: str, tax_value: str) -> str:
    if level == "domain":
        return domain_prefix or tax_value
    if domain_prefix:
        return f"{domain_prefix}__{tax_value}"
    return tax_value


def _lineage_matches_prefix(
    lineage: tuple[str, ...], selected_prefix: tuple[str, ...]
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


def _aggregate_confidence_flags(flags: Iterable[str]) -> str:
    """Join confidence flags in their fixed order, or return ``high`` for none."""
    unique = sorted(set(flags), key=_CONFIDENCE_FLAG_ORDER.index)
    return ",".join(unique) or "high"
