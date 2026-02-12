"""
Full summarization module for GVClass results matching original output format.
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
from collections import Counter

from src.config.marker_sets import (
    MIRUS_MODELS,
    MIRUS_CATEGORY_MODELS,
    BUSCO_MODELS,
    PHAGE_MODELS,
    GVOG4M_MODELS,
    GVOG8M_MODELS,
    UNI56_MODELS,
    MCP_MODELS,
    NCLDV_MCP_MODELS,
    MRYA_MODELS,
    VP_CATEGORY_PREFIXES,
    PLV_PREFIX,
)
from src.core.weighted_completeness import create_weighted_calculator

logger = logging.getLogger(__name__)


class FullSummarizer:
    """Generate complete summary matching original GVClass output."""

    TAX_LEVELS = [
        "domain",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    ]
    TAX_LEVEL_MAPPING = {
        "domain": 0,
        "phylum": 1,
        "class": 2,
        "order": 3,
        "family": 4,
        "genus": 5,
        "species": 6,
    }

    def __init__(self, database_path: Path):
        """Initialize with database path."""
        self.database_path = database_path
        self.labels_file = database_path / "gvclassFeb26_labels.tsv"
        self.completeness_table = database_path / "order_completeness.tab"
        self.labels_dict = self.load_labels()

        # Initialize weighted completeness calculator
        self.weighted_calculator = create_weighted_calculator(database_path)

    def load_labels(self) -> Dict[str, List[str]]:
        """Load taxonomy labels from file."""
        labels = {}
        try:
            with open(self.labels_file, "r") as f:
                for line in f:
                    parsed = self._parse_label_line(line)
                    if parsed is None:
                        continue
                    genome_id, tax_parts = parsed
                    labels[genome_id] = tax_parts
        except Exception as e:
            logger.error(f"Error loading labels: {e}")
        return labels

    def _parse_label_line(self, line: str) -> Optional[Tuple[str, List[str]]]:
        if line.startswith("#"):
            return None

        parts = line.strip().split("\t")
        if len(parts) < 2:
            return None

        genome_id = parts[0]
        tax_parts = parts[1].split("|") if "|" in parts[1] else [parts[1]]
        while len(tax_parts) < 7:
            tax_parts.append("")
        return genome_id, tax_parts

    def calculate_vp_metrics(
        self, marker_counts: Dict[str, int]
    ) -> Tuple[str, int, int, float]:
        """Calculate VP (Virophage) completeness metrics.

        VP has 4 core categories: MCP, Penton, ATPase, Protease.
        Uses prefix matching (e.g., VP_MCP_* matches MCP category).

        Returns:
            Tuple of (completeness_str, vp_mcp_count, plv_count, vp_df)
            - completeness_str: "n/4" format
            - vp_mcp_count: count of unique proteins hitting VP_MCP markers
            - plv_count: count of unique proteins hitting PLV markers
            - vp_df: duplication factor (total VP hits / 4)
        """
        categories_present = set()
        vp_mcp_count = 0
        total_vp_hits = 0
        plv_count = 0

        for marker, count in marker_counts.items():
            if count <= 0:
                continue

            # Check VP categories by prefix
            for category, prefix in VP_CATEGORY_PREFIXES.items():
                if marker.startswith(prefix):
                    categories_present.add(category)
                    total_vp_hits += count
                    if category == "MCP":
                        vp_mcp_count += count

            # Check PLV markers
            if marker.startswith(PLV_PREFIX):
                plv_count += count

        completeness = len(categories_present)
        vp_df = total_vp_hits / 4.0 if total_vp_hits > 0 else 0.0

        return f"{completeness}/4", vp_mcp_count, plv_count, vp_df

    def calculate_mirus_completeness(
        self, marker_counts: Dict[str, int]
    ) -> Tuple[str, float]:
        """Calculate Mirusviricota completeness metrics.

        Mirus has 4 core categories:
        - MCP: Mirus_MCP, Mirus_JellyRoll
        - ATPase: Mirus_Terminase_ATPase, Mirus_Terminase_merged
        - Portal: Mirus_Portal
        - Triplex: Mirus_Triplex1, Mirus_Triplex2

        Returns:
            Tuple of (completeness_str, mirus_df)
            - completeness_str: "n/4" format
            - mirus_df: duplication factor (total Mirus hits / 4)
        """
        categories_present = set()
        total_mirus_hits = 0

        for marker, count in marker_counts.items():
            if count <= 0:
                continue

            # Check each Mirus category
            for category, models in MIRUS_CATEGORY_MODELS.items():
                if marker in models:
                    categories_present.add(category)
                    total_mirus_hits += count

        completeness = len(categories_present)
        mirus_df = total_mirus_hits / 4.0 if total_mirus_hits > 0 else 0.0

        return f"{completeness}/4", mirus_df

    def _extract_genome_id(self, protein_id: str) -> str:
        """
        Extract genome ID from a protein ID for labels lookup.

        Protein IDs can have formats like:
        - VP__IMGVR_UViG_3300044959|000235_9 -> VP__IMGVR_UViG_3300044959|000235
        - NCLDV__GCA_000123456|contig_1_42 -> NCLDV__GCA_000123456|contig_1

        The protein suffix is typically _N where N is a number at the end.

        Args:
            protein_id: Full protein identifier from tree

        Returns:
            Genome/contig ID suitable for labels lookup
        """
        # First try: check if the full ID (minus trailing _digits) is in labels
        # This handles: VP__IMGVR_UViG_3300044959|000235_9 -> VP__IMGVR_UViG_3300044959|000235
        stripped = re.sub(r"_\d+$", "", protein_id)
        if stripped in self.labels_dict:
            return stripped

        # Second try: just the part before | (for older format entries)
        # This handles: genome_id|protein_id -> genome_id
        if "|" in protein_id:
            base_id = protein_id.split("|")[0]
            if base_id in self.labels_dict:
                return base_id

        # Third try: strip protein suffix from full ID even if not in labels
        return stripped

    def format_tax_level_counts(self, tax_counter: Counter) -> str:
        """Format taxonomy counts with percentages, grouping low-frequency taxa."""
        if not tax_counter:
            return ""

        total = sum(tax_counter.values())
        if total <= 0:
            return ""

        grouped_counter = Counter(tax_counter)
        low_threshold = 5.0
        low_counts_by_prefix = Counter()

        for tax, count in tax_counter.items():
            percentage = (count / total) * 100
            if percentage <= low_threshold:
                if "__" in tax:
                    prefix = tax.split("__", 1)[0]
                elif "_" in tax:
                    prefix = tax.split("_", 1)[0]
                else:
                    prefix = tax
                low_counts_by_prefix[prefix] += count
                grouped_counter.pop(tax, None)

        for prefix, count in low_counts_by_prefix.items():
            grouped_counter[f"{prefix}_other"] += count

        sorted_counts = sorted(
            grouped_counter.items(), key=lambda x: x[1], reverse=True
        )

        formatted = []
        for tax, count in sorted_counts:
            percentage = (count / total) * 100
            formatted.append(f"{tax}:{count}({percentage:.2f}%)")

        return ",".join(formatted)

    def get_tax_consensus(self, tax_counter: Counter, level: str) -> Tuple[str, str]:
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

    def calculate_order_metrics(
        self, counts_file: Path, order_tax: str
    ) -> Tuple[float, float, float, float]:
        """Calculate order-specific completeness and duplication metrics.

        Returns:
            Tuple of (completeness, duplication, weighted_completeness, confidence_score)
        """
        try:
            order_ogs = self._load_order_orthogroups(order_tax)
            if not order_ogs:
                return 0.0, 0.0, 0.0, 0.0

            marker_counts = self._load_marker_counts(counts_file)
            completeness, duplication = self._calculate_traditional_order_metrics(
                marker_counts, order_ogs
            )
            weighted_completeness, confidence_score = (
                self._calculate_weighted_order_metrics(
                    marker_counts, order_tax, order_ogs, completeness
                )
            )
            return completeness, duplication, weighted_completeness, confidence_score

        except Exception as e:
            logger.error(f"Error calculating order metrics: {e}")
            return 0.0, 0.0, 0.0, 0.0

    def _load_order_orthogroups(self, order_tax: str) -> List[str]:
        """Load expected orthogroups for a specific taxonomic order."""
        comp_df = pd.read_csv(self.completeness_table, sep="\t")
        order_row = comp_df[comp_df["Order"] == order_tax]
        if order_row.empty:
            return []
        ogs_str = str(order_row["Orthogroups"].values[0])
        return [og.strip() for og in ogs_str.split(",") if og.strip()]

    def _load_marker_counts(self, counts_file: Path) -> Dict[str, int]:
        """Load marker counts from a models.counts file."""
        marker_counts = {}
        if not counts_file.exists():
            return marker_counts

        try:
            with open(counts_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) != 2:
                        continue
                    try:
                        marker_counts[parts[0]] = int(parts[1])
                    except ValueError:
                        logger.warning(
                            "Skipping malformed marker count '%s' in %s",
                            line.strip(),
                            counts_file,
                        )
        except Exception as e:
            logger.error(f"Error reading marker counts from {counts_file}: {e}")
        return marker_counts

    def _calculate_traditional_order_metrics(
        self, marker_counts: Dict[str, int], order_ogs: List[str]
    ) -> Tuple[float, float]:
        """Calculate completeness and duplication from expected order markers."""
        present_ogs = sum(1 for og in order_ogs if marker_counts.get(og, 0) > 0)
        total_hits = sum(marker_counts.get(og, 0) for og in order_ogs)

        completeness = (present_ogs / len(order_ogs)) * 100 if order_ogs else 0.0
        duplication = total_hits / present_ogs if present_ogs > 0 else 0.0
        return completeness, duplication

    def _calculate_weighted_order_metrics(
        self,
        marker_counts: Dict[str, int],
        order_tax: str,
        order_ogs: List[str],
        fallback_completeness: float,
    ) -> Tuple[float, float]:
        """Calculate ML-weighted completeness with fallback to traditional metrics."""
        try:
            weighted_completeness, confidence_score, _ = (
                self.weighted_calculator.calculate_weighted_completeness(
                    marker_counts=marker_counts,
                    taxonomic_order=order_tax,
                    expected_markers=order_ogs,
                )
            )
            logger.debug(
                f"Weighted completeness for {order_tax}: {weighted_completeness:.2f}% "
                f"(traditional: {fallback_completeness:.2f}%)"
            )
            return weighted_completeness, confidence_score
        except Exception as e:
            logger.warning(f"Weighted completeness calculation failed: {e}")
            return fallback_completeness, 50.0

    def _read_stats_tsv(self, stats_tsv_file: Path) -> Dict[str, any]:
        """Read base stats from the reformat stage TSV."""
        basic_stats = {}
        if not stats_tsv_file.exists():
            return basic_stats

        try:
            stats_df = pd.read_csv(stats_tsv_file, sep="\t")
            if not stats_df.empty:
                basic_stats = stats_df.iloc[0].to_dict()
        except Exception as e:
            logger.error(f"Error reading stats TSV: {e}")
        return basic_stats

    def _read_stats_tab(self, stats_tab_file: Path, basic_stats: Dict[str, any]) -> None:
        """Merge coding table/genecount fields from legacy stats.tab file."""
        if not stats_tab_file.exists():
            return

        try:
            with open(stats_tab_file, "r") as f:
                for line in f:
                    self._update_basic_stats_from_tab_line(line, basic_stats)
        except Exception as e:
            logger.error(f"Error reading stats tab: {e}")

    def _update_basic_stats_from_tab_line(
        self, line: str, basic_stats: Dict[str, any]
    ) -> None:
        if line.startswith("ttable\t"):
            basic_stats["ttable"] = line.strip().split("\t")[1]
            return
        if line.startswith("genes\t"):
            basic_stats["genecount"] = int(line.strip().split("\t")[1])
            return
        if line.startswith("coding_density\t"):
            basic_stats["CODINGperc"] = float(line.strip().split("\t")[1])

    def _load_basic_stats(self, query_id: str, query_output_dir: Path) -> Dict[str, any]:
        """Load query summary statistics from available stats files."""
        stats_tsv_file = query_output_dir / "stats" / f"{query_id}_stats.tsv"
        stats_tab_file = query_output_dir / "stats" / f"{query_id}.stats.tab"

        basic_stats = self._read_stats_tsv(stats_tsv_file)
        self._read_stats_tab(stats_tab_file, basic_stats)
        return basic_stats

    def _initialize_summary_result(
        self, query_id: str, basic_stats: Dict[str, any]
    ) -> Dict[str, any]:
        """Create the base result dictionary with normalized defaults."""
        ttable = basic_stats.get("ttable", "unknown")
        return {
            "query": query_id,
            "contigs": basic_stats.get("contigs", 0),
            "LENbp": basic_stats.get("LENbp", 0),
            "GCperc": basic_stats.get("GCperc", 0.0),
            "genecount": basic_stats.get("genecount", 0),
            "CODINGperc": basic_stats.get("CODINGperc", 0.0),
            "ttable": "codemeta" if ttable == "0" else ttable,
        }

    def _collect_taxonomy_counts(
        self, tree_nn_results: Dict[str, Dict[str, Dict[str, float]]]
    ) -> Tuple[Dict[str, Counter], List[float]]:
        """Aggregate taxonomy counts and distance statistics from NN tree hits."""
        tax_counters = {level: Counter() for level in self.TAX_LEVELS}
        distances = []

        for query_neighbors in tree_nn_results.values():
            for neighbors in query_neighbors.values():
                for neighbor, distance in neighbors.items():
                    genome_id = self._extract_genome_id(neighbor)
                    tax_info = self.labels_dict.get(genome_id)
                    if not tax_info:
                        continue

                    self._add_taxonomy_counts_for_neighbor(
                        genome_id, tax_info, tax_counters
                    )
                    distances.append(distance)

        return tax_counters, distances

    def _add_taxonomy_counts_for_neighbor(
        self, genome_id: str, tax_info: List[str], tax_counters: Dict[str, Counter]
    ) -> None:
        domain_prefix = genome_id.split("__")[0] if "__" in genome_id else ""
        for level, idx in self.TAX_LEVEL_MAPPING.items():
            if idx >= len(tax_info):
                continue
            tax_value = tax_info[idx].strip()
            if not tax_value:
                continue
            tax_key = self._build_taxonomy_key(level, domain_prefix, tax_value)
            tax_counters[level][tax_key] += 1

    @staticmethod
    def _build_taxonomy_key(level: str, domain_prefix: str, tax_value: str) -> str:
        if level == "domain":
            return domain_prefix or tax_value
        if domain_prefix:
            return f"{domain_prefix}__{tax_value}"
        return tax_value

    def _build_consensus_taxonomies(
        self, tax_counters: Dict[str, Counter]
    ) -> Tuple[str, str]:
        """Build majority and strict taxonomy strings from per-level counters."""
        majority_parts = []
        strict_parts = []

        for level in reversed(self.TAX_LEVELS):
            majority, strict = self.get_tax_consensus(tax_counters[level], level)
            majority_parts.append(majority)
            strict_parts.append(strict)

        return ";".join(reversed(majority_parts)), ";".join(reversed(strict_parts))

    def _add_taxonomy_summary(
        self,
        result: Dict[str, any],
        tax_counters: Dict[str, Counter],
        distances: List[float],
    ) -> None:
        """Populate taxonomy count and consensus fields on the result payload."""
        for level in self.TAX_LEVELS:
            result[level] = self.format_tax_level_counts(tax_counters[level])

        result["avgdist"] = sum(distances) / len(distances) if distances else 0.0
        taxonomy_majority, taxonomy_strict = self._build_consensus_taxonomies(
            tax_counters
        )
        result["taxonomy_majority"] = taxonomy_majority
        result["taxonomy_strict"] = taxonomy_strict

    def _add_marker_metrics(
        self, result: Dict[str, any], marker_counts: Dict[str, int]
    ) -> None:
        """Populate marker-based metrics for the query."""
        result["gvog4_unique"] = sum(1 for m in GVOG4M_MODELS if marker_counts.get(m, 0) > 0)

        gvog8_unique = sum(1 for m in GVOG8M_MODELS if marker_counts.get(m, 0) > 0)
        gvog8_total = sum(marker_counts.get(m, 0) for m in GVOG8M_MODELS)
        result["gvog8_unique"] = gvog8_unique
        result["gvog8_total"] = gvog8_total
        result["gvog8_dup"] = gvog8_total / gvog8_unique if gvog8_unique > 0 else 0

        result["ncldv_mcp_total"] = sum(marker_counts.get(m, 0) for m in NCLDV_MCP_MODELS)
        result["mcp_total"] = sum(marker_counts.get(m, 0) for m in MCP_MODELS)

        vp_completeness, vp_mcp, plv_count, vp_df = self.calculate_vp_metrics(marker_counts)
        result["vp_completeness"] = vp_completeness
        result["vp_mcp"] = vp_mcp
        result["plv"] = plv_count
        result["vp_df"] = vp_df

        mirus_completeness, mirus_df = self.calculate_mirus_completeness(marker_counts)
        result["mirus_completeness"] = mirus_completeness
        result["mirus_df"] = mirus_df

        mrya_unique = sum(1 for m in MRYA_MODELS if marker_counts.get(m, 0) > 0)
        result["mrya_unique"] = mrya_unique
        result["mrya_total"] = sum(marker_counts.get(m, 0) for m in MRYA_MODELS)

        phage_unique = sum(1 for m in PHAGE_MODELS if marker_counts.get(m, 0) > 0)
        result["phage_unique"] = phage_unique
        result["phage_total"] = sum(marker_counts.get(m, 0) for m in PHAGE_MODELS)

        cellular_models = UNI56_MODELS + BUSCO_MODELS
        cellular_unique = sum(1 for m in cellular_models if marker_counts.get(m, 0) > 0)
        cellular_total = sum(marker_counts.get(m, 0) for m in cellular_models)
        result["cellular_unique"] = cellular_unique
        result["cellular_total"] = cellular_total
        result["cellular_dup"] = cellular_total / cellular_unique if cellular_unique > 0 else 0

    def _set_default_marker_metrics(self, result: Dict[str, any]) -> None:
        """Apply default marker metric values when counts are unavailable."""
        for metric in [
            "gvog4_unique",
            "gvog8_unique",
            "gvog8_total",
            "ncldv_mcp_total",
            "mcp_total",
            "vp_mcp",
            "plv",
            "mrya_unique",
            "mrya_total",
            "phage_unique",
            "phage_total",
            "cellular_unique",
            "cellular_total",
        ]:
            result[metric] = 0

        for metric in ["gvog8_dup", "vp_df", "mirus_df", "cellular_dup"]:
            result[metric] = 0.0

        result["vp_completeness"] = "0/4"
        result["mirus_completeness"] = "0/4"

    def _extract_order_taxonomy(self, tax_counters: Dict[str, Counter]) -> str:
        """Extract most-supported order taxonomy token from counters."""
        if not tax_counters["order"]:
            return ""

        most_common_order = tax_counters["order"].most_common(1)[0][0]
        if "__" not in most_common_order:
            return ""
        parts = most_common_order.split("__")
        return parts[1] if len(parts) > 1 else most_common_order

    def _add_order_metrics(
        self, result: Dict[str, any], counts_file: Path, tax_counters: Dict[str, Counter]
    ) -> None:
        """Populate order-level completeness metrics from order assignment."""
        order_tax = self._extract_order_taxonomy(tax_counters)
        if order_tax:
            completeness, duplication, weighted_completeness, confidence_score = (
                self.calculate_order_metrics(counts_file, order_tax)
            )
            result["order_completeness"] = completeness
            result["order_dup"] = duplication
            result["order_weighted_completeness"] = weighted_completeness
            result["order_confidence_score"] = confidence_score
            return

        result["order_completeness"] = 0.0
        result["order_dup"] = 0.0
        result["order_weighted_completeness"] = 0.0
        result["order_confidence_score"] = 0.0

    def summarize_query_full(
        self,
        query_id: str,
        query_output_dir: Path,
        tree_nn_results: Dict[str, Dict[str, Dict[str, float]]],
        mode_fast: bool = True,
    ) -> Dict[str, any]:
        """Generate full summary for a query matching original format."""
        basic_stats = self._load_basic_stats(query_id, query_output_dir)
        result = self._initialize_summary_result(query_id, basic_stats)

        tax_counters, distances = self._collect_taxonomy_counts(tree_nn_results)
        self._add_taxonomy_summary(result, tax_counters, distances)

        counts_file = query_output_dir / "hmmout" / "models.counts"
        marker_counts = self._load_marker_counts(counts_file)
        if marker_counts:
            self._add_marker_metrics(result, marker_counts)
        else:
            self._set_default_marker_metrics(result)

        self._add_order_metrics(result, counts_file, tax_counters)
        return result
