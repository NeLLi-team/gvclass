"""
Full summarization module for GVClass results matching original output format.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
from collections import Counter

from src.config.marker_sets import (
    MIRUS_MODELS,
    BUSCO_MODELS,
    PHAGE_MODELS,
    GVOG4M_MODELS,
    GVOG8M_MODELS,
    UNI56_MODELS,
    MCP_MODELS,
    MRYA_MODELS,
)
from src.core.weighted_completeness import create_weighted_calculator

logger = logging.getLogger(__name__)


class FullSummarizer:
    """Generate complete summary matching original GVClass output."""

    def __init__(self, database_path: Path):
        """Initialize with database path."""
        self.database_path = database_path
        self.labels_file = database_path / "gvclassSeptember25_labels.tsv"
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
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        genome_id = parts[0]
                        # Parse taxonomy string: Domain|Phylum|Class|Order|Family|Genus|Species
                        tax_str = parts[1]
                        if "|" in tax_str:
                            tax_parts = tax_str.split("|")
                        else:
                            # Handle single value or malformed taxonomy
                            tax_parts = [tax_str]

                        # Pad with empty strings if not enough levels
                        while len(tax_parts) < 7:
                            tax_parts.append("")

                        # Store as list: [domain, phylum, class, order, family, genus, species]
                        labels[genome_id] = tax_parts
        except Exception as e:
            logger.error(f"Error loading labels: {e}")
        return labels

    def format_tax_level_counts(self, tax_counter: Counter) -> str:
        """Format taxonomy counts with percentages."""
        if not tax_counter:
            return ""

        total = sum(tax_counter.values())
        sorted_counts = sorted(tax_counter.items(), key=lambda x: x[1], reverse=True)

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
            # Load order completeness table
            comp_df = pd.read_csv(self.completeness_table, sep="\t")

            # Find the row for this order
            order_row = comp_df[comp_df["Order"] == order_tax]
            if order_row.empty:
                return 0.0, 0.0, 0.0, 0.0

            # Get orthogroups for this order
            ogs_str = order_row["Orthogroups"].values[0]
            order_ogs = [og.strip() for og in ogs_str.split(",")]

            # Load marker counts
            marker_counts = {}
            if counts_file.exists():
                with open(counts_file, "r") as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if len(parts) == 2:
                            marker_counts[parts[0]] = int(parts[1])

            # Calculate traditional metrics
            present_ogs = sum(1 for og in order_ogs if marker_counts.get(og, 0) > 0)
            total_hits = sum(marker_counts.get(og, 0) for og in order_ogs)

            completeness = (present_ogs / len(order_ogs)) * 100 if order_ogs else 0
            duplication = total_hits / present_ogs if present_ogs > 0 else 0

            # Calculate weighted completeness using ML-enhanced approach
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
                    f"(traditional: {completeness:.2f}%)"
                )
            except Exception as e:
                logger.warning(f"Weighted completeness calculation failed: {e}")
                weighted_completeness = completeness  # Fallback to traditional
                confidence_score = 50.0  # Moderate confidence for fallback

            return completeness, duplication, weighted_completeness, confidence_score

        except Exception as e:
            logger.error(f"Error calculating order metrics: {e}")
            return 0.0, 0.0, 0.0, 0.0

    def summarize_query_full(
        self,
        query_id: str,
        query_output_dir: Path,
        tree_nn_results: Dict[str, Dict[str, Dict[str, float]]],
        mode_fast: bool = True,
    ) -> Dict[str, any]:
        """Generate full summary for a query matching original format."""

        # First try to load from original stats file created during reformatting
        stats_tsv_file = query_output_dir / "stats" / f"{query_id}_stats.tsv"
        basic_stats = {}
        if stats_tsv_file.exists():
            try:
                stats_df = pd.read_csv(stats_tsv_file, sep="\t")
                if not stats_df.empty:
                    basic_stats = stats_df.iloc[0].to_dict()
            except Exception as e:
                logger.error(f"Error reading stats TSV: {e}")

        # Then load genetic code info from the stats.tab file
        stats_tab_file = query_output_dir / "stats" / f"{query_id}.stats.tab"
        if stats_tab_file.exists():
            try:
                with open(stats_tab_file, "r") as f:
                    for line in f:
                        if line.startswith("ttable\t"):
                            basic_stats["ttable"] = line.strip().split("\t")[1]
                        elif line.startswith("genes\t"):
                            basic_stats["genecount"] = int(line.strip().split("\t")[1])
                        elif line.startswith("coding_density\t"):
                            basic_stats["CODINGperc"] = float(
                                line.strip().split("\t")[1]
                            )
            except Exception as e:
                logger.error(f"Error reading stats tab: {e}")

        # Initialize result with proper stats
        result = {
            "query": query_id,
            "contigs": basic_stats.get("contigs", 0),
            "LENbp": basic_stats.get("LENbp", 0),
            "GCperc": basic_stats.get("GCperc", 0.0),
            "genecount": basic_stats.get("genecount", 0),
            "CODINGperc": basic_stats.get("CODINGperc", 0.0),
            "ttable": (
                basic_stats.get("ttable", "unknown")
                if basic_stats.get("ttable", "unknown") != "0"
                else "codemeta"
            ),
        }

        # Process tree nearest neighbors for full taxonomy
        tax_levels = [
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
        tax_level_mapping = {
            "domain": 0,  # Domain (e.g., NCLDV)
            "phylum": 1,  # Phylum (e.g., Nucleocytoviricota)
            "class": 2,  # Class (e.g., Megaviricetes)
            "order": 3,  # Order (e.g., Imitervirales)
            "family": 4,  # Family
            "genus": 5,  # Genus
            "species": 6,  # Species
        }

        # Collect all taxonomies by level
        tax_counters = {level: Counter() for level in tax_levels}
        distances = []

        # Process all nearest neighbors
        for marker, query_neighbors in tree_nn_results.items():
            for query_protein, neighbors in query_neighbors.items():
                for neighbor, distance in neighbors.items():
                    # Extract genome ID (handle missing | gracefully)
                    if "|" in neighbor:
                        genome_id = neighbor.split("|")[0]
                    else:
                        genome_id = neighbor

                    if genome_id in self.labels_dict:
                        tax_info = self.labels_dict[genome_id]

                        # Count taxonomies at each level
                        for level, idx in tax_level_mapping.items():
                            if idx < len(tax_info) and tax_info[idx].strip():
                                tax_value = tax_info[idx].strip()
                                # For domain, use the first part of genome_id if it has __
                                if level == "domain":
                                    if "__" in genome_id:
                                        domain_prefix = genome_id.split("__")[0]
                                    else:
                                        domain_prefix = tax_value
                                    tax_counters[level][domain_prefix] += 1
                                else:
                                    # For all other levels, use domain prefix + tax value
                                    if "__" in genome_id:
                                        domain_prefix = genome_id.split("__")[0]
                                        tax_counters[level][
                                            f"{domain_prefix}__{tax_value}"
                                        ] += 1
                                    else:
                                        tax_counters[level][tax_value] += 1

                        distances.append(distance)

        # Generate formatted counts for each level
        for level in tax_levels:
            result[level] = self.format_tax_level_counts(tax_counters[level])

        # Calculate average distance
        result["avgdist"] = sum(distances) / len(distances) if distances else 0.0

        # Generate majority and strict consensus taxonomies
        majority_parts = []
        strict_parts = []

        # Start with species and work up to domain
        for level in reversed(tax_levels):
            majority, strict = self.get_tax_consensus(tax_counters[level], level)
            majority_parts.append(majority)
            strict_parts.append(strict)

        result["taxonomy_majority"] = ";".join(reversed(majority_parts))
        result["taxonomy_strict"] = ";".join(reversed(strict_parts))

        # Calculate marker-based metrics
        counts_file = query_output_dir / "hmmout" / "models.counts"
        if counts_file.exists():
            marker_counts = {}
            with open(counts_file, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        marker_counts[parts[0]] = int(parts[1])

            # GVOG4 metrics
            gvog4_unique = sum(1 for m in GVOG4M_MODELS if marker_counts.get(m, 0) > 0)
            result["gvog4_unique"] = gvog4_unique

            # GVOG8 metrics
            gvog8_unique = sum(1 for m in GVOG8M_MODELS if marker_counts.get(m, 0) > 0)
            gvog8_total = sum(marker_counts.get(m, 0) for m in GVOG8M_MODELS)
            gvog8_dup = gvog8_total / gvog8_unique if gvog8_unique > 0 else 0
            result["gvog8_unique"] = gvog8_unique
            result["gvog8_total"] = gvog8_total
            result["gvog8_dup"] = gvog8_dup

            # MCP metrics
            mcp_total = sum(marker_counts.get(m, 0) for m in MCP_MODELS)
            result["mcp_total"] = mcp_total

            # MIRUS metrics
            mirus_unique = sum(1 for m in MIRUS_MODELS if marker_counts.get(m, 0) > 0)
            mirus_total = sum(marker_counts.get(m, 0) for m in MIRUS_MODELS)
            mirus_dup = mirus_total / mirus_unique if mirus_unique > 0 else 0
            result["mirus_unique"] = mirus_unique
            result["mirus_total"] = mirus_total
            result["mirus_dup"] = mirus_dup

            # MRYA metrics
            mrya_unique = sum(1 for m in MRYA_MODELS if marker_counts.get(m, 0) > 0)
            mrya_total = sum(marker_counts.get(m, 0) for m in MRYA_MODELS)
            result["mrya_unique"] = mrya_unique
            result["mrya_total"] = mrya_total

            # Phage metrics
            phage_unique = sum(1 for m in PHAGE_MODELS if marker_counts.get(m, 0) > 0)
            phage_total = sum(marker_counts.get(m, 0) for m in PHAGE_MODELS)
            result["phage_unique"] = phage_unique
            result["phage_total"] = phage_total

            # Cellular contamination (UNI56 + BUSCO)
            cellular_models = UNI56_MODELS + BUSCO_MODELS
            cellular_unique = sum(
                1 for m in cellular_models if marker_counts.get(m, 0) > 0
            )
            cellular_total = sum(marker_counts.get(m, 0) for m in cellular_models)
            cellular_dup = (
                cellular_total / cellular_unique if cellular_unique > 0 else 0
            )
            result["cellular_unique"] = cellular_unique
            result["cellular_total"] = cellular_total
            result["cellular_dup"] = cellular_dup

        else:
            # Set defaults if no counts file
            for metric in [
                "gvog4_unique",
                "gvog8_unique",
                "gvog8_total",
                "mcp_total",
                "mirus_unique",
                "mirus_total",
                "mrya_unique",
                "mrya_total",
                "phage_unique",
                "phage_total",
                "cellular_unique",
                "cellular_total",
            ]:
                result[metric] = 0
            for metric in ["gvog8_dup", "mirus_dup", "cellular_dup"]:
                result[metric] = 0.0

        # Calculate order-specific metrics if we have an order assignment
        order_tax = None
        if tax_counters["order"]:
            most_common_order = tax_counters["order"].most_common(1)[0][0]
            if "__" in most_common_order:
                parts = most_common_order.split("__")
                order_tax = parts[1] if len(parts) > 1 else most_common_order

        if order_tax:
            completeness, duplication, weighted_completeness, confidence_score = (
                self.calculate_order_metrics(counts_file, order_tax)
            )
            result["order_completeness"] = completeness
            result["order_dup"] = duplication
            result["order_weighted_completeness"] = weighted_completeness
            result["order_confidence_score"] = confidence_score
        else:
            result["order_completeness"] = 0.0
            result["order_dup"] = 0.0
            result["order_weighted_completeness"] = 0.0
            result["order_confidence_score"] = 0.0

        return result
