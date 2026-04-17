"""
Full summarization module for GVClass results matching original output format.
"""

from collections import Counter, defaultdict
import logging
import math
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

from src.config.marker_sets import (
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
from src.core.marker_extraction import (
    count_unique_proteins_by_category_models,
    count_unique_proteins_by_category_prefixes,
    count_unique_proteins_for_markers,
    count_unique_proteins_for_prefixes,
    parse_hmm_output,
)
from src.core.weighted_completeness import create_weighted_calculator
from src.core.novelty_completeness import create_novelty_completeness_scorer
from src.core.contamination_scoring import create_contamination_scorer

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

    def __init__(
        self,
        database_path: Path,
        completeness_mode: str = "legacy",
        sensitive_mode: bool = False,
    ):
        """Initialize with database path.

        ``sensitive_mode`` is propagated to
        :class:`~src.core.contamination_scoring.ContaminationScorer` for
        diagnostics/logging. As of v1.4.3 the bundled contamination model
        is trained on sensitive-mode features
        (``training_profile: sensitive_mode_features`` in the model card),
        so the same trained-model path runs under both sensitive and
        non-sensitive settings.
        """
        self.database_path = database_path
        self.completeness_mode = completeness_mode
        self.sensitive_mode = sensitive_mode
        self.labels_file = database_path / "gvclassFeb26_labels.tsv"
        self.completeness_table = database_path / "order_completeness.tab"
        self.order_stats_df = self._load_order_stats()
        self.labels_dict = self.load_labels()

        # Initialize weighted completeness calculator
        self.weighted_calculator = create_weighted_calculator(database_path)
        self.novelty_scorer = create_novelty_completeness_scorer(database_path)
        self.contamination_scorer = create_contamination_scorer(
            database_path, sensitive_mode=sensitive_mode
        )

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

    def _load_order_stats(self) -> pd.DataFrame:
        """Load order-level marker recovery baselines used for completeness scaling."""
        try:
            return pd.read_csv(self.completeness_table, sep="\t")
        except Exception as e:
            logger.error(f"Error loading order completeness table: {e}")
            return pd.DataFrame()

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
        self,
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
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
        if marker_hits:
            category_counts = count_unique_proteins_by_category_prefixes(
                marker_hits, VP_CATEGORY_PREFIXES
            )
            categories_present = {
                category for category, count in category_counts.items() if count > 0
            }
            vp_mcp_count = category_counts.get("MCP", 0)
            total_vp_hits = sum(category_counts.values())
            plv_count = count_unique_proteins_for_prefixes(marker_hits, [PLV_PREFIX])
            completeness = len(categories_present)
            vp_df = total_vp_hits / 4.0 if total_vp_hits > 0 else 0.0
            return f"{completeness}/4", vp_mcp_count, plv_count, vp_df

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
        self,
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
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
        if marker_hits:
            category_counts = count_unique_proteins_by_category_models(
                marker_hits, MIRUS_CATEGORY_MODELS
            )
            categories_present = {
                category for category, count in category_counts.items() if count > 0
            }
            total_mirus_hits = sum(category_counts.values())
            completeness = len(categories_present)
            mirus_df = total_mirus_hits / 4.0 if total_mirus_hits > 0 else 0.0
            return f"{completeness}/4", mirus_df

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
        self, counts_file: Path, order_tax: str, family_tax: str = ""
    ) -> Dict[str, any]:
        """Calculate order-specific completeness and duplication metrics.

        Returns:
            Dictionary of raw and normalized order completeness metrics.
        """
        try:
            order_ogs = self._load_order_orthogroups(order_tax)
            if not order_ogs:
                return self._default_order_metrics(order_tax)

            marker_counts = self._load_marker_counts(counts_file)
            raw_completeness, duplication = self._calculate_traditional_order_metrics(
                marker_counts, order_ogs
            )
            raw_weighted_completeness, confidence_score = (
                self._calculate_weighted_order_metrics(
                    marker_counts, order_tax, order_ogs, raw_completeness
                )
            )
            normalized_completeness, baseline_mean, baseline_std = (
                self._normalize_order_completeness(raw_completeness, order_tax)
            )
            normalized_weighted_completeness, _, _ = (
                self._normalize_order_completeness(raw_weighted_completeness, order_tax)
            )
            result = {
                "order_completeness": normalized_completeness,
                "order_completeness_raw": raw_completeness,
                "order_dup": duplication,
                "order_weighted_completeness": normalized_weighted_completeness,
                "order_weighted_completeness_raw": raw_weighted_completeness,
                "order_confidence_score": confidence_score,
                "order_completeness_baseline_mean": baseline_mean,
                "order_completeness_baseline_std": baseline_std,
                "order_completeness_reference_order": order_tax,
                "order_completeness_strategy": "order_baseline_ratio_v1",
                "estimated_completeness": normalized_completeness,
                "estimated_completeness_strategy": "order_baseline_ratio_v1",
            }
            if self.novelty_scorer.available:
                result.update(
                    self.novelty_scorer.calculate(
                        marker_counts=marker_counts,
                        order_tax=order_tax,
                        family_tax=family_tax,
                        strategy1_raw=raw_completeness,
                        strategy1_score=normalized_completeness,
                        selected_mode=self.completeness_mode,
                    )
                )
            return result

        except Exception as e:
            logger.error(f"Error calculating order metrics: {e}")
            return self._default_order_metrics(order_tax)

    def _default_order_metrics(self, order_tax: str = "") -> Dict[str, any]:
        """Return default order completeness metrics."""
        return {
            "order_completeness": 0.0,
            "order_completeness_raw": 0.0,
            "order_dup": 0.0,
            "order_weighted_completeness": 0.0,
            "order_weighted_completeness_raw": 0.0,
            "order_confidence_score": 0.0,
            "order_completeness_baseline_mean": 0.0,
            "order_completeness_baseline_std": 0.0,
            "order_completeness_reference_order": order_tax,
            "order_completeness_strategy": "order_baseline_ratio_v1",
            "order_completeness_v2": 0.0,
            "order_completeness_v2_strategy": "novelty_aware_v1",
            "order_completeness_v2_strategy2_raw": 0.0,
            "order_completeness_v2_strategy2_normalized": 0.0,
            "order_completeness_v2_support_score": 0.0,
            "order_completeness_v2_ood_flag": "unassigned",
            "order_completeness_v2_reference_group": "unavailable",
            "order_completeness_v2_validation_mode": "unavailable",
            "order_completeness_v2_informative_fraction": 0.0,
            "estimated_completeness": 0.0,
            "estimated_completeness_strategy": "order_baseline_ratio_v1",
        }

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

    def _load_marker_hits(self, hmmout_file: Path) -> Dict[str, Set[str]]:
        """Load per-marker protein hits from a filtered HMM output file."""
        if not hmmout_file.exists():
            return {}

        try:
            return parse_hmm_output(hmmout_file)
        except Exception as exc:
            logger.error("Error reading marker hits from %s: %s", hmmout_file, exc)
            return {}

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

    def _normalize_order_completeness(
        self, raw_completeness: float, order_tax: str
    ) -> Tuple[float, float, float]:
        """Normalize marker recovery against complete-reference order baselines."""
        if self.order_stats_df.empty:
            return raw_completeness, 0.0, 0.0

        order_row = self.order_stats_df[self.order_stats_df["Order"] == order_tax]
        if order_row.empty:
            return raw_completeness, 0.0, 0.0

        baseline_mean = float(order_row["Average_Percent"].iloc[0])
        baseline_std = float(order_row["Std_Percent"].iloc[0])
        if baseline_mean <= 0:
            return raw_completeness, baseline_mean, baseline_std

        normalized = (raw_completeness / baseline_mean) * 100.0
        normalized = max(0.0, min(100.0, normalized))
        return normalized, baseline_mean, baseline_std

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
    ) -> Tuple[
        Dict[str, Counter],
        Dict[str, Dict[str, Counter]],
        List[float],
    ]:
        """Aggregate taxonomy counts from NN tree hits.

        Returns a tuple of:
        * ``tax_counters`` — flat per-level Counter (legacy shape, used by
          contamination scoring and the per-level detailed columns).
        * ``per_marker_counters`` — ``Dict[level, Dict[marker, Counter]]``
          that retains the per-marker breakdown. The per-marker majority
          rule in :meth:`_build_consensus_taxonomies` uses this to avoid
          the prior bug where a single marker with many paralog hits could
          dominate the majority vote by sheer count.
        * ``distances`` — flat list of patristic distances for the mean.
        """
        tax_counters = {level: Counter() for level in self.TAX_LEVELS}
        per_marker_counters: Dict[str, Dict[str, Counter]] = {
            level: defaultdict(Counter) for level in self.TAX_LEVELS
        }
        distances: List[float] = []

        for marker, query_neighbors in tree_nn_results.items():
            for neighbors in query_neighbors.values():
                for neighbor, distance in neighbors.items():
                    genome_id = self._extract_genome_id(neighbor)
                    tax_info = self.labels_dict.get(genome_id)
                    if not tax_info:
                        continue

                    self._add_taxonomy_counts_for_neighbor(
                        genome_id, tax_info, tax_counters
                    )
                    self._add_taxonomy_counts_for_marker(
                        genome_id, tax_info, per_marker_counters, marker
                    )
                    distances.append(distance)

        return tax_counters, per_marker_counters, distances

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

    def _add_taxonomy_counts_for_marker(
        self,
        genome_id: str,
        tax_info: List[str],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        marker: str,
    ) -> None:
        domain_prefix = genome_id.split("__")[0] if "__" in genome_id else ""
        for level, idx in self.TAX_LEVEL_MAPPING.items():
            if idx >= len(tax_info):
                continue
            tax_value = tax_info[idx].strip()
            if not tax_value:
                continue
            tax_key = self._build_taxonomy_key(level, domain_prefix, tax_value)
            per_marker_counters[level][marker][tax_key] += 1

    @staticmethod
    def _build_taxonomy_key(level: str, domain_prefix: str, tax_value: str) -> str:
        if level == "domain":
            return domain_prefix or tax_value
        if domain_prefix:
            return f"{domain_prefix}__{tax_value}"
        return tax_value

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

    @classmethod
    def _min_markers_for_level(cls, level: str, mode_fast: bool) -> int:
        threshold_map = cls._MIN_MARKERS_FOR_LEVEL_FAST if mode_fast else cls._MIN_MARKERS_FOR_LEVEL
        return threshold_map.get(level, 2)

    @staticmethod
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

    @classmethod
    def _per_marker_majority(
        cls,
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
            winner_pair = cls._stable_top(counter)
            if winner_pair is None:
                continue
            marker_votes[winner_pair[0]] += 1

        if not marker_votes:
            return None, 0, total_markers

        winning_taxon, supporting_markers = cls._stable_top(marker_votes)
        return winning_taxon, supporting_markers, total_markers

    def _build_consensus_taxonomies(
        self,
        tax_counters: Dict[str, Counter],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        mode_fast: bool,
    ) -> Tuple[str, str, str]:
        """Build majority and strict taxonomy strings and a confidence label.

        The majority call now uses a per-marker vote: one winning taxon per
        marker, then the across-marker winner. A minimum number of distinct
        supporting markers is required (see ``_min_markers_for_level``);
        otherwise the level is emitted unassigned. Strict consensus
        continues to require 100% agreement within the flat counter for
        backward compatibility with pre-1.4.3 callers.
        """
        majority_parts = []
        strict_parts = []
        confidence_flags: List[str] = []

        for level in reversed(self.TAX_LEVELS):
            flat_counter = tax_counters[level]
            per_marker_level = per_marker_counters.get(level, {})
            threshold = self._min_markers_for_level(level, mode_fast)
            winning, supporting, total = self._per_marker_majority(per_marker_level)

            if winning is not None and supporting >= threshold:
                majority = self._format_tax_label(level, winning)
                # Only flag reduced_fastmode when the relaxed threshold is
                # what actually allowed the call — i.e. the support falls
                # below the standard-mode threshold but clears the fast-mode
                # one. A fast-mode run that already has 3+ supporting
                # order-level markers should not be downgraded.
                standard_threshold = self._MIN_MARKERS_FOR_LEVEL.get(level, 2)
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

            # Strict keeps the flat-counter 100% rule.
            _, strict = self.get_tax_consensus(flat_counter, level)
            majority_parts.append(majority)
            strict_parts.append(strict)

        majority_str = ";".join(reversed(majority_parts))
        strict_str = ";".join(reversed(strict_parts))
        confidence = self._aggregate_confidence_flags(confidence_flags)
        return majority_str, strict_str, confidence

    @staticmethod
    def _format_tax_label(level: str, tax_key: str) -> str:
        if level == "domain":
            return f"{level[0]}_{tax_key}"
        if "__" in tax_key:
            parts = tax_key.split("__")
            tax_name = parts[-1] if len(parts) > 1 and parts[-1] else tax_key
        else:
            tax_name = tax_key
        return f"{level[0]}_{tax_name}"

    @staticmethod
    def _aggregate_confidence_flags(flags: List[str]) -> str:
        if not flags:
            return "high"
        # Preserve a stable priority so the string always reports the
        # strongest caveat first.
        priority = {"low_support": 0, "reduced_fastmode": 1}
        unique = sorted(set(flags), key=lambda flag: priority.get(flag, 99))
        return ",".join(unique)

    def _add_taxonomy_summary(
        self,
        result: Dict[str, any],
        tax_counters: Dict[str, Counter],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        distances: List[float],
        mode_fast: bool,
    ) -> None:
        """Populate taxonomy count and consensus fields on the result payload."""
        for level in self.TAX_LEVELS:
            result[level] = self.format_tax_level_counts(tax_counters[level])

        result["avgdist"] = sum(distances) / len(distances) if distances else 0.0
        taxonomy_majority, taxonomy_strict, taxonomy_confidence = (
            self._build_consensus_taxonomies(
                tax_counters, per_marker_counters, mode_fast
            )
        )
        result["taxonomy_majority"] = taxonomy_majority
        result["taxonomy_strict"] = taxonomy_strict
        result["taxonomy_confidence"] = taxonomy_confidence

    def _add_marker_metrics(
        self,
        result: Dict[str, any],
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
    ) -> None:
        """Populate marker-based metrics for the query."""
        result["gvog4_unique"] = sum(1 for m in GVOG4M_MODELS if marker_counts.get(m, 0) > 0)

        gvog8_unique = sum(1 for m in GVOG8M_MODELS if marker_counts.get(m, 0) > 0)
        gvog8_total = sum(marker_counts.get(m, 0) for m in GVOG8M_MODELS)
        result["gvog8_unique"] = gvog8_unique
        result["gvog8_total"] = gvog8_total
        result["gvog8_dup"] = gvog8_total / gvog8_unique if gvog8_unique > 0 else 0

        if marker_hits:
            result["ncldv_mcp_total"] = count_unique_proteins_for_markers(
                marker_hits, NCLDV_MCP_MODELS
            )
            result["mcp_total"] = count_unique_proteins_for_markers(
                marker_hits, MCP_MODELS
            )
            result["mrya_total"] = count_unique_proteins_for_markers(
                marker_hits, MRYA_MODELS
            )
            result["phage_total"] = count_unique_proteins_for_markers(
                marker_hits, PHAGE_MODELS
            )
            cellular_models = UNI56_MODELS + BUSCO_MODELS
            cellular_total = count_unique_proteins_for_markers(
                marker_hits, cellular_models
            )
        else:
            result["ncldv_mcp_total"] = sum(marker_counts.get(m, 0) for m in NCLDV_MCP_MODELS)
            result["mcp_total"] = sum(marker_counts.get(m, 0) for m in MCP_MODELS)
            result["mrya_total"] = sum(marker_counts.get(m, 0) for m in MRYA_MODELS)
            result["phage_total"] = sum(marker_counts.get(m, 0) for m in PHAGE_MODELS)
            cellular_models = UNI56_MODELS + BUSCO_MODELS
            cellular_total = sum(marker_counts.get(m, 0) for m in cellular_models)

        vp_completeness, vp_mcp, plv_count, vp_df = self.calculate_vp_metrics(
            marker_counts, marker_hits
        )
        result["vp_completeness"] = vp_completeness
        result["vp_mcp"] = vp_mcp
        result["plv"] = plv_count
        result["vp_df"] = vp_df

        mirus_completeness, mirus_df = self.calculate_mirus_completeness(
            marker_counts, marker_hits
        )
        result["mirus_completeness"] = mirus_completeness
        result["mirus_df"] = mirus_df

        mrya_unique = sum(1 for m in MRYA_MODELS if marker_counts.get(m, 0) > 0)
        result["mrya_unique"] = mrya_unique

        phage_unique = sum(1 for m in PHAGE_MODELS if marker_counts.get(m, 0) > 0)
        result["phage_unique"] = phage_unique

        cellular_unique = sum(1 for m in cellular_models if marker_counts.get(m, 0) > 0)
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

    def _extract_family_taxonomy(self, tax_counters: Dict[str, Counter]) -> str:
        """Extract most-supported family taxonomy token from counters."""
        if not tax_counters["family"]:
            return ""

        most_common_family = tax_counters["family"].most_common(1)[0][0]
        if "__" not in most_common_family:
            return ""
        parts = most_common_family.split("__")
        return parts[1] if len(parts) > 1 else most_common_family

    def _add_order_metrics(
        self, result: Dict[str, any], counts_file: Path, tax_counters: Dict[str, Counter]
    ) -> None:
        """Populate order-level completeness metrics from order assignment."""
        order_tax = self._extract_order_taxonomy(tax_counters)
        family_tax = self._extract_family_taxonomy(tax_counters)
        if order_tax:
            result.update(self.calculate_order_metrics(counts_file, order_tax, family_tax))
            return

        result.update(self._default_order_metrics(order_tax))

    def _add_contamination_metrics(
        self,
        result: Dict[str, any],
        query_id: str,
        query_output_dir: Path,
        marker_counts: Dict[str, int],
        tax_counters: Dict[str, Counter],
        tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None,
    ) -> None:
        """Populate contamination metrics (rule-based, then ML override if available)."""
        try:
            result.update(
                self.contamination_scorer.score_rule_based(
                    result=result,
                    marker_counts=marker_counts,
                    tax_counters=tax_counters,
                    query_output_dir=query_output_dir,
                )
            )
        except Exception as exc:
            logger.error("Error calculating contamination metrics: %s", exc)
            result.update(
                {
                    "contamination_score_v1": 0.0,
                    "contamination_flag_v1": "unknown",
                    "contamination_source_v1": "none",
                    "contamination_cellular_signal_v1": 0.0,
                    "contamination_phage_signal_v1": 0.0,
                    "contamination_duplication_signal_v1": 0.0,
                    "contamination_viral_mixture_signal_v1": 0.0,
                    "contamination_nonviral_hit_fraction_v1": 0.0,
                    "estimated_contamination": 0.0,
                    "estimated_contamination_strategy": "rule_based_v1",
                }
            )

        # Collect contig-level features and apply the trained ML model.
        try:
            primary_order = self._extract_order_taxonomy(tax_counters)
            primary_family = self._extract_family_taxonomy(tax_counters)
            query_fna = self._resolve_query_fna_path(query_id, query_output_dir)
            query_faa = self._resolve_query_faa_path(query_id, query_output_dir)
            contig_features_raw = self.contamination_scorer.collect_contig_features(
                query_fna,
                query_faa,
                query_output_dir / "blastp_out",
                primary_order,
                primary_family,
                tree_nn_results=tree_nn_results or {},
            )
            contig_features = {
                "suspicious_bp_fraction_v2": contig_features_raw["suspicious_bp_fraction"],
                "suspicious_contig_count_v2": contig_features_raw["suspicious_contig_count"],
                # v1.4.3 Phase 2 per-contig purity features.
                "cellular_coherent_contig_count": contig_features_raw.get(
                    "cellular_coherent_contig_count", 0
                ),
                "cellular_coherent_protein_fraction": contig_features_raw.get(
                    "cellular_coherent_protein_fraction", 0.0
                ),
                "cellular_coherent_bp_fraction": contig_features_raw.get(
                    "cellular_coherent_bp_fraction", 0.0
                ),
                "cellular_lineage_purity_median": contig_features_raw.get(
                    "cellular_lineage_purity_median", 0.0
                ),
                "cellular_hit_identity_median": contig_features_raw.get(
                    "cellular_hit_identity_median", 0.0
                ),
                "viral_bearing_contig_count": contig_features_raw.get(
                    "viral_bearing_contig_count", 0
                ),
                "contig_attribution_mode": contig_features_raw.get(
                    "contig_attribution_mode", "fna_gene_calling"
                ),
            }
            result.update(contig_features)

            # v1.4.3: the bundled contamination model is trained on
            # sensitive_mode=true features (see
            # ``src/bundled_models/contamination_model.yaml``
            # ``training_profile: sensitive_mode_features``). It is therefore
            # safe to apply under either sensitive or non-sensitive runs; no
            # gate is required here. The sensitive_mode flag is retained on
            # the class for test coverage and forward compatibility.
            _ = self.sensitive_mode

            if not self.contamination_scorer.ml_available:
                raise RuntimeError(
                    "A trained contamination model is required for production contamination estimates"
                )
            result.update(
                self.contamination_scorer.predict_contamination(result, contig_features)
            )
            result["contamination_type"] = self._classify_contamination_type(result)
            result["_contamination_reporting_threshold"] = self._contamination_reporting_threshold()
            result["_contamination_candidates"] = [
                {
                    **candidate,
                    "candidate_type": self._candidate_type_for_reason(candidate.get("reason", "")),
                }
                for candidate in contig_features_raw.get("suspicious_contigs", [])
            ]
        except Exception as exc:
            raise RuntimeError(
                f"Failed to compute contamination estimate with the trained model: {exc}"
            ) from exc

    def _contamination_reporting_threshold(self) -> float:
        threshold = float(self.contamination_scorer.ml_threshold or 0.0)
        return max(10.0, threshold)

    @staticmethod
    def _candidate_type_for_reason(reason: str) -> str:
        mapping = {
            "cellular_markers": "cellular",
            "cellular_hits": "cellular",
            "phage_markers": "phage",
            "phage_hits": "phage",
            "viral_mixture": "mixed_viral",
        }
        return mapping.get(str(reason), "uncertain")

    def _classify_contamination_type(self, result: Dict[str, any]) -> str:
        raw_estimate = result.get("estimated_contamination", 0.0)
        try:
            estimated = float(raw_estimate) if raw_estimate is not None else 0.0
        except (TypeError, ValueError):
            estimated = 0.0
        # NaN is not a sub-threshold ``clean`` bin — it reflects an
        # unavailable estimate (e.g. a future model-skip path, or a
        # scorer that returned no numeric prediction). Report it as
        # ``uncertain`` so downstream consumers treat it as missing
        # rather than clean.
        if math.isnan(estimated):
            return "uncertain"
        if estimated < self._contamination_reporting_threshold():
            return "clean"

        # v1.4.3 Phase 2: when the per-contig classifier is active and has
        # identified at least one cellular_coherent contig, the cellular
        # path takes precedence over the rule-based source. This is the
        # sharper per-contig signal promised in the plan — it fires on
        # genuine host-contig contamination (consistent cellular lineage
        # with no viral markers on the same contig) rather than on
        # scattered HGT-leaning tree placements.
        cellular_coherent = int(result.get("cellular_coherent_contig_count", 0) or 0)
        if cellular_coherent >= 1:
            return "cellular"

        source = str(result.get("contamination_source_v1", "none") or "none")

        # v1.4.3 Phase 2: downgrade `mixed_viral` to `uncertain` when the
        # per-contig classifier found zero cellular_coherent contigs AND
        # the bin carries a real viral signature (≥ 3 viral_bearing
        # contigs). That pattern is the fingerprint of a novel giant
        # virus whose markers scatter against non-coherent EUK neighbors
        # due to sparse references, rather than a bin that actually
        # carries multiple giant-virus orders. The absolute ML score is
        # still surfaced so users can curate.
        viral_bearing = int(result.get("viral_bearing_contig_count", 0) or 0)
        if source == "viral_mixture" and cellular_coherent == 0 and viral_bearing >= 3:
            return "uncertain"

        mapping = {
            "cellular": "cellular",
            "phage": "phage",
            "duplication": "duplication",
            "viral_mixture": "mixed_viral",
            "none": "uncertain",
        }
        return mapping.get(source, "uncertain")

    @staticmethod
    def _resolve_query_fna_path(query_id: str, query_output_dir: Path) -> Path:
        candidates = [
            query_output_dir / "query_fna" / f"{query_id}.fna",
            query_output_dir / f"{query_id}_reformatted.fna",
            query_output_dir / "gene_calling" / f"{query_id}_reformatted" / f"{query_id}_reformatted.fna",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        return Path()

    @staticmethod
    def _resolve_query_faa_path(query_id: str, query_output_dir: Path) -> Path:
        candidates = [
            query_output_dir / "query_faa" / f"{query_id}.faa",
            query_output_dir / "gene_calling" / f"{query_id}_reformatted" / f"{query_id}_reformatted.faa",
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate
        matches = sorted(query_output_dir.glob("gene_calling/**/*.faa"))
        return matches[0] if matches else Path()

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

        tax_counters, per_marker_counters, distances = self._collect_taxonomy_counts(
            tree_nn_results
        )
        self._add_taxonomy_summary(
            result, tax_counters, per_marker_counters, distances, mode_fast
        )

        counts_file = query_output_dir / "hmmout" / "models.counts"
        marker_hits_file = query_output_dir / "hmmout" / "models.out.filtered"
        marker_counts = self._load_marker_counts(counts_file)
        marker_hits = self._load_marker_hits(marker_hits_file)
        if marker_counts:
            self._add_marker_metrics(result, marker_counts, marker_hits)
        else:
            self._set_default_marker_metrics(result)

        self._add_order_metrics(result, counts_file, tax_counters)
        self._add_contamination_metrics(
            result,
            query_id,
            query_output_dir,
            marker_counts,
            tax_counters,
            tree_nn_results=tree_nn_results,
        )
        return result
