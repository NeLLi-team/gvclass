"""
Full summarization module for GVClass results matching original output format.
"""

from collections import Counter, defaultdict
import logging
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

import pandas as pd

from src.config.marker_sets import CAPS_MARKER_GROUP
from src.core.summarize import (
    consensus,
    labels_taxonomy,
    marker_panels,
    order_completeness,
    paths_stats,
    tax_format,
)
from src.core.weighted_completeness import create_weighted_calculator
from src.core.novelty_completeness import create_novelty_completeness_scorer
from src.core.contamination_scoring import create_contamination_scorer
from src.utils.resource_store import ResourceStore

logger = logging.getLogger(__name__)


# ``_dominant`` and the caps-group whitelist moved to
# src/core/summarize/labels_taxonomy.py together with the tree-NN capsid
# methods that use them; ``derive_capscan_group`` below still reads both, and
# the module-level names stay reachable for existing callers.
_dominant = labels_taxonomy._dominant
CAPS_GROUPS = labels_taxonomy.CAPS_GROUPS


def derive_capscan_group(marker_hits, marker_counts) -> str:
    """Dominant capscan (Bellas&Sommaruga 2026) group for a query.

    Resolves a single best group PER PROTEIN first (a protein cross-hitting several
    caps HMMs of one group votes once), then takes the dominant group across proteins.
    Avoids the profile-multiplicity bias of counting raw per-(protein,model) hits and
    the per-marker tree-vote fragmentation of GVClass's consensus. Falls back to raw
    caps-hit counts when per-protein hits are unavailable. Ties resolve to "".
    """
    caps_groups: Counter = Counter()
    if marker_hits:
        prot_groups: defaultdict = defaultdict(Counter)
        for marker, proteins in marker_hits.items():
            if "_caps_" in marker and marker in CAPS_MARKER_GROUP:
                grp = CAPS_MARKER_GROUP[marker]
                for protein in proteins:
                    prot_groups[protein][grp] += 1
        for grp_counts in prot_groups.values():
            winner = _dominant(grp_counts)
            if winner:
                caps_groups[winner] += 1
    if not caps_groups:
        for marker, count in marker_counts.items():
            if count > 0 and "_caps_" in marker and marker in CAPS_MARKER_GROUP:
                caps_groups[CAPS_MARKER_GROUP[marker]] += count
    return _dominant(caps_groups)


class FullSummarizer:
    """Generate complete summary matching original GVClass output."""

    TAX_LEVELS = consensus.TAX_LEVELS
    TAX_LEVEL_MAPPING = labels_taxonomy.TAX_LEVEL_MAPPING

    # Empty label index until ``__init__`` loads one. The wrappers below read
    # ``labels_dict`` before delegating, whereas the bodies they replaced read
    # it only on the branches that need it; the default keeps the no-label
    # branches (e.g. a PLV count with no marker tree) working as they did.
    labels_dict: Dict[str, List[str]] = {}

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
        self.labels_file = ResourceStore(database_path).label_path("labels.tsv")
        self.completeness_table = database_path / "markers" / "order_completeness.tab"
        self.order_stats_df = self._load_order_stats()
        self.labels_dict = self.load_labels()

        # Initialize weighted completeness calculator
        self.weighted_calculator = create_weighted_calculator(database_path)
        self.novelty_scorer = create_novelty_completeness_scorer(database_path)
        self.contamination_scorer = create_contamination_scorer(
            database_path, sensitive_mode=sensitive_mode
        )

    # The label index and the taxonomy counting built on it live in
    # src/core/summarize/labels_taxonomy.py; the methods below delegate so
    # existing callers keep working unchanged. ``labels_dict`` is read here and
    # passed in, so the module functions stay stateless.

    def load_labels(self) -> Dict[str, List[str]]:
        return labels_taxonomy.load_labels(self.labels_file)

    def _load_order_stats(self) -> pd.DataFrame:
        """Load order-level marker recovery baselines used for completeness scaling."""
        try:
            return pd.read_csv(self.completeness_table, sep="\t")
        except Exception as e:
            logger.error(f"Error loading order completeness table: {e}")
            return pd.DataFrame()

    def _parse_label_line(self, line: str) -> Optional[Tuple[str, List[str]]]:
        return labels_taxonomy._parse_label_line(line)

    def _plv_count_tree_aware(
        self,
        marker_hits: Optional[Dict[str, Set[str]]],
        tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]],
    ) -> int:
        return labels_taxonomy._plv_count_tree_aware(
            self.labels_dict, marker_hits, tree_nn_results
        )

    def _protein_nearest_domain(
        self,
        tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]],
    ) -> Dict[str, str]:
        return labels_taxonomy._protein_nearest_domain(
            self.labels_dict, tree_nn_results
        )

    # Marker panel metrics live in src/core/summarize/marker_panels.py; the
    # methods below delegate so existing callers keep working unchanged. The
    # panel functions never touch ``labels_dict``: the label-dependent values
    # (``_protein_nearest_domain``, ``_plv_count_tree_aware``,
    # ``_capsid_group_counter``) are computed here and passed in.

    @staticmethod
    def _panel_has_lineage_protein(
        proteins: Set[str],
        protein_domain: Optional[Dict[str, str]],
        expected_domain: str,
    ) -> bool:
        return marker_panels._panel_has_lineage_protein(
            proteins, protein_domain, expected_domain
        )

    def calculate_vp_metrics(
        self,
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
        tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None,
        protein_domain: Optional[Dict[str, str]] = None,
    ) -> Tuple[str, int, int, float]:
        # The tree-aware PLV count needs ``labels_dict``; it is computed here
        # under the same ``marker_hits`` guard as the branch that consumes it.
        plv_tree_count = (
            self._plv_count_tree_aware(marker_hits, tree_nn_results)
            if marker_hits
            else 0
        )
        return marker_panels.calculate_vp_metrics(
            marker_counts, marker_hits, protein_domain, plv_tree_count
        )

    def calculate_mirus_completeness(
        self,
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
        protein_domain: Optional[Dict[str, str]] = None,
    ) -> Tuple[str, float]:
        return marker_panels.calculate_mirus_completeness(
            marker_counts, marker_hits, protein_domain
        )

    def _extract_genome_id(self, protein_id: str) -> str:
        return labels_taxonomy._extract_genome_id(self.labels_dict, protein_id)

    # Taxonomy formatting lives in src/core/summarize/tax_format.py; the
    # methods below delegate so existing callers keep working unchanged.

    def format_tax_level_counts(
        self, tax_counter: Counter, level: str | None = None
    ) -> str:
        return tax_format.format_tax_level_counts(tax_counter, level)

    @staticmethod
    def _format_tax_count_key(tax_key: str, level: str | None = None) -> str:
        return tax_format._format_tax_count_key(tax_key, level)

    @staticmethod
    def _low_frequency_tax_group_key(tax_key: str, level: str | None = None) -> str:
        return tax_format._low_frequency_tax_group_key(tax_key, level)

    def get_tax_consensus(self, tax_counter: Counter, level: str) -> Tuple[str, str]:
        return tax_format.get_tax_consensus(tax_counter, level)

    # Order-level completeness lives in
    # src/core/summarize/order_completeness.py; the methods below delegate so
    # existing callers keep working unchanged. The reference tables and scorers
    # are read here and passed in, so the module functions stay stateless.

    def calculate_order_metrics(
        self, counts_file: Path, order_tax: str, family_tax: str = ""
    ) -> Dict[str, Any]:
        return order_completeness.calculate_order_metrics(
            self.completeness_table,
            self.order_stats_df,
            self.weighted_calculator,
            self.novelty_scorer,
            self.completeness_mode,
            counts_file,
            order_tax,
            family_tax,
        )

    def _default_order_metrics(self, order_tax: str = "") -> Dict[str, Any]:
        return order_completeness._default_order_metrics(order_tax)

    def _load_order_orthogroups(self, order_tax: str) -> List[str]:
        return order_completeness._load_order_orthogroups(
            self.completeness_table, order_tax
        )

    # File/path helpers live in src/core/summarize/paths_stats.py; the methods
    # below delegate so existing callers keep working unchanged.

    def _load_marker_counts(self, counts_file: Path) -> Dict[str, int]:
        return paths_stats._load_marker_counts(counts_file)

    def _load_marker_hits(self, hmmout_file: Path) -> Dict[str, Set[str]]:
        return paths_stats._load_marker_hits(hmmout_file)

    def _load_marker_scores(self, hmmout_file: Path) -> Dict[str, Dict[str, float]]:
        return paths_stats._load_marker_scores(hmmout_file)

    def _calculate_traditional_order_metrics(
        self, marker_counts: Dict[str, int], order_ogs: List[str]
    ) -> Tuple[float, float]:
        return order_completeness._calculate_traditional_order_metrics(
            marker_counts, order_ogs
        )

    def _calculate_weighted_order_metrics(
        self,
        marker_counts: Dict[str, int],
        order_tax: str,
        order_ogs: List[str],
        fallback_completeness: float,
    ) -> Tuple[float, float]:
        return order_completeness._calculate_weighted_order_metrics(
            self.weighted_calculator,
            marker_counts,
            order_tax,
            order_ogs,
            fallback_completeness,
        )

    def _normalize_order_completeness(
        self, raw_completeness: float, order_tax: str
    ) -> Tuple[float, float, float]:
        return order_completeness._normalize_order_completeness(
            self.order_stats_df, raw_completeness, order_tax
        )

    def _read_stats_tsv(self, stats_tsv_file: Path) -> Dict[str, Any]:
        return paths_stats._read_stats_tsv(stats_tsv_file)

    def _read_stats_tab(
        self, stats_tab_file: Path, basic_stats: Dict[str, Any]
    ) -> None:
        paths_stats._read_stats_tab(stats_tab_file, basic_stats)

    def _update_basic_stats_from_tab_line(
        self, line: str, basic_stats: Dict[str, Any]
    ) -> None:
        paths_stats._update_basic_stats_from_tab_line(line, basic_stats)

    def _load_basic_stats(
        self, query_id: str, query_output_dir: Path
    ) -> Dict[str, Any]:
        return paths_stats._load_basic_stats(query_id, query_output_dir)

    def _initialize_summary_result(
        self, query_id: str, basic_stats: Dict[str, Any]
    ) -> Dict[str, Any]:
        return paths_stats._initialize_summary_result(query_id, basic_stats)

    def _collect_taxonomy_counts(
        self, tree_nn_results: Dict[str, Dict[str, Dict[str, float]]]
    ) -> Tuple[
        Dict[str, Counter],
        Dict[str, Dict[str, Counter]],
        Dict[str, Counter[Tuple[str, ...]]],
        List[float],
    ]:
        return labels_taxonomy._collect_taxonomy_counts(
            self.labels_dict, tree_nn_results
        )

    def _build_taxonomy_lineage(
        self, genome_id: str, tax_info: List[str]
    ) -> Tuple[str, ...]:
        return labels_taxonomy._build_taxonomy_lineage(genome_id, tax_info)

    def _add_taxonomy_counts_for_neighbor(
        self, genome_id: str, tax_info: List[str], tax_counters: Dict[str, Counter]
    ) -> None:
        labels_taxonomy._add_taxonomy_counts_for_neighbor(
            genome_id, tax_info, tax_counters
        )

    def _add_taxonomy_counts_for_marker(
        self,
        genome_id: str,
        tax_info: List[str],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        marker: str,
    ) -> None:
        labels_taxonomy._add_taxonomy_counts_for_marker(
            genome_id, tax_info, per_marker_counters, marker
        )

    @staticmethod
    def _build_taxonomy_key(level: str, domain_prefix: str, tax_value: str) -> str:
        return tax_format._build_taxonomy_key(level, domain_prefix, tax_value)

    # Taxonomy consensus lives in src/core/summarize/consensus.py; the
    # thresholds and methods below delegate so existing callers keep working
    # unchanged.

    _MIN_MARKERS_FOR_LEVEL: Dict[str, int] = consensus._MIN_MARKERS_FOR_LEVEL
    _MIN_MARKERS_FOR_LEVEL_FAST: Dict[str, int] = consensus._MIN_MARKERS_FOR_LEVEL_FAST

    @classmethod
    def _min_markers_for_level(cls, level: str, mode_fast: bool) -> int:
        return consensus._min_markers_for_level(level, mode_fast)

    @staticmethod
    def _stable_top(counter: Counter) -> Optional[Tuple[str, int]]:
        return tax_format._stable_top(counter)

    @classmethod
    def _per_marker_majority(
        cls,
        per_marker_level_counters: Dict[str, Counter],
    ) -> Tuple[Optional[str], int, int]:
        return consensus._per_marker_majority(per_marker_level_counters)

    @staticmethod
    def _tax_key_matches_domain(tax_key: str, domain: Optional[str]) -> bool:
        return tax_format._tax_key_matches_domain(tax_key, domain)

    @classmethod
    def _filter_counter_to_domain(
        cls, counter: Counter, domain: Optional[str]
    ) -> Counter:
        return consensus._filter_counter_to_domain(counter, domain)

    @classmethod
    def _filter_per_marker_counters_to_domain(
        cls,
        per_marker_level_counters: Dict[str, Counter],
        domain: Optional[str],
    ) -> Dict[str, Counter]:
        return consensus._filter_per_marker_counters_to_domain(
            per_marker_level_counters, domain
        )

    def _build_consensus_taxonomies(
        self,
        tax_counters: Dict[str, Counter],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        mode_fast: bool,
        per_marker_lineage_counters: Optional[
            Dict[str, Counter[Tuple[str, ...]]]
        ] = None,
    ) -> Tuple[str, str, str]:
        return consensus._build_consensus_taxonomies(
            tax_counters, per_marker_counters, mode_fast, per_marker_lineage_counters
        )

    @staticmethod
    def _lineage_matches_prefix(
        lineage: Tuple[str, ...], selected_prefix: Tuple[str, ...]
    ) -> bool:
        return tax_format._lineage_matches_prefix(lineage, selected_prefix)

    @classmethod
    def _lineage_rank_counter(
        cls,
        per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
        rank_index: int,
        selected_prefix: Tuple[str, ...],
    ) -> Counter:
        return consensus._lineage_rank_counter(
            per_marker_lineage_counters, rank_index, selected_prefix
        )

    @classmethod
    def _per_marker_lineage_majority(
        cls,
        per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
        rank_index: int,
        selected_prefix: Tuple[str, ...],
    ) -> Tuple[Optional[str], int, int]:
        return consensus._per_marker_lineage_majority(
            per_marker_lineage_counters, rank_index, selected_prefix
        )

    def _format_lineage_tax_label(
        self, level: str, selected_prefix: Tuple[str, ...], tax_value: str
    ) -> str:
        return consensus._format_lineage_tax_label(level, selected_prefix, tax_value)

    def _build_consensus_taxonomies_from_lineages(
        self,
        per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
        mode_fast: bool,
    ) -> Tuple[str, str, str]:
        return consensus._build_consensus_taxonomies_from_lineages(
            per_marker_lineage_counters, mode_fast
        )

    @staticmethod
    def _format_tax_label(level: str, tax_key: str) -> str:
        return tax_format._format_tax_label(level, tax_key)

    @staticmethod
    def _aggregate_confidence_flags(flags: List[str]) -> str:
        return tax_format._aggregate_confidence_flags(flags)

    def _add_taxonomy_summary(
        self,
        result: Dict[str, Any],
        tax_counters: Dict[str, Counter],
        per_marker_counters: Dict[str, Dict[str, Counter]],
        per_marker_lineage_counters: Dict[str, Counter[Tuple[str, ...]]],
        distances: List[float],
        mode_fast: bool,
    ) -> None:
        """Populate taxonomy count and consensus fields on the result payload."""
        for level in self.TAX_LEVELS:
            result[level] = self.format_tax_level_counts(
                tax_counters[level], level=level
            )

        result["avgdist"] = sum(distances) / len(distances) if distances else 0.0
        taxonomy_majority, _taxonomy_strict, taxonomy_confidence = (
            self._build_consensus_taxonomies(
                tax_counters,
                per_marker_counters,
                mode_fast,
                per_marker_lineage_counters,
            )
        )
        result["taxonomy_majority"] = taxonomy_majority
        result["taxonomy_confidence"] = taxonomy_confidence

    def _caps_group_from_tree(self, tree_nn_results) -> str:
        return labels_taxonomy._caps_group_from_tree(self.labels_dict, tree_nn_results)

    def _capsid_group_counter(self, tree_nn_results) -> str:
        return labels_taxonomy._capsid_group_counter(self.labels_dict, tree_nn_results)

    def _add_marker_metrics(
        self,
        result: Dict[str, Any],
        marker_counts: Dict[str, int],
        marker_hits: Optional[Dict[str, Set[str]]] = None,
        tree_nn_results: Optional[Dict[str, Dict[str, Dict[str, float]]]] = None,
        marker_scores: Optional[Dict[str, Dict[str, float]]] = None,
    ) -> None:
        """Populate marker-based metrics for the query."""
        # Per-protein nearest-neighbour domain (from the marker trees) buckets
        # cross-group panel hits: a VP/Mirus marker on a protein that actually
        # places with NCLDV is not credited to that panel. None on a counts-only
        # resume, where vp/mirus fall back to raw category presence.
        protein_domain = (
            self._protein_nearest_domain(tree_nn_results) if tree_nn_results else None
        )
        # The tree-aware PLV count and the capsid tally need ``labels_dict``, so
        # they are computed here and passed in, each under the same guard as the
        # branch that consumed it.
        plv_tree_count = (
            self._plv_count_tree_aware(marker_hits, tree_nn_results)
            if marker_hits
            else 0
        )
        capsid_group = (
            self._capsid_group_counter(tree_nn_results)
            if tree_nn_results
            and any(m in tree_nn_results for m in ("mcp_ncldv", "mcp_mirus", "mcp_plv"))
            else ""
        )
        result.update(
            marker_panels._add_marker_metrics(
                marker_counts,
                marker_hits,
                marker_scores,
                protein_domain,
                plv_tree_count,
                capsid_group,
            )
        )

    def _set_default_marker_metrics(self, result: Dict[str, Any]) -> None:
        marker_panels._set_default_marker_metrics(result)

    def _extract_order_taxonomy(self, tax_counters: Dict[str, Counter]) -> str:
        return order_completeness._extract_order_taxonomy(tax_counters)

    def _extract_family_taxonomy(self, tax_counters: Dict[str, Counter]) -> str:
        return order_completeness._extract_family_taxonomy(tax_counters)

    def _add_order_metrics(
        self,
        result: Dict[str, Any],
        counts_file: Path,
        tax_counters: Dict[str, Counter],
    ) -> None:
        order_completeness._add_order_metrics(
            self.completeness_table,
            self.order_stats_df,
            self.weighted_calculator,
            self.novelty_scorer,
            self.completeness_mode,
            result,
            counts_file,
            tax_counters,
        )

    def _add_contamination_metrics(
        self,
        result: Dict[str, Any],
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
                "suspicious_bp_fraction_v2": contig_features_raw[
                    "suspicious_bp_fraction"
                ],
                "suspicious_contig_count_v2": contig_features_raw[
                    "suspicious_contig_count"
                ],
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
            # ``resources/contamination/model.yaml``
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
            result["_contamination_reporting_threshold"] = (
                self._contamination_reporting_threshold()
            )
            result["_contamination_candidates"] = [
                {
                    **candidate,
                    "candidate_type": self._candidate_type_for_reason(
                        candidate.get("reason", "")
                    ),
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
        return paths_stats._candidate_type_for_reason(reason)

    def _classify_contamination_type(self, result: Dict[str, Any]) -> str:
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
        return paths_stats._resolve_query_fna_path(query_id, query_output_dir)

    @staticmethod
    def _resolve_query_faa_path(query_id: str, query_output_dir: Path) -> Path:
        return paths_stats._resolve_query_faa_path(query_id, query_output_dir)

    def summarize_query_full(
        self,
        query_id: str,
        query_output_dir: Path,
        tree_nn_results: Dict[str, Dict[str, Dict[str, float]]],
        mode_fast: bool = True,
    ) -> Dict[str, Any]:
        """Generate full summary for a query matching original format."""
        basic_stats = self._load_basic_stats(query_id, query_output_dir)
        result = self._initialize_summary_result(query_id, basic_stats)

        (
            tax_counters,
            per_marker_counters,
            per_marker_lineage_counters,
            distances,
        ) = self._collect_taxonomy_counts(tree_nn_results)
        self._add_taxonomy_summary(
            result,
            tax_counters,
            per_marker_counters,
            per_marker_lineage_counters,
            distances,
            mode_fast,
        )

        counts_file = query_output_dir / "hmmout" / "models.counts"
        marker_hits_file = query_output_dir / "hmmout" / "models.out.filtered"
        marker_counts = self._load_marker_counts(counts_file)
        marker_hits = self._load_marker_hits(marker_hits_file)
        marker_scores = self._load_marker_scores(marker_hits_file)
        if marker_counts:
            self._add_marker_metrics(
                result,
                marker_counts,
                marker_hits,
                tree_nn_results=tree_nn_results,
                marker_scores=marker_scores,
            )
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
