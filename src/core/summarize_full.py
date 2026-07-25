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
# functions that use them; ``derive_capscan_group`` below still reads both.
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
    """Generate complete summary matching original GVClass output.

    The stateless steps live in src/core/summarize/ (paths_stats, tax_format,
    consensus, marker_panels, labels_taxonomy, order_completeness); the methods
    below hold the orchestration and the instance state those functions are
    called with.
    """

    # Empty label index until ``__init__`` loads one. ``_add_marker_metrics``
    # reads it before delegating, so a ``__new__``-constructed instance keeps
    # the no-label branches (e.g. a PLV count with no marker tree) working.
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
        self.labels_dict = labels_taxonomy.load_labels(self.labels_file)

        # Initialize weighted completeness calculator
        self.weighted_calculator = create_weighted_calculator(database_path)
        self.novelty_scorer = create_novelty_completeness_scorer(database_path)
        self.contamination_scorer = create_contamination_scorer(
            database_path, sensitive_mode=sensitive_mode
        )

    def _load_order_stats(self) -> pd.DataFrame:
        """Load order-level marker recovery baselines used for completeness scaling."""
        try:
            return pd.read_csv(self.completeness_table, sep="\t")
        except Exception as e:
            logger.error(f"Error loading order completeness table: {e}")
            return pd.DataFrame()

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
        for level in consensus.TAX_LEVELS:
            result[level] = tax_format.format_tax_level_counts(
                tax_counters[level], level=level
            )

        result["avgdist"] = sum(distances) / len(distances) if distances else 0.0
        taxonomy_majority, _taxonomy_strict, taxonomy_confidence = (
            consensus._build_consensus_taxonomies(
                tax_counters,
                per_marker_counters,
                mode_fast,
                per_marker_lineage_counters,
            )
        )
        result["taxonomy_majority"] = taxonomy_majority
        result["taxonomy_confidence"] = taxonomy_confidence

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
            labels_taxonomy._protein_nearest_domain(self.labels_dict, tree_nn_results)
            if tree_nn_results
            else None
        )
        # The tree-aware PLV count and the capsid tally need ``labels_dict``, so
        # they are computed here and passed in, each under the same guard as the
        # branch that consumed it.
        plv_tree_count = (
            labels_taxonomy._plv_count_tree_aware(
                self.labels_dict, marker_hits, tree_nn_results
            )
            if marker_hits
            else 0
        )
        capsid_group = (
            labels_taxonomy._capsid_group_counter(self.labels_dict, tree_nn_results)
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
            primary_order = order_completeness._extract_order_taxonomy(tax_counters)
            primary_family = order_completeness._extract_family_taxonomy(tax_counters)
            query_fna = paths_stats._resolve_query_fna_path(query_id, query_output_dir)
            query_faa = paths_stats._resolve_query_faa_path(query_id, query_output_dir)
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
                    "candidate_type": paths_stats._candidate_type_for_reason(
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

    def summarize_query_full(
        self,
        query_id: str,
        query_output_dir: Path,
        tree_nn_results: Dict[str, Dict[str, Dict[str, float]]],
        mode_fast: bool = True,
    ) -> Dict[str, Any]:
        """Generate full summary for a query matching original format."""
        basic_stats = paths_stats._load_basic_stats(query_id, query_output_dir)
        result = paths_stats._initialize_summary_result(query_id, basic_stats)

        (
            tax_counters,
            per_marker_counters,
            per_marker_lineage_counters,
            distances,
        ) = labels_taxonomy._collect_taxonomy_counts(self.labels_dict, tree_nn_results)
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
        marker_counts = paths_stats._load_marker_counts(counts_file)
        marker_hits = paths_stats._load_marker_hits(marker_hits_file)
        marker_scores = paths_stats._load_marker_scores(marker_hits_file)
        if marker_counts:
            self._add_marker_metrics(
                result,
                marker_counts,
                marker_hits,
                tree_nn_results=tree_nn_results,
                marker_scores=marker_scores,
            )
        else:
            marker_panels._set_default_marker_metrics(result)

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
        self._add_contamination_metrics(
            result,
            query_id,
            query_output_dir,
            marker_counts,
            tax_counters,
            tree_nn_results=tree_nn_results,
        )
        return result
