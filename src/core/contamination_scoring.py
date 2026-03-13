"""Contamination scoring utilities for GVClass.

Method 1 is a production genome-level rule-based scorer.
Method 2 helpers expose contig-aware features for offline prototyping.
"""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List

from Bio import SeqIO

import numpy as np

from src.config.marker_sets import BUSCO_MODELS, PHAGE_MODELS, UNI56_MODELS
from src.core.marker_extraction import (
    count_unique_proteins_for_markers,
    parse_hmm_output,
)

logger = logging.getLogger(__name__)

CELLULAR_PREFIXES = {"BAC", "ARC", "EUK", "MITO", "PLASTID"}
PHAGE_PREFIXES = {"PHAGE", "VP", "PLV"}
GIANT_VIRUS_PREFIXES = {"NCLDV", "MIRUS", "MRYA"}
CELLULAR_MODELS = set(BUSCO_MODELS + UNI56_MODELS)
PHAGE_MODELS_SET = set(PHAGE_MODELS)

CONTAMINATION_MODEL_FILE = "contamination_model.joblib"
CONTAMINATION_MODEL_FEATURES = [
    "contamination_score_v1",
    "contamination_cellular_signal_v1",
    "contamination_phage_signal_v1",
    "contamination_duplication_signal_v1",
    "contamination_viral_mixture_signal_v1",
    "contamination_nonviral_hit_fraction_v1",
    "order_dup",
    "gvog8_dup",
    "cellular_unique",
    "cellular_total",
    "phage_unique",
    "phage_total",
    "contigs",
    "suspicious_bp_fraction_v2",
    "suspicious_contig_count_v2",
    "estimated_completeness",
]


@dataclass(frozen=True)
class BlastHit:
    query_id: str
    subject_id: str
    identity: float
    score: float
    subject_domain: str
    subject_order: str
    subject_family: str


class ContaminationScorer:
    """Runtime rule-based contamination scorer with reusable feature extractors."""

    def __init__(self, database_path: Path):
        self.database_path = database_path
        self.labels_file = database_path / "gvclassFeb26_labels.tsv"
        self.labels = self._load_labels()
        self.ml_model = None
        self.ml_model_name = "hist_gbm"
        self.ml_threshold = 0.0
        self.ml_available = False
        self._load_ml_model()

    def _load_labels(self) -> Dict[str, Dict[str, str]]:
        labels: Dict[str, Dict[str, str]] = {}
        try:
            with self.labels_file.open() as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 2:
                        continue
                    genome_id = parts[0]
                    tax = parts[1].split("|")
                    while len(tax) < 7:
                        tax.append("")
                    labels[genome_id] = {
                        "domain": genome_id.split("__", 1)[0],
                        "phylum": tax[1],
                        "class": tax[2],
                        "order": tax[3],
                        "family": tax[4],
                        "genus": tax[5],
                        "species": tax[6],
                    }
        except Exception as exc:
            logger.error("Failed to load contamination labels: %s", exc)
        return labels

    def _load_ml_model(self) -> None:
        bundled_model_path = (
            Path(__file__).resolve().parents[1] / "bundled_models" / CONTAMINATION_MODEL_FILE
        )
        candidate_paths = [bundled_model_path, self.database_path / CONTAMINATION_MODEL_FILE]
        for model_path in candidate_paths:
            if not model_path.exists():
                continue
            try:
                import joblib

                bundle = joblib.load(model_path)
                self.ml_model = bundle["model"]
                self.ml_model_name = str(bundle.get("model_name", "hist_gbm"))
                self.ml_threshold = bundle.get("threshold", 0.0)
                self.ml_available = True
                logger.info("Loaded ML contamination model from %s", model_path)
                return
            except Exception as exc:
                logger.warning(
                    "Failed to load ML contamination model from %s: %s",
                    model_path,
                    exc,
                )

    def predict_contamination(
        self, result: Dict[str, Any], contig_features: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Predict contamination using the ML model."""
        feature_values = []
        for feat in CONTAMINATION_MODEL_FEATURES:
            if feat in contig_features:
                feature_values.append(float(contig_features[feat]))
            else:
                feature_values.append(float(result.get(feat, 0.0) or 0.0))
        X = np.array(feature_values, dtype=np.float64).reshape(1, -1)
        prediction = float(np.clip(self.ml_model.predict(X)[0], 0.0, 100.0))
        return {
            "estimated_contamination": round(prediction, 2),
            "estimated_contamination_strategy": f"{self.ml_model_name}_v1",
        }

    @staticmethod
    def _extract_tax_token(counter: Counter) -> str:
        if not counter:
            return ""
        value = counter.most_common(1)[0][0]
        if "__" not in value:
            return ""
        parts = value.split("__")
        return parts[1] if len(parts) > 1 else value

    @staticmethod
    def _parse_subject_genome_id(subject_id: str) -> str:
        if "|" in subject_id:
            return subject_id.split("|", 1)[0]
        return subject_id

    @staticmethod
    def protein_to_contig_id(protein_id: str) -> str:
        token = protein_id.split("|", 1)[1] if "|" in protein_id else protein_id
        return token.rsplit("_", 1)[0] if "_" in token else token

    def _subject_taxonomy(self, subject_id: str) -> Dict[str, str]:
        genome_id = self._parse_subject_genome_id(subject_id)
        return self.labels.get(
            genome_id,
            {
                "domain": genome_id.split("__", 1)[0] if "__" in genome_id else "NA",
                "order": "",
                "family": "",
            },
        )

    def collect_best_blast_hits(self, blast_dir: Path) -> Dict[str, BlastHit]:
        best_hits: Dict[str, BlastHit] = {}
        if not blast_dir.exists():
            return best_hits

        for blast_file in sorted(blast_dir.glob("*.m8")):
            try:
                with blast_file.open() as handle:
                    for line in handle:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 12:
                            continue
                        query_id = parts[0]
                        subject_id = parts[1]
                        try:
                            identity = float(parts[2])
                            score = float(parts[11])
                        except ValueError:
                            continue
                        tax = self._subject_taxonomy(subject_id)
                        hit = BlastHit(
                            query_id=query_id,
                            subject_id=subject_id,
                            identity=identity,
                            score=score,
                            subject_domain=tax.get("domain", "NA"),
                            subject_order=tax.get("order", ""),
                            subject_family=tax.get("family", ""),
                        )
                        current = best_hits.get(query_id)
                        if current is None or (hit.score, hit.identity) > (current.score, current.identity):
                            best_hits[query_id] = hit
            except Exception as exc:
                logger.warning("Failed to parse blast file %s: %s", blast_file, exc)
        return best_hits

    @staticmethod
    def _majority_fraction(counter: Counter) -> float:
        if not counter:
            return 0.0
        total = sum(counter.values())
        if total == 0:
            return 0.0
        return counter.most_common(1)[0][1] / total

    @staticmethod
    def _second_fraction(counter: Counter) -> float:
        if len(counter) < 2:
            return 0.0
        total = sum(counter.values())
        if total == 0:
            return 0.0
        return counter.most_common(2)[1][1] / total

    def compute_blast_features(
        self,
        best_hits: Dict[str, BlastHit],
        primary_order: str,
        primary_family: str,
    ) -> Dict[str, float | int | str]:
        strong_hits = [hit for hit in best_hits.values() if hit.identity >= 60.0 and hit.score >= 80.0]
        if not strong_hits:
            return {
                "strong_best_hit_count": 0,
                "nonviral_best_hit_count": 0,
                "nonviral_best_hit_fraction": 0.0,
                "strong_cellular_hit_count": 0,
                "strong_phage_hit_count": 0,
                "foreign_viral_order_hit_count": 0,
                "foreign_viral_family_hit_count": 0,
                "dominant_nonviral_lineage_fraction": 0.0,
                "dominant_foreign_viral_order_fraction": 0.0,
            }

        nonviral = []
        cellular = []
        phage = []
        foreign_viral_order = []
        foreign_viral_family = []
        nonviral_lineages = Counter()
        foreign_orders = Counter()

        for hit in strong_hits:
            domain = hit.subject_domain
            if domain in CELLULAR_PREFIXES:
                cellular.append(hit)
                nonviral.append(hit)
                nonviral_lineages[hit.subject_order or domain] += 1
            elif domain in PHAGE_PREFIXES:
                phage.append(hit)
                nonviral.append(hit)
                nonviral_lineages[hit.subject_family or domain] += 1
            elif domain in GIANT_VIRUS_PREFIXES:
                if primary_order and hit.subject_order and hit.subject_order != primary_order:
                    foreign_viral_order.append(hit)
                    foreign_orders[hit.subject_order] += 1
                if primary_family and hit.subject_family and hit.subject_family != primary_family:
                    foreign_viral_family.append(hit)

        strong_count = len(strong_hits)
        nonviral_count = len(nonviral)
        nonviral_fraction = nonviral_count / strong_count
        dominant_nonviral_fraction = 0.0
        if nonviral_lineages and nonviral_count > 0:
            dominant_nonviral_fraction = nonviral_lineages.most_common(1)[0][1] / nonviral_count
        dominant_foreign_order_fraction = 0.0
        if foreign_orders and foreign_viral_order:
            dominant_foreign_order_fraction = foreign_orders.most_common(1)[0][1] / len(foreign_viral_order)

        return {
            "strong_best_hit_count": strong_count,
            "nonviral_best_hit_count": nonviral_count,
            "nonviral_best_hit_fraction": round(nonviral_fraction * 100.0, 2),
            "strong_cellular_hit_count": len(cellular),
            "strong_phage_hit_count": len(phage),
            "foreign_viral_order_hit_count": len(foreign_viral_order),
            "foreign_viral_family_hit_count": len(foreign_viral_family),
            "dominant_nonviral_lineage_fraction": round(dominant_nonviral_fraction * 100.0, 2),
            "dominant_foreign_viral_order_fraction": round(dominant_foreign_order_fraction * 100.0, 2),
        }

    def compute_tree_features(self, tax_counters: Dict[str, Counter]) -> Dict[str, float]:
        order_majority = self._majority_fraction(tax_counters.get("order", Counter()))
        family_majority = self._majority_fraction(tax_counters.get("family", Counter()))
        order_second = self._second_fraction(tax_counters.get("order", Counter()))
        family_second = self._second_fraction(tax_counters.get("family", Counter()))
        return {
            "order_majority_fraction": round(order_majority * 100.0, 2),
            "family_majority_fraction": round(family_majority * 100.0, 2),
            "order_secondary_fraction": round(order_second * 100.0, 2),
            "family_secondary_fraction": round(family_second * 100.0, 2),
        }

    @staticmethod
    def _clip(value: float) -> float:
        return max(0.0, min(100.0, value))

    def score_rule_based(
        self,
        result: Dict[str, Any],
        marker_counts: Dict[str, int],
        tax_counters: Dict[str, Counter],
        query_output_dir: Path,
    ) -> Dict[str, Any]:
        marker_hits = {}
        hmmout_file = query_output_dir / "hmmout" / "models.out.filtered"
        if hmmout_file.exists():
            try:
                marker_hits = parse_hmm_output(hmmout_file)
            except Exception as exc:
                logger.warning("Failed to parse marker hits from %s: %s", hmmout_file, exc)

        cellular_unique = sum(1 for m in CELLULAR_MODELS if marker_counts.get(m, 0) > 0)
        phage_unique = sum(1 for m in PHAGE_MODELS_SET if marker_counts.get(m, 0) > 0)
        if marker_hits:
            cellular_total = count_unique_proteins_for_markers(
                marker_hits, CELLULAR_MODELS
            )
            phage_total = count_unique_proteins_for_markers(
                marker_hits, PHAGE_MODELS_SET
            )
        else:
            cellular_total = sum(marker_counts.get(m, 0) for m in CELLULAR_MODELS)
            phage_total = sum(marker_counts.get(m, 0) for m in PHAGE_MODELS_SET)
        order_dup = float(result.get("order_dup", 0.0) or 0.0)
        gvog8_dup = float(result.get("gvog8_dup", 0.0) or 0.0)

        primary_order = self._extract_tax_token(tax_counters.get("order", Counter()))
        primary_family = self._extract_tax_token(tax_counters.get("family", Counter()))
        best_hits = self.collect_best_blast_hits(query_output_dir / "blastp_out")
        blast_features = self.compute_blast_features(best_hits, primary_order, primary_family)
        tree_features = self.compute_tree_features(tax_counters)

        cellular_signal = 0.0
        if cellular_unique >= 1:
            cellular_signal += max(0.0, (cellular_unique - 1) * 9.0)
            cellular_signal += max(0.0, cellular_total - cellular_unique) * 1.5
        cellular_signal += float(blast_features["strong_cellular_hit_count"]) * 3.0
        cellular_signal += float(blast_features["nonviral_best_hit_fraction"]) * 0.20
        if float(blast_features["strong_cellular_hit_count"]) >= 3:
            cellular_signal += float(blast_features["dominant_nonviral_lineage_fraction"]) * 0.10
        cellular_signal = self._clip(cellular_signal)

        phage_signal = 0.0
        if phage_unique >= 2:
            phage_signal += (phage_unique - 1) * 10.0
        phage_signal += max(0.0, phage_total - phage_unique) * 2.0
        phage_signal += float(blast_features["strong_phage_hit_count"]) * 4.0
        phage_signal = self._clip(phage_signal)

        duplication_signal = self._clip(max(0.0, order_dup - 1.20) * 45.0 + max(0.0, gvog8_dup - 1.30) * 30.0)

        viral_mixture_signal = 0.0
        viral_mixture_signal += max(0.0, float(tree_features["order_secondary_fraction"]) - 10.0) * 1.3
        viral_mixture_signal += max(0.0, float(tree_features["family_secondary_fraction"]) - 15.0) * 0.8
        viral_mixture_signal += float(blast_features["foreign_viral_order_hit_count"]) * 5.0
        viral_mixture_signal += float(blast_features["foreign_viral_family_hit_count"]) * 2.5
        if float(blast_features["foreign_viral_order_hit_count"]) >= 3:
            viral_mixture_signal += float(blast_features["dominant_foreign_viral_order_fraction"]) * 0.15
        viral_mixture_signal = self._clip(viral_mixture_signal)

        contamination_score = self._clip(
            (0.38 * cellular_signal)
            + (0.12 * phage_signal)
            + (0.20 * duplication_signal)
            + (0.30 * viral_mixture_signal)
        )

        if contamination_score >= 45.0:
            flag = "high"
        elif contamination_score >= 20.0:
            flag = "moderate"
        elif contamination_score >= 8.0:
            flag = "low"
        else:
            flag = "clean"

        components = {
            "cellular": cellular_signal,
            "phage": phage_signal,
            "duplication": duplication_signal,
            "viral_mixture": viral_mixture_signal,
        }
        source = max(components.items(), key=lambda item: item[1])[0]
        if contamination_score < 8.0:
            source = "none"

        return {
            "contamination_score_v1": round(contamination_score, 2),
            "contamination_flag_v1": flag,
            "contamination_source_v1": source,
            "contamination_cellular_signal_v1": round(cellular_signal, 2),
            "contamination_phage_signal_v1": round(phage_signal, 2),
            "contamination_duplication_signal_v1": round(duplication_signal, 2),
            "contamination_viral_mixture_signal_v1": round(viral_mixture_signal, 2),
            "contamination_nonviral_hit_fraction_v1": float(blast_features["nonviral_best_hit_fraction"]),
            "estimated_contamination": round(contamination_score, 2),
            "estimated_contamination_strategy": "rule_based_v1",
        }

    def collect_contig_features(
        self,
        query_fna: Path,
        query_faa: Path,
        blast_dir: Path,
        primary_order: str,
        primary_family: str,
    ) -> Dict[str, Any]:
        contig_lengths: Dict[str, int] = {}
        if query_fna and query_fna.exists() and query_fna.is_file():
            for record in SeqIO.parse(str(query_fna), "fasta"):
                contig_lengths[record.id] = len(record.seq)

        protein_to_contig: Dict[str, str] = {}
        if query_faa.exists():
            for record in SeqIO.parse(str(query_faa), "fasta"):
                contig_id = self.protein_to_contig_id(record.id)
                protein_to_contig[record.id] = contig_id
                if contig_id not in contig_lengths:
                    contig_lengths[contig_id] = 0
                contig_lengths[contig_id] += len(record.seq) * 3

        contig_markers: Dict[str, set[str]] = defaultdict(set)
        contig_cellular_markers: Dict[str, set[str]] = defaultdict(set)
        contig_phage_markers: Dict[str, set[str]] = defaultdict(set)
        best_hits: Dict[str, BlastHit] = {}

        for blast_file in sorted(blast_dir.glob("*.m8")):
            marker = blast_file.stem
            try:
                with blast_file.open() as handle:
                    for line in handle:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 12:
                            continue
                        query_id = parts[0]
                        subject_id = parts[1]
                        try:
                            identity = float(parts[2])
                            score = float(parts[11])
                        except ValueError:
                            continue
                        contig_id = protein_to_contig.get(query_id)
                        if contig_id:
                            contig_markers[contig_id].add(marker)
                            if marker in CELLULAR_MODELS:
                                contig_cellular_markers[contig_id].add(marker)
                            if marker in PHAGE_MODELS_SET:
                                contig_phage_markers[contig_id].add(marker)
                        tax = self._subject_taxonomy(subject_id)
                        hit = BlastHit(
                            query_id=query_id,
                            subject_id=subject_id,
                            identity=identity,
                            score=score,
                            subject_domain=tax.get("domain", "NA"),
                            subject_order=tax.get("order", ""),
                            subject_family=tax.get("family", ""),
                        )
                        current = best_hits.get(query_id)
                        if current is None or (hit.score, hit.identity) > (current.score, current.identity):
                            best_hits[query_id] = hit
            except Exception as exc:
                logger.warning("Failed to parse contig blast file %s: %s", blast_file, exc)

        contig_query_hits: Dict[str, List[BlastHit]] = defaultdict(list)
        for query_id, hit in best_hits.items():
            contig_id = protein_to_contig.get(query_id)
            if contig_id:
                contig_query_hits[contig_id].append(hit)

        total_bp = sum(contig_lengths.values()) or 1
        suspicious_bp = 0
        suspicious_contigs = []
        for contig_id, length in contig_lengths.items():
            hits = contig_query_hits.get(contig_id, [])
            strong = [hit for hit in hits if hit.identity >= 60.0 and hit.score >= 80.0]
            cellular_hits = [hit for hit in strong if hit.subject_domain in CELLULAR_PREFIXES]
            phage_hits = [hit for hit in strong if hit.subject_domain in PHAGE_PREFIXES]
            viral_hits = [hit for hit in strong if hit.subject_domain in GIANT_VIRUS_PREFIXES]
            foreign_viral = [
                hit for hit in viral_hits if primary_order and hit.subject_order and hit.subject_order != primary_order
            ]
            markers = contig_markers.get(contig_id, set())
            viral_marker_count = len(markers - CELLULAR_MODELS - PHAGE_MODELS_SET)
            cellular_marker_count = len(contig_cellular_markers.get(contig_id, set()))
            phage_marker_count = len(contig_phage_markers.get(contig_id, set()))
            nonviral_fraction = 0.0 if not strong else (len(cellular_hits) + len(phage_hits)) / len(strong)
            foreign_viral_fraction = 0.0 if not viral_hits else len(foreign_viral) / len(viral_hits)

            suspicious = False
            reason = ""
            if cellular_marker_count >= 2:
                suspicious = True
                reason = "cellular_markers"
            elif len(cellular_hits) >= 3 and nonviral_fraction >= 0.5 and viral_marker_count <= 1:
                suspicious = True
                reason = "cellular_hits"
            elif phage_marker_count >= 2 and viral_marker_count <= 1:
                suspicious = True
                reason = "phage_markers"
            elif len(phage_hits) >= 3 and nonviral_fraction >= 0.5 and viral_marker_count <= 1:
                suspicious = True
                reason = "phage_hits"
            elif len(foreign_viral) >= 2 and foreign_viral_fraction >= 0.6 and viral_marker_count >= 2:
                suspicious = True
                reason = "viral_mixture"

            if suspicious:
                suspicious_bp += length
                suspicious_contigs.append(
                    {
                        "contig_id": contig_id,
                        "length_bp": length,
                        "reason": reason,
                        "cellular_marker_count": cellular_marker_count,
                        "phage_marker_count": phage_marker_count,
                        "viral_marker_count": viral_marker_count,
                        "nonviral_fraction": round(nonviral_fraction * 100.0, 2),
                        "foreign_viral_fraction": round(foreign_viral_fraction * 100.0, 2),
                    }
                )

        suspicious_bp_fraction = (suspicious_bp / total_bp) * 100.0
        return {
            "contig_count": len(contig_lengths),
            "suspicious_contig_count": len(suspicious_contigs),
            "suspicious_bp_fraction": round(suspicious_bp_fraction, 2),
            "suspicious_contigs": suspicious_contigs,
        }


def create_contamination_scorer(database_path: Path) -> ContaminationScorer:
    return ContaminationScorer(database_path)
