"""Build the frozen experimental quality-model features from query evidence."""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from src.core.contamination_scoring import (
    CELLULAR_PREFIXES,
    GIANT_VIRUS_PREFIXES,
    PHAGE_PREFIXES,
    BlastHit,
    ContaminationScorer,
)

__all__ = [
    "COMPLETENESS_MODEL_FEATURES",
    "COMPLETENESS_NUMERIC_INPUTS",
    "COMPLETENESS_RATIO_INPUTS",
    "COMPLETENESS_RAW_INPUTS",
    "COMPOSITION_FEATURES",
    "CONTAMINATION_FEATURES",
    "CompositionWindows",
    "ContaminationFeatures",
    "build_completeness_features",
    "build_completeness_features_from_counts",
    "build_composition_windows",
    "build_contamination_features",
    "build_family_copy_counts",
]

COMPLETENESS_NUMERIC_INPUTS = (
    "runtime_LENbp",
    "runtime_genecount",
    "runtime_CODINGperc",
    "runtime_contigs",
    "runtime_avgdist",
    "runtime_gvog4_dup",
    "runtime_gvog8_dup",
    "runtime_vp_df",
    "runtime_mirus_df",
    "runtime_mrya_dup",
    "runtime_phage_dup",
    "runtime_busco_dup",
    "runtime_cog_dup",
    "runtime_ncldv_mcp_total",
    "runtime_vp_mcp",
)
COMPLETENESS_RATIO_INPUTS = (
    "runtime_gvog4_completeness",
    "runtime_gvog8_completeness",
    "runtime_vp_completeness",
    "runtime_mirus_completeness",
    "runtime_mrya_completeness",
    "runtime_phage_completeness",
    "runtime_busco_completeness",
    "runtime_cog_completeness",
)
COMPLETENESS_RAW_INPUTS = COMPLETENESS_NUMERIC_INPUTS + COMPLETENESS_RATIO_INPUTS
COMPLETENESS_MODEL_FEATURES = COMPLETENESS_NUMERIC_INPUTS + tuple(
    f"{name}_fraction" for name in COMPLETENESS_RATIO_INPUTS
)
CONTAMINATION_FEATURES = (
    "spatial_cellular_vote_fraction",
    "spatial_cellular_window_bp_fraction",
    "spatial_longest_cellular_region_bp_fraction",
    "spatial_marker_uncovered_bp_fraction",
    "order_dup",
    "gvog8_dup",
    "contamination_nonviral_hit_fraction_v1",
    "cellular_coherent_protein_fraction",
    "cellular_lineage_purity_median",
    "cellular_hit_identity_median",
    "runtime_CODINGperc",
    "spatial_viral_family_mismatch_vote_fraction",
    "spatial_viral_genus_mismatch_vote_fraction",
    "composition_cellular_score_mean",
)

REGIONAL_WINDOW_BP = 10_000
COMPOSITION_WINDOW_BP = 2_000
KMER_SIZE = 4
MINIMUM_VALID_KMER_FRACTION = 0.9
MINIMUM_COMPOSITION_COVERAGE_PCT = 80.0
IUPAC = frozenset("ACGTRYSWKMBDHVN")
IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


@dataclass(frozen=True, slots=True)
class Coordinate:
    """One protein's one-based inclusive coordinates on a query contig."""

    contig: str
    start: int
    end: int


@dataclass(frozen=True, slots=True)
class CompositionWindows:
    """Supported 2-kb composition windows and query-level support."""

    features: pd.DataFrame
    weights_bp: np.ndarray
    status: str
    total_bp: int
    scored_bp: int
    tail_bp: int
    unscored_full_window_bp: int
    scored_window_count: int
    full_window_count: int
    coverage_pct: float


@dataclass(frozen=True, slots=True)
class ContaminationFeatures:
    """Regional model inputs plus composition windows and lineage support."""

    base_features: pd.DataFrame
    composition: CompositionWindows
    lineage_status: str
    family_vote_coverage_pct: float
    genus_vote_coverage_pct: float


def build_family_copy_counts(
    marker_hits: Mapping[str, set[str]],
    model_to_family: Mapping[str, str],
    review_models: Sequence[str],
) -> dict[str, int]:
    """Count unique protein IDs within each trusted marker family."""
    review = frozenset(review_models)
    family_hits: dict[str, set[str]] = defaultdict(set)
    for model, proteins in marker_hits.items():
        if model in review:
            continue
        family = model_to_family.get(model)
        if family is None:
            raise ValueError(
                f"Marker model is absent from the frozen catalogue: {model}"
            )
        family_hits[family].update(proteins)
    return {family: len(proteins) for family, proteins in sorted(family_hits.items())}


def build_completeness_features(
    summary: Mapping[str, object],
    marker_hits: Mapping[str, set[str]],
    *,
    family_vocabulary: Sequence[str],
    model_to_family: Mapping[str, str],
    review_models: Sequence[str],
) -> pd.DataFrame:
    """Return one ordered completeness feature row from observed query evidence."""
    family_counts = build_family_copy_counts(
        marker_hits,
        model_to_family,
        review_models,
    )
    return build_completeness_features_from_counts(
        summary,
        family_counts,
        family_vocabulary=family_vocabulary,
    )


def build_completeness_features_from_counts(
    summary: Mapping[str, object],
    family_counts: Mapping[str, int],
    *,
    family_vocabulary: Sequence[str],
) -> pd.DataFrame:
    """Return one ordered row from prevalidated exact family-union counts."""
    row: dict[str, float] = {}
    for name in COMPLETENESS_NUMERIC_INPUTS:
        row[name] = _numeric_summary_value(summary, name)
    for name in COMPLETENESS_RATIO_INPUTS:
        row[f"{name}_fraction"] = _parse_ratio(_summary_value(summary, name))
    for family in family_vocabulary:
        count = family_counts.get(family, 0)
        row[f"family_present::{family}"] = float(count > 0)
    for family in family_vocabulary:
        row[f"family_count::{family}"] = float(family_counts.get(family, 0))
    return pd.DataFrame([row])


def build_composition_windows(query_fna: Path) -> CompositionWindows:
    """Build supported nonoverlapping 2-kb composition windows from nucleotide input."""
    sequences = _read_fasta(query_fna)
    rows: list[np.ndarray] = []
    total_bp = sum(len(sequence) for sequence in sequences.values())
    full_window_count = 0
    tail_bp = 0
    unscored_bp = 0
    for sequence in sequences.values():
        full_stop = len(sequence) - len(sequence) % COMPOSITION_WINDOW_BP
        tail_bp += len(sequence) - full_stop
        for start in range(0, full_stop, COMPOSITION_WINDOW_BP):
            full_window_count += 1
            features = _composition_vector(
                sequence[start : start + COMPOSITION_WINDOW_BP]
            )
            if features is None:
                unscored_bp += COMPOSITION_WINDOW_BP
                continue
            rows.append(features)

    scored_bp = len(rows) * COMPOSITION_WINDOW_BP
    coverage_pct = 100.0 * scored_bp / total_bp
    status = "available"
    if not rows:
        status = "no_supported_windows"
    elif coverage_pct < MINIMUM_COMPOSITION_COVERAGE_PCT:
        status = "coverage_below_80_pct"
    matrix = pd.DataFrame(rows, columns=COMPOSITION_FEATURES, dtype=np.float32)
    return CompositionWindows(
        features=matrix,
        weights_bp=np.full(len(rows), COMPOSITION_WINDOW_BP, dtype=float),
        status=status,
        total_bp=total_bp,
        scored_bp=scored_bp,
        tail_bp=tail_bp,
        unscored_full_window_bp=unscored_bp,
        scored_window_count=len(rows),
        full_window_count=full_window_count,
        coverage_pct=coverage_pct,
    )


def build_contamination_features(
    summary: Mapping[str, object],
    *,
    query_fna: Path,
    query_output_dir: Path,
    tree_nn_results: Mapping[str, Mapping[str, Mapping[str, float]]],
    accepted_lineage: Mapping[str, str],
    contamination_scorer: ContaminationScorer,
) -> ContaminationFeatures:
    """Return the frozen regional, mismatch, and composition feature inputs."""
    sequences = _read_fasta(query_fna)
    query_id = str(_summary_value(summary, "query"))
    gff_path = query_output_dir / "query_gff" / f"{query_id}.gff"
    coordinates = _read_gff_coordinates(gff_path)
    normalized_trees = {
        str(marker): {
            str(protein): {
                str(name): float(distance) for name, distance in rows.items()
            }
            for protein, rows in proteins.items()
        }
        for marker, proteins in tree_nn_results.items()
    }
    best_hits = contamination_scorer.collect_best_blast_hits(
        query_output_dir / "blastp_out"
    )
    regional, lineages, markers_by_window = _regional_features(
        sequences,
        coordinates,
        normalized_trees,
        best_hits,
        contamination_scorer,
    )
    mismatch, coverage = _mismatch_features(
        sequences,
        markers_by_window,
        lineages,
        normalized_trees,
        accepted_lineage,
        contamination_scorer,
    )
    values = {
        **regional,
        "order_dup": _numeric_summary_value(summary, "order_dup"),
        "gvog8_dup": _numeric_summary_value(summary, "gvog8_dup"),
        "contamination_nonviral_hit_fraction_v1": _numeric_summary_value(
            summary,
            "contamination_nonviral_hit_fraction_v1",
        ),
        "cellular_coherent_protein_fraction": _numeric_summary_value(
            summary,
            "cellular_coherent_protein_fraction",
        ),
        "cellular_lineage_purity_median": _numeric_summary_value(
            summary,
            "cellular_lineage_purity_median",
        ),
        "cellular_hit_identity_median": _numeric_summary_value(
            summary,
            "cellular_hit_identity_median",
        ),
        "runtime_CODINGperc": _numeric_summary_value(summary, "runtime_CODINGperc"),
        **mismatch,
    }
    base_names = CONTAMINATION_FEATURES[:-1]
    base_features = pd.DataFrame(
        [[values[name] for name in base_names]], columns=base_names
    )
    lineage_status = _lineage_status(accepted_lineage, coverage)
    return ContaminationFeatures(
        base_features=base_features,
        composition=build_composition_windows(query_fna),
        lineage_status=lineage_status,
        family_vote_coverage_pct=coverage["family"],
        genus_vote_coverage_pct=coverage["genus"],
    )


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(IUPAC_COMPLEMENT)[::-1]


def _canonical_4mers() -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                min(word, _reverse_complement(word))
                for bases in itertools.product("ACGT", repeat=KMER_SIZE)
                for word in ("".join(bases),)
            }
        )
    )


CANONICAL_4MERS = _canonical_4mers()
COMPOSITION_FEATURES = (
    "gc_acgt_fraction",
    *(f"k4_{word}_proportion" for word in CANONICAL_4MERS),
)


def _kmer_lookup() -> tuple[np.ndarray, np.ndarray]:
    base_codes = np.full(256, -1, dtype=np.int8)
    for code, base in enumerate("ACGT"):
        base_codes[ord(base)] = code
    family_by_word = {word: index for index, word in enumerate(CANONICAL_4MERS)}
    family_codes = np.empty(4**KMER_SIZE, dtype=np.int16)
    for bases in itertools.product("ACGT", repeat=KMER_SIZE):
        word = "".join(bases)
        code = sum(
            "ACGT".index(base) * (4 ** (3 - index)) for index, base in enumerate(word)
        )
        family_codes[code] = family_by_word[min(word, _reverse_complement(word))]
    return base_codes, family_codes


BASE_CODES, KMER_FAMILY_CODES = _kmer_lookup()


def _summary_value(summary: Mapping[str, object], name: str) -> object:
    runtime_name = name.removeprefix("runtime_")
    if runtime_name not in summary:
        raise KeyError(f"Summary lacks experimental quality feature: {runtime_name}")
    return summary[runtime_name]


def _numeric_summary_value(summary: Mapping[str, object], name: str) -> float:
    value = _summary_value(summary, name)
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _parse_ratio(value: object) -> float:
    try:
        numerator, denominator = str(value).split("/", 1)
        return float(numerator) / float(denominator)
    except (ValueError, ZeroDivisionError):
        return math.nan


def _read_fasta(path: Path) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Nucleotide query is missing: {path}")
    sequences: dict[str, str] = {}
    record_id: str | None = None
    parts: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if stripped.startswith(">"):
                if record_id is not None:
                    _store_sequence(sequences, record_id, parts, path)
                record_id = stripped[1:].split(maxsplit=1)[0]
                parts = []
                continue
            if stripped:
                if record_id is None:
                    raise ValueError(
                        f"Sequence precedes FASTA header at {path}:{line_number}"
                    )
                parts.append(stripped)
    if record_id is not None:
        _store_sequence(sequences, record_id, parts, path)
    if not sequences:
        raise ValueError(f"Nucleotide query has no FASTA records: {path}")
    return sequences


def _store_sequence(
    sequences: dict[str, str],
    record_id: str,
    parts: Sequence[str],
    path: Path,
) -> None:
    sequence = "".join(parts).upper()
    if record_id in sequences or not sequence:
        raise ValueError(f"Nucleotide query has an empty or duplicate record in {path}")
    if not set(sequence) <= IUPAC:
        raise ValueError(f"Nucleotide query contains non-IUPAC DNA in {record_id}")
    sequences[record_id] = sequence


def _composition_vector(sequence: str) -> np.ndarray | None:
    raw = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    encoded = BASE_CODES[raw]
    valid = (
        (encoded[:-3] >= 0)
        & (encoded[1:-2] >= 0)
        & (encoded[2:-1] >= 0)
        & (encoded[3:] >= 0)
    )
    valid_count = int(np.count_nonzero(valid))
    possible = len(sequence) - KMER_SIZE + 1
    if valid_count / possible < MINIMUM_VALID_KMER_FRACTION:
        return None
    acgt = encoded >= 0
    acgt_count = int(np.count_nonzero(acgt))
    codes = (
        encoded[:-3].astype(np.int16) * 64
        + encoded[1:-2].astype(np.int16) * 16
        + encoded[2:-1].astype(np.int16) * 4
        + encoded[3:].astype(np.int16)
    )
    counts = np.bincount(
        KMER_FAMILY_CODES[codes[valid]],
        minlength=len(CANONICAL_4MERS),
    )
    gc_fraction = float(np.count_nonzero((encoded == 1) | (encoded == 2))) / acgt_count
    proportions = counts.astype(np.float32) / np.float32(valid_count)
    return np.concatenate((np.array([gc_fraction], dtype=np.float32), proportions))


def _read_gff_coordinates(path: Path) -> dict[str, Coordinate]:
    if not path.is_file():
        raise FileNotFoundError(f"Query GFF is missing: {path}")
    coordinates: dict[str, Coordinate] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise ValueError(f"Malformed GFF row at {path}:{line_number}")
            attributes = dict(
                item.split("=", 1) for item in fields[8].split(";") if "=" in item
            )
            protein_id = attributes.get("ID", "")
            if not protein_id or protein_id in coordinates:
                raise ValueError(
                    f"Missing or duplicate GFF protein ID at {path}:{line_number}"
                )
            coordinates[protein_id] = Coordinate(
                contig=fields[0],
                start=int(fields[3]),
                end=int(fields[4]),
            )
    return coordinates


def _regional_features(
    sequences: Mapping[str, str],
    coordinates: Mapping[str, Coordinate],
    tree_results: dict[str, dict[str, dict[str, float]]],
    best_hits: Mapping[str, BlastHit],
    scorer: ContaminationScorer,
) -> tuple[
    dict[str, float],
    dict[str, dict[str, str] | None],
    dict[tuple[str, int], list[str]],
]:
    marker_ids = set().union(*(set(rows) for rows in tree_results.values()))
    _validate_coordinates(marker_ids, coordinates, sequences)
    markers_by_window: dict[tuple[str, int], list[str]] = defaultdict(list)
    lineages: dict[str, dict[str, str] | None] = {}
    for protein_id in sorted(marker_ids):
        coordinate = coordinates[protein_id]
        midpoint = (coordinate.start + coordinate.end) // 2
        window = (coordinate.contig, (midpoint - 1) // REGIONAL_WINDOW_BP)
        markers_by_window[window].append(protein_id)
        lineages[protein_id] = scorer._per_protein_dominant_lineage(  # noqa: SLF001  # frozen runtime rule
            protein_id,
            tree_results,
        )

    total_bp = sum(map(len, sequences.values()))
    cellular_vote_bp = 0.0
    cellular_bp = 0
    uncovered_bp = 0
    longest_region_bp = 0
    for contig, sequence in sequences.items():
        run_bp = 0
        for index, start in enumerate(range(0, len(sequence), REGIONAL_WINDOW_BP)):
            window_bp = min(REGIONAL_WINDOW_BP, len(sequence) - start)
            proteins = markers_by_window.get((contig, index), [])
            purity = scorer._classify_contig_purity(  # noqa: SLF001  # frozen runtime rule
                f"{contig}:{index}",
                proteins,
                tree_results,
                dict(best_hits),
            )
            domains = [
                lineages[protein]["domain"]
                for protein in proteins
                if lineages[protein] is not None
            ]
            counted = [
                domain
                for domain in domains
                if domain in CELLULAR_PREFIXES | GIANT_VIRUS_PREFIXES | PHAGE_PREFIXES
            ]
            cellular_votes = sum(domain in CELLULAR_PREFIXES for domain in counted)
            if counted:
                cellular_vote_bp += window_bp * cellular_votes / len(counted)
            is_cellular = purity["classification"] in {
                "cellular_coherent",
                "mixed_evidence",
            }
            cellular_bp += window_bp * int(is_cellular)
            uncovered_bp += window_bp * int(not domains)
            run_bp = run_bp + window_bp if is_cellular else 0
            longest_region_bp = max(longest_region_bp, run_bp)
    return (
        {
            "spatial_cellular_vote_fraction": 100.0 * cellular_vote_bp / total_bp,
            "spatial_cellular_window_bp_fraction": 100.0 * cellular_bp / total_bp,
            "spatial_longest_cellular_region_bp_fraction": (
                100.0 * longest_region_bp / total_bp
            ),
            "spatial_marker_uncovered_bp_fraction": 100.0 * uncovered_bp / total_bp,
        },
        lineages,
        dict(markers_by_window),
    )


def _validate_coordinates(
    marker_ids: set[str],
    coordinates: Mapping[str, Coordinate],
    sequences: Mapping[str, str],
) -> None:
    missing = marker_ids - set(coordinates)
    if missing:
        raise ValueError(f"Tree proteins lack GFF coordinates: {len(missing)}")
    invalid = [
        protein_id
        for protein_id in marker_ids
        if coordinates[protein_id].contig not in sequences
        or not (
            1
            <= coordinates[protein_id].start
            <= coordinates[protein_id].end
            <= len(sequences[coordinates[protein_id].contig])
        )
    ]
    if invalid:
        raise ValueError(
            f"Tree protein coordinates fall outside query contigs: {len(invalid)}"
        )


def _mismatch_features(
    sequences: Mapping[str, str],
    markers_by_window: Mapping[tuple[str, int], list[str]],
    lineages: dict[str, dict[str, str] | None],
    tree_results: dict[str, dict[str, dict[str, float]]],
    accepted: Mapping[str, str],
    scorer: ContaminationScorer,
) -> tuple[dict[str, float], dict[str, float]]:
    for protein_id, lineage in lineages.items():
        if lineage is not None:
            lineage["family"] = _winning_family(
                protein_id,
                lineage,
                tree_results,
                scorer,
            )
    values: dict[str, float] = {}
    coverage: dict[str, float] = {}
    for rank in ("family", "genus"):
        signal, covered = _rank_mismatch(
            rank,
            str(accepted.get(rank, "")),
            sequences,
            markers_by_window,
            lineages,
        )
        values[f"spatial_viral_{rank}_mismatch_vote_fraction"] = signal
        coverage[rank] = covered
    return values, coverage


def _winning_family(
    protein_id: str,
    lineage: Mapping[str, str],
    tree_results: Mapping[str, Mapping[str, Mapping[str, float]]],
    scorer: ContaminationScorer,
) -> str:
    families: Counter[str] = Counter()
    for marker_rows in tree_results.values():
        for neighbor in marker_rows.get(protein_id, {}):
            genome_id = scorer._parse_subject_genome_id(neighbor)  # noqa: SLF001  # frozen runtime rule
            taxonomy = scorer.labels.get(genome_id)
            if taxonomy is None:
                continue
            if all(
                taxonomy.get(rank, "") == lineage.get(rank, "")
                for rank in ("domain", "order", "genus")
            ):
                family = taxonomy.get("family", "")
                if _usable_rank(family):
                    families[family] += 1
    if not families:
        return ""
    return min(families, key=lambda family: (-families[family], family))


def _rank_mismatch(
    rank: str,
    accepted_rank: str,
    sequences: Mapping[str, str],
    markers_by_window: Mapping[tuple[str, int], list[str]],
    lineages: Mapping[str, Mapping[str, str] | None],
) -> tuple[float, float]:
    total_bp = sum(map(len, sequences.values()))
    covered_bp = 0
    mismatch_bp = 0.0
    for (contig, window_index), proteins in markers_by_window.items():
        assignments = [
            lineage[rank]
            for protein in proteins
            if (lineage := lineages[protein]) is not None
            and lineage.get("domain") == "NCLDV"
            and _usable_rank(lineage.get(rank, ""))
        ]
        if not assignments:
            continue
        start = window_index * REGIONAL_WINDOW_BP
        window_bp = min(REGIONAL_WINDOW_BP, len(sequences[contig]) - start)
        covered_bp += window_bp
        if _usable_rank(accepted_rank):
            mismatches = sum(value != accepted_rank for value in assignments)
            mismatch_bp += window_bp * mismatches / len(assignments)
    signal = 100.0 * mismatch_bp / total_bp if _usable_rank(accepted_rank) else math.nan
    return signal, 100.0 * covered_bp / total_bp


def _lineage_status(accepted: Mapping[str, str], coverage: Mapping[str, float]) -> str:
    if not _usable_rank(accepted.get("family", "")):
        return "accepted_family_missing"
    if not _usable_rank(accepted.get("genus", "")):
        return "accepted_genus_missing"
    if not math.isfinite(coverage["family"]) or coverage["family"] <= 0:
        return "family_vote_coverage_missing"
    if not math.isfinite(coverage["genus"]) or coverage["genus"] <= 0:
        return "genus_vote_coverage_missing"
    return "supported"


def _usable_rank(value: object) -> bool:
    normalized = str(value).strip().lower()
    return normalized not in {"", "nan", "none"} and not any(
        token in normalized for token in ("unclassified", "singleton", "other")
    )
