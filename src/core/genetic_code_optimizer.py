"""Genetic-code selection for nucleotide queries.

Every code in the panel gene-calls the query with pyrodigal, and the proteins
of each call are searched against the marker models. Pyrodigal's metagenomic
mode (code 0, "meta") is the baseline: another code replaces it only with
marker evidence, never with coding density alone.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import pyrodigal  # Using tomasbruna's fork with codes 106, 129
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pyhmmer.plan7 import HMM

from src.core.hmm_search import (
    DomainHit,
    digitize_proteins,
    load_profiles,
    search_default_thresholds,
)

logger = logging.getLogger(__name__)

__all__ = [
    "META_CODE",
    "AssemblyStats",
    "CodeEvaluation",
    "GeneCallingResult",
    "optimize_genetic_code",
    "select_best_code",
]

META_CODE = 0

#: Pyrodigal cannot train a single-genome model on less sequence.
MIN_TRAINING_BP = 20_000

#: A protein is a complete hit when its best domain covers more than this
#: fraction of the model.
COMPLETE_COVERAGE = 0.66

#: A code qualifies over meta with this many more complete hits, or ...
COMPLETE_HITS_MARGIN = 2
#: ... with at least as many complete hits and an average best-hit score this
#: many percent above meta.
SCORE_MARGIN_PERCENT = 10

# Preferred first. The standard codes come before the exotic ones.
_PREFERENCE = (11, 1, 4, 6, 15, 29, 106, 129)

_STATS_HEADER = f"""\
# Genetic code optimization results
# * marks the selected code
# Code 0 (meta) is kept unless another code has marker evidence over it:
#   at least 1 complete_hits, and
#   at least {COMPLETE_HITS_MARGIN} more complete_hits than meta, OR
#   at least as many complete_hits and an avg_best_hit_score at least \
{SCORE_MARGIN_PERCENT}% above meta.
# Qualifying codes rank by complete_hits, avg_best_hit_score, then code preference.
# coding_density is reported but does not decide.
# Input under {MIN_TRAINING_BP} bp is evaluated with meta only.
#
"""

_STATS_COLUMNS = (
    "query",
    "genetic_code",
    "ttable",
    "complete_hits",
    "avg_best_hit_score",
    "coding_density",
    "gene_count",
    "contigs",
    "LENbp",
    "GCperc",
)


@dataclass(frozen=True, slots=True)
class AssemblyStats:
    """Size and composition of the query, computed once."""

    contigs: int
    length_bp: int
    gc_percent: float


@dataclass(frozen=True, slots=True)
class CodeEvaluation:
    """Gene calling and marker support of one genetic code."""

    code: int
    gene_count: int
    coding_density: float
    complete_hits: int
    avg_best_hit_score: float


@dataclass(frozen=True, slots=True)
class GeneCallingResult:
    """The selected code, every evaluation, and the written gene files."""

    best_code: int
    evaluations: dict[int, CodeEvaluation]
    assembly: AssemblyStats
    faa: Path
    fna: Path
    gff: Path


def optimize_genetic_code(
    query_file: Path,
    output_dir: Path,
    hmm_path: Path,
    genetic_codes: Sequence[int],
    *,
    threads: int,
) -> GeneCallingResult:
    """Select the genetic code of a nucleotide query and write its gene files.

    Evaluations run in parallel, one single-threaded evaluation per worker. The
    protein, gene and GFF files of the selected code go to
    ``<output_dir>/gene_calling/<query stem>/``, together with a ``.statsout``
    table of all evaluations.

    Raises:
        RuntimeError: If no evaluated code predicts a gene.
    """
    records = list(SeqIO.parse(str(query_file), "fasta"))
    assembly = _assembly_stats(records)
    codes = list(genetic_codes)
    if assembly.length_bp < MIN_TRAINING_BP:
        logger.info(
            "Input has %d bp, under %d: evaluating meta mode only",
            assembly.length_bp,
            MIN_TRAINING_BP,
        )
        codes = [META_CODE]
    profiles = load_profiles(hmm_path)
    logger.info("Evaluating %d genetic codes with %d threads", len(codes), threads)
    with ThreadPoolExecutor(max_workers=max(1, min(len(codes), threads))) as pool:
        evaluations = {
            evaluation.code: evaluation
            for evaluation in pool.map(
                lambda code: _evaluate_code(records, code, assembly, profiles), codes
            )
        }
    if all(evaluation.gene_count == 0 for evaluation in evaluations.values()):
        raise RuntimeError(f"No genes predicted for {query_file} under any code")
    best_code = select_best_code(evaluations)
    best = evaluations[best_code]
    logger.info(
        "Selected genetic code %d: complete_hits=%d avg_best_hit_score=%.2f "
        "coding_density=%.2f%%",
        best_code,
        best.complete_hits,
        best.avg_best_hit_score,
        best.coding_density,
    )
    gene_dir = output_dir / "gene_calling" / query_file.stem
    gene_dir.mkdir(parents=True, exist_ok=True)
    result = GeneCallingResult(
        best_code=best_code,
        evaluations=evaluations,
        assembly=assembly,
        faa=gene_dir / f"{query_file.stem}.faa",
        fna=gene_dir / f"{query_file.stem}.fna",
        gff=gene_dir / f"{query_file.stem}.gff",
    )
    _write_gene_files(records, result)
    _write_stats(result, query_file.stem, gene_dir / f"{query_file.stem}.statsout")
    return result


def select_best_code(evaluations: Mapping[int, CodeEvaluation]) -> int:
    """Return meta unless another code has marker evidence over it.

    A nonzero code qualifies with at least one complete hit and either
    ``COMPLETE_HITS_MARGIN`` more complete hits than meta, or at least as many
    and an average best-hit score ``SCORE_MARGIN_PERCENT`` percent above meta.
    Qualifying codes rank by complete hits, then score, then ``_PREFERENCE``.
    """
    meta = evaluations[META_CODE]
    qualifying = [
        evaluation
        for code, evaluation in evaluations.items()
        if code != META_CODE and _has_evidence_over_meta(evaluation, meta)
    ]
    if not qualifying:
        return META_CODE
    return max(qualifying, key=_ranking_key).code


def _has_evidence_over_meta(candidate: CodeEvaluation, meta: CodeEvaluation) -> bool:
    if candidate.complete_hits == 0:
        return False
    if candidate.complete_hits >= meta.complete_hits + COMPLETE_HITS_MARGIN:
        return True
    return (
        candidate.complete_hits >= meta.complete_hits
        and candidate.avg_best_hit_score * 100
        >= meta.avg_best_hit_score * (100 + SCORE_MARGIN_PERCENT)
    )


def _ranking_key(evaluation: CodeEvaluation) -> tuple[int, float, int]:
    preference = (
        _PREFERENCE.index(evaluation.code)
        if evaluation.code in _PREFERENCE
        else len(_PREFERENCE)
    )
    return (
        evaluation.complete_hits,
        round(evaluation.avg_best_hit_score, 2),
        -preference,
    )


def _assembly_stats(records: Sequence[SeqRecord]) -> AssemblyStats:
    length_bp = sum(len(record.seq) for record in records)
    gc_bases = sum(
        str(record.seq).upper().count("G") + str(record.seq).upper().count("C")
        for record in records
    )
    return AssemblyStats(
        contigs=len(records),
        length_bp=length_bp,
        gc_percent=round(gc_bases / length_bp * 100, 2),
    )


def _call_genes(
    records: Sequence[SeqRecord], code: int
) -> list[tuple[str, pyrodigal.Genes]]:
    """Gene-call every record under one code; meta mode needs no training."""
    if code == META_CODE:
        finder = pyrodigal.GeneFinder(meta=True)
    else:
        finder = pyrodigal.GeneFinder(meta=False)
        finder.train(
            *(str(record.seq) for record in records),
            force_nonsd=True,
            translation_table=code,
        )
    return [(record.id, finder.find_genes(str(record.seq))) for record in records]


def _evaluate_code(
    records: Sequence[SeqRecord],
    code: int,
    assembly: AssemblyStats,
    profiles: Sequence[HMM],
) -> CodeEvaluation:
    proteins = [
        (f"{contig_id}_{number}", gene.translate())
        for contig_id, genes in _call_genes(records, code)
        for number, gene in enumerate(genes, start=1)
    ]
    domains = search_default_thresholds(profiles, digitize_proteins(proteins))
    complete_hits, avg_best_hit_score = _marker_support(domains)
    coding_bases = 3 * sum(len(sequence) for _, sequence in proteins)
    evaluation = CodeEvaluation(
        code=code,
        gene_count=len(proteins),
        coding_density=round(coding_bases / assembly.length_bp * 100, 2),
        complete_hits=complete_hits,
        avg_best_hit_score=avg_best_hit_score,
    )
    logger.info(
        "Code %d: genes=%d complete_hits=%d avg_best_hit_score=%.2f "
        "coding_density=%.2f%%",
        code,
        evaluation.gene_count,
        evaluation.complete_hits,
        evaluation.avg_best_hit_score,
        evaluation.coding_density,
    )
    return evaluation


def _marker_support(domains: Iterable[DomainHit]) -> tuple[int, float]:
    """Return the complete hits and the mean best-domain score per hit protein."""
    best_domain: dict[str, DomainHit] = {}
    for domain in domains:
        current = best_domain.get(domain.protein)
        if current is None or domain.domain_score > current.domain_score:
            best_domain[domain.protein] = domain
    if not best_domain:
        return 0, 0.0
    complete_hits = sum(
        1
        for domain in best_domain.values()
        if (domain.hmm_to - domain.hmm_from + 1) / domain.model_length
        > COMPLETE_COVERAGE
    )
    mean_score = sum(domain.domain_score for domain in best_domain.values()) / len(
        best_domain
    )
    return complete_hits, round(mean_score, 2)


def _write_gene_files(records: Sequence[SeqRecord], result: GeneCallingResult) -> None:
    with (
        result.faa.open("w", encoding="utf-8") as faa,
        result.fna.open("w", encoding="utf-8") as fna,
        result.gff.open("w", encoding="utf-8") as gff,
    ):
        for contig_id, genes in _call_genes(records, result.best_code):
            genes.write_translations(faa, contig_id)
            genes.write_genes(fna, contig_id)
            genes.write_gff(gff, contig_id, include_translation_table=True)


def _write_stats(result: GeneCallingResult, query_name: str, path: Path) -> None:
    ranked = sorted(
        result.evaluations.values(),
        key=lambda item: (
            -item.complete_hits,
            -item.avg_best_hit_score,
            -item.coding_density,
        ),
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(_STATS_HEADER)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(_STATS_COLUMNS)
        for item in ranked:
            selected = "*" if item.code == result.best_code else ""
            writer.writerow(
                (
                    query_name,
                    f"{item.code}{selected}",
                    "codemeta" if item.code == META_CODE else f"code{item.code}",
                    item.complete_hits,
                    item.avg_best_hit_score,
                    item.coding_density,
                    item.gene_count,
                    result.assembly.contigs,
                    result.assembly.length_bp,
                    result.assembly.gc_percent,
                )
            )
