"""Marker HMM search with pyhmmer.

The search keeps only the hits and domains that pass the reporting thresholds,
parses each model file once per process, and writes the three marker tables the
rest of the pipeline reads: ``models.out`` (one row per reported domain),
``models.out.filtered`` (the best row per protein and model), and
``models.counts`` (proteins per model).
"""

from __future__ import annotations

import logging
import threading
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock
from pyhmmer.plan7 import HMM, Hit

logger = logging.getLogger(__name__)

__all__ = [
    "COMBINED_HMM",
    "DomainHit",
    "best_domain_per_pair",
    "count_proteins_per_model",
    "digitize_proteins",
    "load_profiles",
    "read_proteins",
    "run_marker_search",
    "search",
    "search_default_thresholds",
    "write_counts",
    "write_domain_table",
]

#: Location of the marker models inside a GVClass database.
COMBINED_HMM = "hmm/combined.hmm"

#: Sequence and domain E-value threshold of sensitive mode, and of models without
#: gathering cutoffs in the other mode.
EVALUE_THRESHOLD = 1e-5

RAW_TABLE = "models.out"
FILTERED_TABLE = "models.out.filtered"
COUNTS_TABLE = "models.counts"

# Downstream readers take the protein from column 1, the model from column 4,
# the sequence score from column 8 and the domain score from column 14.
_TABLE_HEADER = (
    "# target name\taccession\ttlen\tquery name\taccession\tqlen\tE-value\tscore\t"
    "bias\t#\tof\tc-Evalue\ti-Evalue\tscore\tbias\tfrom\tto\tfrom\tto\tfrom\tto\t"
    "acc\n"
)

_AMINO = Alphabet.amino()
_profile_cache: dict[Path, tuple[tuple[int, int], tuple[HMM, ...]]] = {}
_profile_cache_lock = threading.Lock()


@dataclass(frozen=True, slots=True)
class DomainHit:
    """One reported domain of one reported hit, in domain-table terms.

    Coordinates are 1-based and inclusive. ``domain_index`` counts the reported
    domains of the hit, as does ``domain_count``.
    """

    protein: str
    protein_length: int
    model: str
    model_accession: str
    model_length: int
    evalue: float
    score: float
    bias: float
    domain_index: int
    domain_count: int
    c_evalue: float
    i_evalue: float
    domain_score: float
    domain_bias: float
    hmm_from: int
    hmm_to: int
    ali_from: int
    ali_to: int
    env_from: int
    env_to: int


def load_profiles(hmm_path: Path | str) -> tuple[HMM, ...]:
    """Return the models of an HMM file, parsing the file once per process.

    The cache is keyed by the resolved path and replaced when the file size or
    modification time changes. The call is safe to make from several threads.

    Raises:
        FileNotFoundError: If ``hmm_path`` is not a file.
    """
    path = Path(hmm_path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"HMM file not found: {path}")
    stat = path.stat()
    signature = (stat.st_size, stat.st_mtime_ns)
    with _profile_cache_lock:
        cached = _profile_cache.get(path)
        if cached is not None and cached[0] == signature:
            return cached[1]
        profiles = _parse_profiles(path)
        _profile_cache[path] = (signature, profiles)
        return profiles


def read_proteins(faa_path: Path | str) -> DigitalSequenceBlock:
    """Read a protein FASTA file, forcing the amino alphabet."""
    with pyhmmer.easel.SequenceFile(
        str(faa_path), digital=True, alphabet=_AMINO
    ) as handle:
        return handle.read_block()


def digitize_proteins(proteins: Iterable[tuple[str, str]]) -> DigitalSequenceBlock:
    """Build a search target from ``(name, sequence)`` pairs."""
    sequences = [
        pyhmmer.easel.TextSequence(name=name.encode(), sequence=sequence).digitize(
            _AMINO
        )
        for name, sequence in proteins
    ]
    return DigitalSequenceBlock(_AMINO, sequences)


def search(
    profiles: Sequence[HMM],
    sequences: DigitalSequenceBlock,
    *,
    threads: int,
    sensitive: bool,
) -> list[DomainHit]:
    """Search the proteins and return the reported domains of reported hits.

    With ``sensitive`` every model uses ``E = domE = 1e-5``. Otherwise a model
    with gathering cutoffs uses them and a model without them uses the E-value
    threshold. Rows follow the search calls (models with gathering cutoffs
    first), then the model, the hit and the domain.
    """
    if sensitive:
        return _search(
            profiles, sequences, threads, E=EVALUE_THRESHOLD, domE=EVALUE_THRESHOLD
        )
    with_cutoffs = [hmm for hmm in profiles if hmm.cutoffs.gathering is not None]
    without_cutoffs = [hmm for hmm in profiles if hmm.cutoffs.gathering is None]
    return _search(with_cutoffs, sequences, threads, bit_cutoffs="gathering") + _search(
        without_cutoffs, sequences, threads, E=EVALUE_THRESHOLD, domE=EVALUE_THRESHOLD
    )


def search_default_thresholds(
    profiles: Sequence[HMM], sequences: DigitalSequenceBlock, *, threads: int = 1
) -> list[DomainHit]:
    """Search with pyhmmer's default thresholds (``E = domE = 10``)."""
    return _search(profiles, sequences, threads)


def best_domain_per_pair(domains: Iterable[DomainHit]) -> list[DomainHit]:
    """Keep the best row per protein and model.

    The best row has the highest sequence score, then the highest domain score;
    the first row wins a full tie. Rows keep the order of their first appearance.
    """
    best: dict[tuple[str, str], DomainHit] = {}
    for domain in domains:
        key = (domain.protein, domain.model)
        current = best.get(key)
        if current is None or (domain.score, domain.domain_score) > (
            current.score,
            current.domain_score,
        ):
            best[key] = domain
    return list(best.values())


def count_proteins_per_model(best_rows: Iterable[DomainHit]) -> dict[str, int]:
    """Count proteins per model from the output of :func:`best_domain_per_pair`."""
    return dict(Counter(row.model for row in best_rows))


def write_domain_table(domains: Iterable[DomainHit], path: Path) -> None:
    """Write domain rows in the 22-column layout downstream readers parse."""
    with path.open("w", encoding="utf-8") as handle:
        handle.write(_TABLE_HEADER)
        for domain in domains:
            handle.write(_format_row(domain))


def write_counts(counts: dict[str, int], path: Path) -> None:
    """Write ``model<TAB>count`` lines sorted by model."""
    path.write_text(
        "".join(f"{model}\t{count}\n" for model, count in sorted(counts.items())),
        encoding="utf-8",
    )


def run_marker_search(
    protein_faa: Path,
    hmm_path: Path,
    out_dir: Path,
    *,
    threads: int,
    sensitive: bool,
) -> Path:
    """Search the proteins and write the marker tables into ``out_dir``.

    Returns:
        The path of ``models.out.filtered``.
    """
    profiles = load_profiles(hmm_path)
    sequences = read_proteins(protein_faa)
    logger.info(
        "Searching %d proteins against %d models (sensitive=%s)",
        len(sequences),
        len(profiles),
        sensitive,
    )
    domains = search(profiles, sequences, threads=threads, sensitive=sensitive)
    best = best_domain_per_pair(domains)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_domain_table(domains, out_dir / RAW_TABLE)
    write_domain_table(best, out_dir / FILTERED_TABLE)
    write_counts(count_proteins_per_model(best), out_dir / COUNTS_TABLE)
    logger.info(
        "Marker search kept %d domains of %d protein-model pairs",
        len(domains),
        len(best),
    )
    return out_dir / FILTERED_TABLE


def _parse_profiles(path: Path) -> tuple[HMM, ...]:
    logger.info("Parsing HMM file %s", path)
    with pyhmmer.plan7.HMMFile(str(path)) as handle:
        return tuple(handle)


def _search(
    profiles: Sequence[HMM],
    sequences: DigitalSequenceBlock,
    threads: int,
    **thresholds: object,
) -> list[DomainHit]:
    if not profiles:
        return []
    results = pyhmmer.hmmsearch(profiles, sequences, cpus=threads, **thresholds)
    domains: list[DomainHit] = []
    for hmm, top_hits in zip(profiles, results, strict=True):
        for hit in top_hits:
            if hit.reported:
                domains.extend(_reported_domains(hmm, hit))
    return domains


def _reported_domains(hmm: HMM, hit: Hit) -> list[DomainHit]:
    reported = [domain for domain in hit.domains if domain.reported]
    return [
        DomainHit(
            protein=hit.name,
            protein_length=hit.length,
            model=hmm.name,
            model_accession=hmm.accession or "-",
            model_length=hmm.M,
            evalue=hit.evalue,
            score=hit.score,
            bias=hit.bias,
            domain_index=index,
            domain_count=len(reported),
            c_evalue=domain.c_evalue,
            i_evalue=domain.i_evalue,
            domain_score=domain.score,
            domain_bias=domain.bias,
            hmm_from=domain.alignment.hmm_from,
            hmm_to=domain.alignment.hmm_to,
            ali_from=domain.alignment.target_from,
            ali_to=domain.alignment.target_to,
            env_from=domain.env_from,
            env_to=domain.env_to,
        )
        for index, domain in enumerate(reported, start=1)
    ]


def _format_row(domain: DomainHit) -> str:
    fields = (
        domain.protein,
        "-",
        str(domain.protein_length),
        domain.model,
        domain.model_accession,
        str(domain.model_length),
        f"{domain.evalue:.2e}",
        f"{domain.score:.1f}",
        f"{domain.bias:.1f}",
        str(domain.domain_index),
        str(domain.domain_count),
        f"{domain.c_evalue:.2e}",
        f"{domain.i_evalue:.2e}",
        f"{domain.domain_score:.1f}",
        f"{domain.domain_bias:.1f}",
        str(domain.hmm_from),
        str(domain.hmm_to),
        str(domain.ali_from),
        str(domain.ali_to),
        str(domain.env_from),
        str(domain.env_to),
        "-",
    )
    return "\t".join(fields) + "\n"
