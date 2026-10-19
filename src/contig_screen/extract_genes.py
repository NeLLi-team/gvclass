"""Cut one bin into windows, call genes, search the HMM catalog and the references.

For a bin, the input contigs are cut into windows of ``tile_kb`` (every
``step_kb``) and each window is a contig of its own. Genes come from pyrodigal in
metagenomic mode, so a window does not depend on its neighbours. Unique proteins
are searched once with pyhmmer against the GVClass catalog and, per marker, with
pyswrd against the marker references (the two searches GVClass runs). Three
tables are written: windows, genes and hits.
"""

from __future__ import annotations

import functools
import gzip
import hashlib
import json
import logging
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import TypeAlias

import pandas as pd
import pyhmmer
import pyrodigal
import pyswrd
from Bio import SeqIO

from src.contig_screen import layout

logger = logging.getLogger(__name__)

MIN_SOURCE_CONTIG = 3000
MAX_N_FRACTION = 0.05
MIN_PROTEIN_AA = 30
E_VALUE = 1e-5
TOP_HITS = 100
SEED = 20261015
CODONS = tuple(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT")
HIT_COLUMNS = [
    "prot_hash",
    "marker",
    "rank",
    "subject",
    "identity",
    "score",
    "evalue",
    "aln_len",
]

HmmHits: TypeAlias = dict[str, list[tuple[str, float, float]]]


def extract_bin(
    genome: dict[str, str],
    resources: layout.ScreenResources,
    out_dir: Path,
    *,
    tile_kb: int,
    step_kb: int | None,
    threads: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Write the window, gene and hit tables of a bin; return its windows and genes.

    Raises:
        ValueError: The bin has no contig of at least 3 kb.
    """
    started = time.monotonic()
    contigs = tile_contigs(genome, tile_kb, step_kb)
    if not contigs:
        msg = f"{genome['genome_id']}: no contig of at least 3 kb"
        raise ValueError(msg)
    genes, contig_rows = call_genes(genome, contigs)
    proteins = unique_proteins(genes)
    models = hmm_search(proteins, resources.path(layout.HMM_FILE), threads)
    hits = reference_search(proteins, models, resources, threads)
    write_tables(out_dir, genome["genome_id"], contig_rows, genes, models, hits)
    logger.info(
        "%s: %d windows, %d genes, %d unique proteins, %d with a model, "
        "%d hit rows, %.0f s",
        genome["genome_id"],
        len(contig_rows),
        len(genes),
        len(proteins),
        len(models),
        len(hits),
        time.monotonic() - started,
    )
    return contigs, genes


def read_fasta(path: str) -> dict[str, str]:
    """Return contig name to upper-case sequence for contigs of at least 3 kb."""
    opener = gzip.open if path.endswith(".gz") else open
    contigs: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq) >= MIN_SOURCE_CONTIG:
                contigs[record.id] = str(record.seq).upper()
    return contigs


def tile_contigs(
    genome: dict[str, str],
    tile_kb: int,
    step_kb: int | None = None,
    offset_kb: float = 0.0,
) -> list[dict[str, object]]:
    """Keep the input contigs (tile_kb 0) or cut them into windows of tile_kb.

    Windows start every ``step_kb`` (default: the window length, no overlap), the
    first at ``offset_kb`` (to test junctions that do not sit on the window grid).
    """
    records: list[dict[str, object]] = []
    for name, sequence in read_fasta(genome["fasta_path"]).items():
        size = tile_kb * 1000 if tile_kb else len(sequence)
        step = (step_kb * 1000) if (step_kb and tile_kb) else size
        first = int(offset_kb * 1000) if tile_kb else 0
        for start in range(first, len(sequence), step):
            piece = sequence[start : start + size]
            if len(piece) < MIN_SOURCE_CONTIG:
                continue
            if piece.count("N") / len(piece) > MAX_N_FRACTION:
                continue
            records.append(
                {
                    "contig_id": f"{genome['genome_id']}|{name}|{start}",
                    "length_kb": len(piece) // 1000,
                    "source_contig": name,
                    "source_start": start,
                    "sequence": piece,
                }
            )
    return records


def call_genes(
    genome: dict[str, str], contigs: list[dict[str, object]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Call genes on every contig and summarize each contig."""
    finder = pyrodigal.GeneFinder(meta=True)
    gene_rows: list[dict[str, object]] = []
    contig_rows = []
    for contig in contigs:
        sequence = str(contig["sequence"])
        genes = finder.find_genes(sequence.encode("ascii"))
        covered = 0
        last_end = 0
        codons: Counter[str] = Counter()
        for index, gene in enumerate(genes):
            gene_rows.append(gene_row(str(contig["contig_id"]), index, gene))
            begin = max(gene.begin, last_end + 1)
            covered += max(0, gene.end - begin + 1)
            last_end = max(last_end, gene.end)
            cds = gene.sequence().upper()
            codons.update(cds[i : i + 3] for i in range(0, len(cds) - 2, 3))
        row: dict[str, object] = {
            "contig_id": contig["contig_id"],
            "genome_id": genome["genome_id"],
            "class": genome["class"],
            "length_kb": contig["length_kb"],
            "length_bp": len(sequence),
            "source_contig": contig["source_contig"],
            "source_start": contig["source_start"],
            "gc": (sequence.count("G") + sequence.count("C")) / len(sequence),
            "n_genes": len(genes),
            "coding_bp": covered,
        }
        row.update({f"cod_{codon}": codons.get(codon, 0) for codon in CODONS})
        contig_rows.append(row)
    return gene_rows, contig_rows


def gene_row(contig_id: str, index: int, gene: pyrodigal.Gene) -> dict[str, object]:
    """Return the features of one called gene."""
    protein = gene.translate().rstrip("*")
    return {
        "contig_id": contig_id,
        "gene_idx": index,
        "begin": gene.begin,
        "end": gene.end,
        "strand": gene.strand,
        "partial": int(gene.partial_begin) + int(gene.partial_end),
        "start_type": gene.start_type,
        "rbs_motif": gene.rbs_motif or "None",
        "rbs_spacer": gene.rbs_spacer or "None",
        "gc_cont": gene.gc_cont,
        "score": gene.score,
        "cscore": gene.cscore,
        "rscore": gene.rscore,
        "sscore": gene.sscore,
        "tscore": gene.tscore,
        "uscore": gene.uscore,
        "conf": gene.confidence(),
        "ttable": gene.translation_table,
        "aa_len": len(protein),
        "prot_hash": hashlib.sha1(protein.encode(), usedforsecurity=False).hexdigest(),
        "protein": protein,
    }


def unique_proteins(genes: list[dict[str, object]]) -> dict[str, str]:
    """Return protein hash to sequence for the unique proteins."""
    proteins: dict[str, str] = {}
    for gene in genes:
        if int(str(gene["aa_len"])) >= MIN_PROTEIN_AA:
            proteins.setdefault(str(gene["prot_hash"]), str(gene["protein"]))
    return proteins


def text(name: bytes | str) -> str:
    """Return a pyhmmer name as text."""
    return name.decode() if isinstance(name, bytes) else str(name)


@functools.cache
def load_models(hmm_path: Path) -> list[pyhmmer.plan7.HMM]:
    """Load the catalog once per process."""
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as handle:
        return list(handle)


def hmm_search(proteins: dict[str, str], hmm_path: Path, threads: int) -> HmmHits:
    """Search the catalog; return protein hash to (model, score, evalue) hits."""
    if not proteins:
        return {}
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = [
        pyhmmer.easel.TextSequence(name=key.encode(), sequence=seq).digitize(alphabet)
        for key, seq in proteins.items()
    ]
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet, sequences)
    models = load_models(hmm_path)
    found: HmmHits = defaultdict(list)
    for top_hits in pyhmmer.hmmsearch(
        models, block, cpus=threads, E=E_VALUE, domE=E_VALUE
    ):
        model = text(top_hits.query.name)
        for hit in top_hits:
            found[text(hit.name)].append((model, float(hit.score), float(hit.evalue)))
    return dict(found)


def reference_search(
    proteins: dict[str, str],
    models: HmmHits,
    resources: layout.ScreenResources,
    threads: int,
) -> list[dict[str, object]]:
    """Search every marker's references with pyswrd, as GVClass does."""
    raw = resources.read_bytes(layout.MODEL_TO_GROUP_FILE)
    to_group = json.loads(raw.decode("utf-8"))
    by_marker: dict[str, list[str]] = defaultdict(list)
    for key, found in models.items():
        for marker in {to_group.get(model, model) for model, _, _ in found}:
            by_marker[marker].append(key)
    rows: list[dict[str, object]] = []
    for marker, keys in sorted(by_marker.items()):
        relative = f"{layout.MARKER_REFS_DIR}/{marker}.faa"
        if not resources.has(relative):
            continue
        ids, targets = read_reference(resources.path(relative))
        rows.extend(search_marker(marker, keys, proteins, ids, targets, threads))
    return rows


def search_marker(
    marker: str,
    keys: list[str],
    proteins: dict[str, str],
    ids: list[str],
    targets: list[str],
    threads: int,
) -> list[dict[str, object]]:
    """Return the top subjects of every query protein of one marker."""
    per_query: dict[int, list[tuple[str, float, float, float, int]]] = defaultdict(list)
    queries = [proteins[key] for key in keys]
    for hit in pyswrd.search(
        queries,
        targets,
        scorer_name="BLOSUM62",
        score_threshold=0,
        max_alignments=TOP_HITS,
        threads=threads,
    ):
        if hit.evalue > E_VALUE:
            continue
        result = hit.result
        per_query[hit.query_index].append(
            (
                ids[hit.target_index],
                result.identity() * 100.0,
                float(hit.score),
                float(hit.evalue),
                result.query_end - result.query_start + 1,
            )
        )
    rows: list[dict[str, object]] = []
    for query_index, found in per_query.items():
        found.sort(key=lambda item: (-item[2], item[0]))
        for rank, (subject, identity, score, evalue, aln) in enumerate(
            found[:TOP_HITS], 1
        ):
            rows.append(
                {
                    "prot_hash": keys[query_index],
                    "marker": marker,
                    "rank": rank,
                    "subject": subject,
                    "identity": identity,
                    "score": score,
                    "evalue": evalue,
                    "aln_len": aln,
                }
            )
    return rows


@functools.lru_cache(maxsize=96)
def read_reference(path: Path) -> tuple[list[str], list[str]]:
    """Read a marker reference FASTA into ids and sequences."""
    ids: list[str] = []
    seqs: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id)
            seqs.append(str(record.seq).replace("*", ""))
    return ids, seqs


def write_tables(
    out: Path,
    name: str,
    contig_rows: list[dict[str, object]],
    genes: list[dict[str, object]],
    models: HmmHits,
    hits: list[dict[str, object]],
) -> None:
    """Write the contig, gene and hit tables of one genome."""
    pd.DataFrame(contig_rows).to_parquet(out / f"{name}.contigs.parquet", index=False)
    gene_frame = pd.DataFrame(genes).drop(columns=["protein"])
    best = {key: max(found, key=lambda item: item[1]) for key, found in models.items()}
    counts = {key: len({hit[0] for hit in found}) for key, found in models.items()}
    gene_frame["best_model"] = gene_frame["prot_hash"].map(
        lambda key: best[key][0] if key in best else None
    )
    gene_frame["best_hmm_score"] = gene_frame["prot_hash"].map(
        lambda key: best[key][1] if key in best else float("nan")
    )
    gene_frame["n_models"] = gene_frame["prot_hash"].map(lambda key: counts.get(key, 0))
    gene_frame.to_parquet(out / f"{name}.genes.parquet", index=False)
    hit_frame = pd.DataFrame(hits, columns=HIT_COLUMNS)
    hit_frame.to_parquet(out / f"{name}.hits.parquet", index=False)
