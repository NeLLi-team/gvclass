"""Find bacterial, archaeal and eukaryotic sequence in a giant-virus bin.

Usage: python -m src.contig_screen.cellular <bin.fna> --resources DIR --out out.tsv
                                            [--tile-kb 10 --step-kb 5] [--threads N]

The bin is cut into windows. Each window gets genes (pyrodigal), catalog and marker
searches, and features (``extract_genes``, ``build_features``); the best hits of
its proteins to the whole-proteome reference (DIAMOND against ``ref7`` and the
catalog giant viruses ``ref7x``, ``ref_hits``); and the mean ESM-2 embedding of up
to eight proteins (``extra_features``). Stage 1 is the stack of the boosting model
and of the set model, which reads the per-gene hit tokens and protein embeddings
of every window (``hit_tokens``, ``gene_embeddings``, ``set_model``): the mean of
the logits of the two. Stage 2 rescores every window in the context of the bin
(``context_stage``). The output has one row per window with P(cellular) of stage 1
(``p_cellular``), its two parts (``p_gbm``, ``set_logit``), the stage 2 score
(``p_context``) and a call at the threshold of stage 2. The share printed at the
end is the share of base pairs called. A call is a review flag, not a
contamination estimate.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd

from src.contig_screen import layout, torch_stack
from src.contig_screen.build_features import load_labels, unit_features
from src.contig_screen.context_stage import context_table, tilings
from src.contig_screen.extra_features import Embedder, draw_plan, window_embeddings
from src.contig_screen.extract_genes import extract_bin, read_fasta
from src.contig_screen.gene_embeddings import extend_tokens
from src.contig_screen.hit_tokens import SCALARS, window_tokens
from src.contig_screen.ref_hits import load_reference, reference_features
from src.contig_screen.stack_context import stack_scores
from src.contig_screen.thresholds import sliding_calls

logger = logging.getLogger(__name__)

ESM_GENES = 8
SET_VARIANT = "hits_struct"
MOSTLY_CELLULAR_PCT = 80.0
SKIPPED_CLASSES = ["NCLDV"]
DIAMOND_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue"]
OUTPUT_COLUMNS = [
    "contig_id",
    "source_contig",
    "source_start",
    "length_bp",
    "s_n_genes",
    "c_gc",
    "s_coding_density",
    "p_cellular",
    "p_gbm",
    "set_logit",
    "p_context",
]


@dataclass(frozen=True, slots=True)
class Settings:
    """How a bin is cut and how many threads the searches and models use."""

    tile_kb: int = 10
    step_kb: int = 5
    threads: int = 4
    diamond: Path | None = None
    sample_seed: str | None = None


def build_parser() -> argparse.ArgumentParser:
    """Return the command line parser."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("fasta", type=Path, help="bin (FASTA, optionally gzipped)")
    parser.add_argument(
        "--resources",
        type=Path,
        required=True,
        help="resource bundle root (holds contig_screen/, hmm/, parquet/)",
    )
    parser.add_argument("--out", type=Path, required=True, help="output table (TSV)")
    parser.add_argument("--tile-kb", type=int, default=10, help="window length (kb)")
    parser.add_argument(
        "--step-kb",
        type=int,
        default=5,
        help="window step in kb (0: no overlap)",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="stage 1 threshold (default: thresholds.json)",
    )
    parser.add_argument(
        "--threshold2",
        type=float,
        default=None,
        help="stage 2 threshold (default: thresholds.json)",
    )
    parser.add_argument(
        "--keep", type=Path, default=None, help="keep the intermediate files here"
    )
    parser.add_argument(
        "--diamond", type=Path, default=None, help="DIAMOND binary (default: on PATH)"
    )
    parser.add_argument(
        "--sample-seed",
        default=None,
        metavar="NAME",
        help=(
            "reproduce numbers of the research prototype; NAME is the bin file "
            "stem (default: a seed from the sequences, independent of file name "
            "and contig order)"
        ),
    )
    parser.add_argument(
        "--manifest-sha256",
        default=None,
        help="digest of manifest.json, for a bundle other than the release bundle",
    )
    parser.add_argument(
        "--verify-all",
        action="store_true",
        help="hash every file of the bundle, DIAMOND databases included",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Score one bin, write the table and print the summary."""
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    try:
        torch_stack.require_torch()
        resources = layout.ScreenResources.open(
            args.resources,
            manifest_sha256=args.manifest_sha256,
            verify_all=args.verify_all,
        )
    except (torch_stack.MissingTorchError, layout.IntegrityError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    cuts = layout.load_thresholds(resources)
    stage1 = cuts["stage1"] if args.threshold is None else args.threshold
    stage2 = cuts["stage2"] if args.threshold2 is None else args.threshold2
    settings = Settings(
        args.tile_kb, args.step_kb, args.threads, args.diamond, args.sample_seed
    )
    with torch_stack.limit_threads(args.threads):
        result = score_bin(args.fasta, resources, settings, keep=args.keep)
    result["called_cellular"] = result["p_context"] >= stage2
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    total = fasta_bp(args.fasta)
    lengths = contig_lengths(args.fasta) if settings.step_kb else None
    summarize(result, stage2, total, lengths, "p_context")
    print("stage 1 alone (no bin context):")
    alone = result.assign(called_cellular=result["p_cellular"] >= stage1)
    summarize(alone, stage1, total, lengths, "p_cellular")
    return 0


def score_bin(
    fasta: Path,
    resources: layout.ScreenResources,
    settings: Settings,
    *,
    keep: Path | None = None,
) -> pd.DataFrame:
    """Return the windows of a bin with their stage 1 and stage 2 scores.

    The choice of the proteins for the embedding is seeded from the sequences of
    the bin, so the file name and the contig order do not change the scores (see
    ``extra_features.draw_plan``). The intermediate files go to ``keep`` or a
    temporary directory.
    """
    boosting = layout.load_joblib(resources, layout.BOOSTING_FILE)
    with tempfile.TemporaryDirectory() as scratch:
        work = keep if keep is not None else Path(scratch)
        work.mkdir(parents=True, exist_ok=True)
        features = window_features(fasta, resources, settings, boosting, work)
    for column in boosting["columns"]:
        if column not in features:
            features[column] = float("nan")
    features = stage1_scores(features, boosting)
    features["p_context"] = context_scores(features, resources, settings)
    return features[OUTPUT_COLUMNS].copy()


def window_features(
    fasta: Path,
    resources: layout.ScreenResources,
    settings: Settings,
    boosting: dict[str, object],
    work: Path,
) -> pd.DataFrame:
    """Return the feature table of the windows, with the set model logit joined."""
    genome = {
        "genome_id": fasta.stem,
        "class": "BIN",
        "fasta_path": str(fasta.resolve()),
    }
    extract = work / "extract"
    extract.mkdir(parents=True, exist_ok=True)
    contigs, genes = extract_bin(
        genome,
        resources,
        extract,
        tile_kb=settings.tile_kb,
        step_kb=settings.step_kb or None,
        threads=settings.threads,
    )
    labels = load_labels(
        resources.path(layout.LABELS_FILE), resources.path(layout.LABEL_DOMAIN_FILE)
    )
    features = unit_features(extract, genome["genome_id"], genome["class"], labels)
    queries = [(f"{g['contig_id']}|{g['gene_idx']}", str(g["protein"])) for g in genes]
    index = query_index(queries)
    features, best = add_reference_hits(
        features, queries, index, resources, settings, work
    )
    if any(c.startswith("e_") for c in boosting["columns"]):
        pooling = str(boosting.get("embedding_pooling", "legacy"))
        embedder = make_embedder(resources, settings, pooling)
        seed, drawn = draw_plan(
            contigs,
            read_fasta(genome["fasta_path"]),
            sample_seed=settings.sample_seed,
        )
        embedded = window_embeddings(seed, drawn, genes, embedder, ESM_GENES)
        keep = ["contig_id", *[c for c in embedded.columns if c.startswith("e_")]]
        features = features.merge(embedded[keep], on="contig_id", how="left")
    return add_set_logits(features, best, queries, index, resources, settings)


def make_embedder(
    resources: layout.ScreenResources, settings: Settings, pooling: str
) -> Embedder:
    """Return an ESM-2 embedder over the verified local model files."""
    return Embedder(resources.directory(layout.ESM_DIR), settings.threads, pooling)


def stage1_scores(features: pd.DataFrame, boosting: dict[str, object]) -> pd.DataFrame:
    """Return the windows with ``p_gbm`` (boosting) and ``p_cellular`` (the stack)."""
    x = features[boosting["columns"]].to_numpy(dtype=float)
    gbm = boosting["model"].predict_proba(x)[:, 1]
    missing = int(features["set_logit"].isna().sum())
    if missing:
        msg = f"{missing} windows have no set logit"
        raise ValueError(msg)
    stack = stack_scores(gbm, features["set_logit"].to_numpy())
    return features.assign(p_gbm=gbm, p_cellular=stack)


def context_scores(
    features: pd.DataFrame, resources: layout.ScreenResources, settings: Settings
) -> np.ndarray:
    """Return stage 2 scores of the windows of the bin (one bin)."""
    stage2 = layout.load_joblib(resources, layout.STAGE2_FILE)
    norm = layout.load_npz(resources, layout.STAGE2_NORM_FILE)
    frame = features.rename(columns={"p_cellular": "p"}).copy()
    phase = tilings(frame, settings.tile_kb, settings.step_kb)
    frame["bin_id"] = "bin#" + phase.astype(str)
    table = context_table(frame, list(norm["comp"]), norm["mu"], norm["sd"])
    x = table[stage2["columns"]].to_numpy(dtype=float)
    return stage2["model"].predict_proba(x)[:, 1]


def query_index(queries: Sequence[tuple[str, str]]) -> pd.DataFrame:
    """Return the number of proteins of every window that has one, by window id."""
    counts: dict[str, int] = {}
    for name, _ in queries:
        contig = name.rpartition("|")[0]
        counts[contig] = counts.get(contig, 0) + 1
    return pd.DataFrame(sorted(counts.items()), columns=["contig_id", "n_genes"])


def write_queries(queries: Sequence[tuple[str, str]], path: Path) -> None:
    """Write the proteins of the windows as a FASTA file."""
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in queries:
            handle.write(f">{name}\n{sequence}\n")


def diamond_binary(requested: Path | None) -> str:
    """Return the DIAMOND binary: the requested one, the environment's or on PATH."""
    if requested is not None:
        return str(requested)
    beside = Path(sys.executable).parent / "diamond"
    if beside.is_file():
        return str(beside)
    found = shutil.which("diamond")
    if found is None:
        msg = "diamond was not found next to the interpreter or on PATH; use --diamond"
        raise FileNotFoundError(msg)
    return found


def diamond_search(
    binary: str, db: Path, queries: Path, out: Path, targets: int, threads: int
) -> None:
    """Run DIAMOND with the settings of the training searches."""
    command = [
        binary,
        "blastp",
        "--db",
        str(db),
        "--query",
        str(queries),
        "--out",
        str(out),
        "--outfmt",
        "6",
        *DIAMOND_COLUMNS,
        "bitscore",
        "--max-target-seqs",
        str(targets),
        "--evalue",
        "1e-5",
        "--mid-sensitive",
        "--threads",
        str(threads),
        "--tmpdir",
        str(out.parent),
    ]
    done = subprocess.run(command, check=False, capture_output=True, text=True)
    if done.returncode != 0:
        msg = f"diamond failed ({done.returncode}) for {db}: {done.stderr.strip()}"
        raise RuntimeError(msg)


def add_reference_hits(
    features: pd.DataFrame,
    queries: Sequence[tuple[str, str]],
    index: pd.DataFrame,
    resources: layout.ScreenResources,
    settings: Settings,
    work: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join the best-hit features against the whole-proteome reference (no hold-out).

    Also returns the best hit of every gene and class, which the set model reads.
    """
    faa = work / "queries.faa"
    write_queries(queries, faa)
    binary = diamond_binary(settings.diamond)
    hits, hits_x = work / "hits.tsv", work / "hitsx.tsv"
    ref_db = resources.path(layout.REF_DMND).with_suffix("")
    refx_db = resources.path(layout.REFX_DMND).with_suffix("")
    diamond_search(binary, ref_db, faa, hits, 30, settings.threads)
    diamond_search(binary, refx_db, faa, hits_x, 100, settings.threads)
    ref = load_reference(
        resources.path(layout.REF_GENOMES_FILE),
        resources.path(layout.REFX_GENOMES_FILE),
        SKIPPED_CLASSES,
    )
    r_features, best = reference_features(ref, [hits, hits_x], index)
    return features.merge(r_features, on="contig_id", how="left"), best


def token_features(features: pd.DataFrame) -> pd.DataFrame:
    """Return the columns of the windows that ``hit_tokens`` reads."""
    names = ["contig_id", "genome_id", "class", "length_kb", "holdout", *SCALARS[2:]]
    return features[[name for name in names if name in features.columns]]


def add_set_logits(
    features: pd.DataFrame,
    best: pd.DataFrame,
    queries: Sequence[tuple[str, str]],
    index: pd.DataFrame,
    resources: layout.ScreenResources,
    settings: Settings,
) -> pd.DataFrame:
    """Join ``set_logit``: the logit of the set model of the hit tokens and embeddings.

    The tokens come from the best hits of ``add_reference_hits``, the embeddings are
    those of every protein of the windows, and the models are the fitted set models
    of the bundle.
    """
    from src.contig_screen import set_model  # noqa: PLC0415  # imports torch

    set_model.set_threads(settings.threads)
    genes, offsets, raw, pieces = window_tokens(token_features(features), best, index)
    embedder = make_embedder(resources, settings, "mean")
    tokens = extend_tokens(queries, genes, offsets, pieces, embedder)
    logits = set_model.score_contigs(resources, SET_VARIANT, tokens, offsets, raw)
    table = pieces[["contig_id"]].assign(set_logit=logits)
    return features.merge(table, on="contig_id", how="left")


def open_text(path: Path) -> TextIO:
    """Open a FASTA file for reading text, gzipped or not."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def fasta_bp(path: Path) -> int:
    """Return the number of sequence characters of the whole input file."""
    total = 0
    with open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def contig_lengths(path: Path) -> dict[str, int]:
    """Return contig name to length for every contig of the FASTA file."""
    lengths: dict[str, int] = {}
    name = ""
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            else:
                lengths[name] += len(line.strip())
    return lengths


def summarize(
    result: pd.DataFrame,
    threshold: float,
    total: int,
    lengths: dict[str, int] | None,
    score: str,
) -> None:
    """Print the scored share of the input and the called share.

    Contigs under 3 kb are not scored; they count as not called, so the share is
    relative to the whole input. With ``lengths`` the rows are overlapping windows:
    a base pair is called when the mean score of the windows over it reaches the
    threshold, so overlaps are not counted twice.
    """
    if lengths is None:
        scored = int(result["length_bp"].sum())
        print(f"{len(result)} contigs scored, {scored} of {total} bp")
        called = int(result.loc[result["called_cellular"], "length_bp"].sum())
    else:
        print(f"{len(result)} overlapping windows scored, input {total} bp")
        windows = result.rename(columns={score: "p"})
        calls = sliding_calls(windows, lengths, threshold)
        called = sum(int(c.sum()) for c in calls.values())
    share = 100.0 * called / total
    print(f"called cellular: {called} bp ({share:.1f}% of the input)")
    if share >= MOSTLY_CELLULAR_PCT:
        print(
            f"note: {MOSTLY_CELLULAR_PCT:.0f}% or more of the input is called "
            "cellular; this is a review flag, not a contamination estimate (the "
            "bin may not be a giant-virus bin, or its family may be absent from "
            "training)"
        )


if __name__ == "__main__":
    raise SystemExit(main())
