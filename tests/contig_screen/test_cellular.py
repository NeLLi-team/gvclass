"""Tests of the screen driver: scores, searches, summary and command line."""

from __future__ import annotations

import contextlib
import gzip
import importlib.util
import io
import stat
import sys
from collections.abc import Iterator
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest

from src.contig_screen import cellular, layout
from src.contig_screen.stack_context import stack_scores
from tests.contig_screen.conftest import BundleFactory, write_bundle


class SumModel:
    """Stand-in for the boosting model: P(cellular) is the first feature."""

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        return np.column_stack([1 - x[:, 0], x[:, 0]])


class SigmoidOfZ:
    """Stand-in for stage 2: the sigmoid of the first context feature (``z``)."""

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        p = 1 / (1 + np.exp(-x[:, 0]))
        return np.column_stack([1 - p, p])


def dump(obj: object) -> bytes:
    buffer = io.BytesIO()
    joblib.dump(obj, buffer)
    return buffer.getvalue()


def windows_frame(scores: list[float]) -> pd.DataFrame:
    n = len(scores)
    frame = pd.DataFrame(
        {
            "contig_id": [f"b|c1|{i * 5000}" for i in range(n)],
            "source_contig": "c1",
            "source_start": [i * 5000 for i in range(n)],
            "length_bp": 10_000,
            "s_n_genes": 8,
            "c_gc": np.linspace(0.4, 0.6, n),
            "s_coding_density": 0.85,
            "f1": scores,
            "set_logit": 0.0,
        }
    )
    for name in ("s_genes_per_kb", "s_gene_len_mean", "h_frac_nohit"):
        frame[name] = 1.0
    for name in ("h_share_cellular", "h_share_giant", "h_share_PHAGE"):
        frame[name] = 0.0
    return frame


def test_stage_1_is_the_stack_of_the_boosting_probability_and_the_set_logit() -> None:
    frame = windows_frame([0.2, 0.9, 0.5]).assign(set_logit=[-1.0, 2.0, 0.0])
    boosting = {"model": SumModel(), "columns": ["f1"]}

    scored = cellular.stage1_scores(frame, boosting)

    expected = stack_scores(np.array([0.2, 0.9, 0.5]), np.array([-1.0, 2.0, 0.0]))
    assert scored["p_cellular"].tolist() == pytest.approx(expected.tolist())
    assert scored["p_gbm"].tolist() == pytest.approx([0.2, 0.9, 0.5])


def test_stage_1_refuses_a_window_without_a_set_logit() -> None:
    frame = windows_frame([0.2, 0.9]).assign(set_logit=[0.0, float("nan")])

    with pytest.raises(ValueError, match="1 windows have no set logit"):
        cellular.stage1_scores(frame, {"model": SumModel(), "columns": ["f1"]})


def screen_bundle(tmp_path: Path) -> layout.ScreenResources:
    """Return a bundle with stand-in models and the thresholds of the tests."""
    norm = io.BytesIO()
    np.savez(norm, comp=np.array(["c_gc"]), mu=np.array([0.5]), sd=np.array([0.1]))
    columns = ["z", "log_len"]
    files = {
        layout.BOOSTING_FILE: dump({"model": SumModel(), "columns": ["f1"]}),
        layout.STAGE2_FILE: dump({"model": SigmoidOfZ(), "columns": columns}),
        layout.STAGE2_NORM_FILE: norm.getvalue(),
        layout.THRESHOLDS_FILE: b'{"stage1": 0.5, "stage2": 0.5}',
    }
    digest = write_bundle(tmp_path, files)
    return layout.ScreenResources.open(tmp_path, manifest_sha256=digest)


def test_score_bin_stacks_stage_1_and_rescores_in_the_context_of_the_bin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    resources = screen_bundle(tmp_path)
    frame = windows_frame([0.01, 0.5, 0.99])
    monkeypatch.setattr(cellular, "window_features", lambda *args: frame.copy())

    result = cellular.score_bin(tmp_path / "b.fna", resources, cellular.Settings())

    assert result.columns.tolist() == cellular.OUTPUT_COLUMNS
    assert result["p_gbm"].tolist() == pytest.approx([0.01, 0.5, 0.99])
    assert result["p_context"].tolist() == pytest.approx(
        result["p_cellular"].tolist(), abs=1e-5
    )


def test_score_bin_adds_a_missing_model_column_as_missing_values(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    resources = screen_bundle(tmp_path)

    def features(*_: object) -> pd.DataFrame:
        return windows_frame([0.3, 0.6]).rename(columns={"f1": "f1_renamed"})

    monkeypatch.setattr(cellular, "window_features", features)

    result = cellular.score_bin(tmp_path / "b.fna", resources, cellular.Settings())

    assert result["p_gbm"].isna().all()


def test_query_index_counts_the_proteins_of_each_window_sorted_by_id() -> None:
    queries = [("b|c|5000|0", "A"), ("b|c|0|0", "A"), ("b|c|0|1", "C")]

    index = cellular.query_index(queries)

    assert index["contig_id"].tolist() == ["b|c|0", "b|c|5000"]
    assert index["n_genes"].tolist() == [2, 1]


def test_write_queries_writes_one_fasta_record_per_protein(tmp_path: Path) -> None:
    path = tmp_path / "q.faa"

    cellular.write_queries([("w|0", "MKT"), ("w|1", "AAA")], path)

    assert path.read_text() == ">w|0\nMKT\n>w|1\nAAA\n"


def test_fasta_helpers_count_bases_and_contigs_of_plain_and_gzipped_files(
    tmp_path: Path,
) -> None:
    text = ">c1 first\nACGT\nAC\n>c2\nGGG\n"
    plain, zipped = tmp_path / "b.fna", tmp_path / "b.fna.gz"
    plain.write_text(text)
    with gzip.open(zipped, "wt") as handle:
        handle.write(text)

    for path in (plain, zipped):
        assert cellular.fasta_bp(path) == 9
        assert cellular.contig_lengths(path) == {"c1": 6, "c2": 3}


def called_frame(called: list[bool], column: str) -> pd.DataFrame:
    frame = windows_frame([0.5] * len(called))
    frame[column] = [0.9 if c else 0.1 for c in called]
    return frame.assign(called_cellular=called)


def test_summarize_counts_the_called_base_pairs_of_whole_contigs(
    capsys: pytest.CaptureFixture[str],
) -> None:
    frame = called_frame([True, False], "p_context")

    cellular.summarize(frame, 0.5, 40_000, None, "p_context")

    out = capsys.readouterr().out
    assert "2 contigs scored, 20000 of 40000 bp" in out
    assert "called cellular: 10000 bp (25.0% of the input)" in out


def test_summarize_does_not_count_the_overlap_of_windows_twice(
    capsys: pytest.CaptureFixture[str],
) -> None:
    frame = called_frame([True, True], "p_context")

    cellular.summarize(frame, 0.5, 20_000, {"c1": 15_000}, "p_context")

    out = capsys.readouterr().out
    assert "2 overlapping windows scored, input 20000 bp" in out
    assert "called cellular: 15000 bp (75.0% of the input)" in out


def test_summarize_flags_a_bin_that_is_mostly_called(
    capsys: pytest.CaptureFixture[str],
) -> None:
    frame = called_frame([True], "p_context")

    cellular.summarize(frame, 0.5, 10_000, None, "p_context")

    assert "review flag, not a contamination estimate" in capsys.readouterr().out


def write_fake_diamond(path: Path, body: str) -> Path:
    path.write_text(f"#!/bin/bash\n{body}\n")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def test_diamond_binary_prefers_the_request_then_the_interpreters_folder(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "bin").mkdir()
    beside = write_fake_diamond(tmp_path / "bin" / "diamond", "exit 0")
    monkeypatch.setattr(sys, "executable", str(tmp_path / "bin" / "python"))

    assert cellular.diamond_binary(Path("/x/diamond")) == "/x/diamond"
    assert cellular.diamond_binary(None) == str(beside)


def test_diamond_binary_fails_clearly_when_there_is_none(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python"))
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(FileNotFoundError, match="--diamond"):
        cellular.diamond_binary(None)


def test_a_failing_diamond_run_raises_with_its_message(tmp_path: Path) -> None:
    binary = write_fake_diamond(tmp_path / "diamond", "echo 'bad database' >&2; exit 3")

    with pytest.raises(RuntimeError, match="bad database"):
        cellular.diamond_search(
            str(binary), tmp_path / "db", tmp_path / "q.faa", tmp_path / "o.tsv", 5, 1
        )


def test_reference_hits_join_the_best_hit_features_and_return_the_best_table(
    tmp_path: Path, make_resources: BundleFactory
) -> None:
    ref = "index\tclass\tlineage\n0\tBAC\tf__A\n1\tEUK\tf__B\n"
    refx = "index\tclass\tlineage\n2\tNCLDV\tf__C\n"
    resources = make_resources(
        tmp_path / "bundle",
        {
            layout.REF_DMND: b"db",
            layout.REFX_DMND: b"dbx",
            layout.REF_GENOMES_FILE: ref.encode(),
            layout.REFX_GENOMES_FILE: refx.encode(),
        },
    )
    hits = tmp_path / "hits"
    hits.mkdir()
    (hits / "ref7.tsv").write_text("w|0\tr0_x\t50.0\t100\t100\t100\t1e-9\t120.0\n")
    (hits / "ref7x.tsv").write_text("w|1\tr2_x\t60.0\t100\t100\t100\t1e-9\t150.0\n")
    diamond = write_fake_diamond(
        tmp_path / "diamond",
        "while [ $# -gt 0 ]; do\n"
        "  case $1 in --db) db=$2;; --out) out=$2;; esac; shift\n"
        "done\n"
        f'cp "{hits}/$(basename "$db").tsv" "$out"',
    )
    queries = [("w|0", "MKT"), ("w|1", "AAA")]
    features = pd.DataFrame({"contig_id": ["w"], "x": [1.0]})
    work = tmp_path / "work"
    work.mkdir()

    joined, best = cellular.add_reference_hits(
        features,
        queries,
        cellular.query_index(queries),
        resources,
        cellular.Settings(diamond=diamond, threads=1),
        work,
    )

    assert joined["r_share_BAC"].tolist() == [0.5]
    assert joined["r_share_NCLDV"].tolist() == [0.5]
    assert best["cls"].tolist() == ["BAC", "NCLDV"]
    assert (work / "queries.faa").read_text() == ">w|0\nMKT\n>w|1\nAAA\n"


def run_main(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *extra: str
) -> tuple[int, pd.DataFrame | None]:
    """Run ``main`` on a bundle with scored windows from a stand-in ``score_bin``."""
    root = tmp_path / "bundle"
    root.mkdir()
    digest = write_bundle(
        root, {layout.THRESHOLDS_FILE: b'{"stage1": 0.5, "stage2": 0.6}'}
    )
    fasta = tmp_path / "b.fna"
    fasta.write_text(">c1\n" + "A" * 20_000 + "\n")
    scored = windows_frame([0.2, 0.4, 0.6]).assign(
        p_gbm=0.5, p_cellular=[0.2, 0.55, 0.9], p_context=[0.1, 0.7, 0.6]
    )
    monkeypatch.setattr(
        cellular, "score_bin", lambda *a, **k: scored[cellular.OUTPUT_COLUMNS].copy()
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
    out = tmp_path / "out" / "result.tsv"

    code = cellular.main(
        [
            str(fasta),
            "--resources",
            str(root),
            "--manifest-sha256",
            digest,
            "--out",
            str(out),
            *extra,
        ]
    )

    return code, pd.read_csv(out, sep="\t") if out.exists() else None


def test_main_calls_windows_at_the_stage_2_threshold_of_the_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    code, table = run_main(tmp_path, monkeypatch)

    assert code == 0
    assert table is not None
    assert table["called_cellular"].tolist() == [False, True, True]
    out = capsys.readouterr().out
    assert "stage 1 alone (no bin context):" in out


def test_main_takes_a_threshold_from_the_command_line_before_the_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, table = run_main(tmp_path, monkeypatch, "--threshold2", "0.65")

    assert table is not None
    assert table["called_cellular"].tolist() == [False, True, False]


def test_main_writes_the_documented_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _, table = run_main(tmp_path, monkeypatch)

    assert table is not None
    assert table.columns.tolist() == [*cellular.OUTPUT_COLUMNS, "called_cellular"]


def test_main_reports_a_bundle_that_does_not_match_its_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(layout, "MANIFEST_SHA256", None)
    (tmp_path / "bundle").mkdir()
    write_bundle(tmp_path / "bundle", {layout.THRESHOLDS_FILE: b"{}"})

    code = cellular.main(
        [
            "x.fna",
            "--resources",
            str(tmp_path / "bundle"),
            "--out",
            str(tmp_path / "o.tsv"),
        ]
    )

    assert code == 2
    assert "no manifest digest" in capsys.readouterr().err


def test_main_reports_missing_torch_before_it_reads_anything(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: None)

    code = cellular.main(
        ["x.fna", "--resources", str(tmp_path), "--out", str(tmp_path / "o.tsv")]
    )

    assert code == 2
    assert "pixi" in capsys.readouterr().err


class SigmoidOfFirst:
    """Stand-in boosting model: the sigmoid of the first feature."""

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        p = 1 / (1 + np.exp(-x[:, 0]))
        return np.column_stack([1 - p, p])


def window_scores(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    contigs: dict[str, str],
    stem: str,
    sample_seed: str | None = None,
) -> pd.Series:
    """Return ``p_gbm`` per window content, with gene calls and sampling real.

    The searches, the label table and the ESM model are replaced by stand-ins; the
    protein sample, which decides the ``e_`` features, is the code under test.
    """
    from src.contig_screen import extract_genes as eg  # noqa: PLC0415
    from src.contig_screen.extra_features import sequence_hash  # noqa: PLC0415
    from tests.contig_screen.test_extra_features import HashEmbedder  # noqa: PLC0415

    work = tmp_path / stem
    work.mkdir()
    fasta = work / f"{stem}.fna"
    fasta.write_text("".join(f">{n}\n{s}\n" for n, s in contigs.items()))
    resources = layout.ScreenResources.open(
        work,
        manifest_sha256=write_bundle(
            work, {layout.LABELS_FILE: b"l", layout.LABEL_DOMAIN_FILE: b"d"}
        ),
    )
    windows: list[dict[str, object]] = []

    def extract(genome, _resources, _out, *, tile_kb, step_kb, threads):
        windows.extend(eg.tile_contigs(genome, tile_kb, step_kb))
        return windows, eg.call_genes(genome, windows)[0]

    def unit_features(*_: object) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "contig_id": [w["contig_id"] for w in windows],
                "content": [sequence_hash(str(w["sequence"])) for w in windows],
            }
        )

    monkeypatch.setattr(cellular, "extract_bin", extract)
    monkeypatch.setattr(cellular, "load_labels", lambda *_: None)
    monkeypatch.setattr(cellular, "unit_features", unit_features)
    monkeypatch.setattr(
        cellular, "add_reference_hits", lambda features, *_: (features, None)
    )
    monkeypatch.setattr(cellular, "make_embedder", lambda *_: HashEmbedder())
    monkeypatch.setattr(
        cellular, "add_set_logits", lambda features, *_: features.assign(set_logit=0.0)
    )
    boosting = {
        "model": SigmoidOfFirst(),
        "columns": ["e_0"],
        "embedding_pooling": "mean",
    }
    settings = cellular.Settings(threads=1, sample_seed=sample_seed)
    features = cellular.window_features(fasta, resources, settings, boosting, work)
    scored = cellular.stage1_scores(features, boosting)
    return scored.set_index("content")["p_gbm"].sort_index()


def test_p_gbm_is_the_same_for_a_renamed_and_reordered_copy_of_a_bin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from tests.contig_screen.test_extract_genes import coding_dna  # noqa: PLC0415

    contigs = {"c1": coding_dna(35, 1), "c2": coding_dna(35, 2)}
    copy = {"zz": contigs["c2"], "aa": contigs["c1"]}

    original = window_scores(tmp_path, monkeypatch, contigs, "first")
    renamed = window_scores(tmp_path, monkeypatch, copy, "second")
    legacy_a = window_scores(tmp_path, monkeypatch, contigs, "third", "NAME")
    legacy_b = window_scores(tmp_path, monkeypatch, copy, "fourth", "NAME")

    pd.testing.assert_series_equal(original, renamed, check_names=False)
    assert len(set(original.round(9))) > 1
    assert not np.allclose(legacy_a.to_numpy(), legacy_b.to_numpy())


def test_main_scores_inside_the_thread_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[int] = []

    @contextlib.contextmanager
    def spy(threads: int) -> Iterator[None]:
        seen.append(threads)
        yield

    monkeypatch.setattr(cellular.torch_stack, "limit_threads", spy)

    code, _ = run_main(tmp_path, monkeypatch, "--threads", "3")

    assert code == 0
    assert seen == [3]


def test_main_refuses_a_bundle_whose_manifest_is_not_the_pinned_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    monkeypatch.setattr(importlib.util, "find_spec", lambda name: object())
    (tmp_path / "bundle").mkdir()
    write_bundle(tmp_path / "bundle", {layout.THRESHOLDS_FILE: b"{}"})

    code = cellular.main(
        [
            "x.fna",
            "--resources",
            str(tmp_path / "bundle"),
            "--out",
            str(tmp_path / "o.tsv"),
        ]
    )

    assert code == 2
    assert "manifest SHA-256 mismatch" in capsys.readouterr().err
