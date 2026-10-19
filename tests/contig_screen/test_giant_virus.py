"""Tests of the giant-virus screen entry point with a small synthetic bin."""

from __future__ import annotations

import contextlib
import hashlib
import io
from collections.abc import Iterator
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import HistGradientBoostingClassifier

from src.contig_screen import giant_virus, layout
from src.contig_screen.gv_genes import GENE_COLUMNS, GeneCalls
from src.contig_screen.gv_order import N_CLASSES, Priors
from tests.contig_screen.conftest import write_bundle

COLUMNS = ["dup_frac", "pair_shared"]
CONTIG_BP = 6000
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"


def model_bytes() -> bytes:
    """Return a bundle whose score is high for contigs that share their markers."""
    features = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 1.0], [1.0, 0.0]] * 3)
    model = HistGradientBoostingClassifier(min_samples_leaf=1, max_iter=20)
    model.fit(features, [0, 1, 0, 1] * 3)
    priors = Priors(
        {"OG1": 0.01, "OG2": 0.01},
        0.1,
        np.zeros(N_CLASSES + 1),
        np.ones(N_CLASSES + 1),
        {"codon": (np.zeros(65), np.ones(65)), "struct": (np.zeros(4), np.ones(4))},
    )
    buffer = io.BytesIO()
    joblib.dump({"model": model, "columns": COLUMNS, "priors": priors}, buffer)
    return buffer.getvalue()


def bundle_files(*, with_cellular_cuts: bool = False) -> dict[str, bytes]:
    files = {
        layout.GV_MODEL_FILE: model_bytes(),
        layout.GV_THRESHOLDS_FILE: b'{"0.01": 0.5}',
        layout.HMM_FILE: b"hmm",
    }
    if with_cellular_cuts:
        files[layout.THRESHOLDS_FILE] = b'{"stage1": 0.4, "stage2": 0.3}'
    return files


def write_fasta(path: Path, names: list[str]) -> dict[str, str]:
    rng = np.random.default_rng(3)
    sequences = {n: "".join(rng.choice(list("ACGT"), size=CONTIG_BP)) for n in names}
    path.write_text("".join(f">{n}\n{s}\n" for n, s in sequences.items()))
    return sequences


def fake_calls(models: dict[str, list[str]]) -> GeneCalls:
    """Return a gene table with marker genes on the named contigs."""
    rows = []
    for contig, contig_models in models.items():
        for index, model in enumerate(contig_models):
            begin = 100 + 600 * index
            gene = (contig, index, begin, begin + 299, 1, 100, f"h-{model}", model)
            rows.append((*gene, 50.0, 1e-9, PROTEIN))
    return GeneCalls(pd.DataFrame(rows, columns=GENE_COLUMNS), [PROTEIN] * len(rows))


@pytest.fixture
def screen(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> tuple[Path, str, Path]:
    """Return the bundle root, its digest and a three-contig bin."""
    root = tmp_path / "bundle"
    digest = write_bundle(root, bundle_files())
    fasta = tmp_path / "bin.fna"
    write_fasta(fasta, ["c1", "c2", "c3"])
    calls = fake_calls({"c1": ["OG1", "OG2"], "c2": ["OG1", "OG2"], "c3": ["OG3"]})
    monkeypatch.setattr(giant_virus, "call_genes", lambda *args: calls)
    return root, digest, fasta


def run(root: Path, digest: str, fasta: Path, out: Path, *extra: str) -> int:
    argv = [str(fasta), "--resources", str(root), "--manifest-sha256", digest]
    return giant_virus.main([*argv, "--out", str(out), *extra])


def test_the_contig_table_has_the_prototype_columns_and_one_row_per_contig(
    screen: tuple[Path, str, Path], tmp_path: Path
) -> None:
    root, digest, fasta = screen
    out = tmp_path / "out.tsv"

    code = run(root, digest, fasta, out)

    table = pd.read_csv(out, sep="\t")
    assert code == 0
    assert list(table.columns) == giant_virus.KEPT_COLUMNS
    assert table["contig"].tolist() == ["c1", "c2", "c3"]
    assert table["p_foreign"].between(0, 1).all()


def test_contigs_that_share_markers_score_higher_and_are_called_at_the_threshold(
    screen: tuple[Path, str, Path],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root, digest, fasta = screen
    out = tmp_path / "out.tsv"

    run(root, digest, fasta, out)

    table = pd.read_csv(out, sep="\t").set_index("contig")
    assert table.loc["c1", "p_foreign"] > table.loc["c3", "p_foreign"]
    assert table["called_foreign"].tolist() == (table["p_foreign"] >= 0.5).tolist()
    assert "contigs scored: 3 of 3" in capsys.readouterr().out


def test_a_threshold_argument_replaces_the_one_in_the_bundle(
    screen: tuple[Path, str, Path], tmp_path: Path
) -> None:
    root, digest, fasta = screen
    out = tmp_path / "out.tsv"

    run(root, digest, fasta, out, "--threshold", "1.1")

    assert not pd.read_csv(out, sep="\t")["called_foreign"].any()


def test_a_bin_with_one_scored_contig_writes_no_table(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root = tmp_path / "bundle"
    digest = write_bundle(root, bundle_files())
    fasta = tmp_path / "bin.fna"
    write_fasta(fasta, ["c1"])
    calls = fake_calls({"c1": ["OG1"]})
    monkeypatch.setattr(giant_virus, "call_genes", lambda *args: calls)
    out = tmp_path / "out.tsv"

    code = run(root, digest, fasta, out)

    assert code == 0
    assert not out.exists()
    assert "fewer than two contigs" in capsys.readouterr().out


def test_a_wrong_manifest_digest_stops_the_run_with_exit_code_2(
    screen: tuple[Path, str, Path],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root, _, fasta = screen
    wrong = hashlib.sha256(b"other").hexdigest()

    code = run(root, wrong, fasta, tmp_path / "out.tsv")

    assert code == 2
    assert "mismatch" in capsys.readouterr().err


def test_the_union_adds_the_cellular_windows_to_the_called_contigs(
    screen: tuple[Path, str, Path],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root, digest, fasta = screen
    cellular = tmp_path / "cellular.tsv"
    windows = pd.DataFrame(
        {
            "source_contig": ["c3", "c3"],
            "source_start": [0, 3000],
            "length_bp": [3000, 3000],
            "p_context": [0.9, 0.1],
        }
    )
    windows.to_csv(cellular, sep="\t", index=False)

    run(
        root,
        digest,
        fasta,
        tmp_path / "out.tsv",
        "--threshold",
        "0.0",
        "--cellular",
        str(cellular),
        "--cellular-threshold",
        "0.5",
    )

    shown = capsys.readouterr().out
    assert "cellular calls: 3000 bp = 16.7%" in shown
    assert "union with giant-virus contigs: 18000 bp = 100.0%" in shown


def test_the_cellular_threshold_comes_from_the_bundle_when_it_has_one(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root = tmp_path / "bundle"
    digest = write_bundle(root, bundle_files(with_cellular_cuts=True))
    fasta = tmp_path / "bin.fna"
    write_fasta(fasta, ["c1", "c2"])
    calls = fake_calls({"c1": ["OG1"], "c2": ["OG2"]})
    monkeypatch.setattr(giant_virus, "call_genes", lambda *args: calls)
    cellular = tmp_path / "cellular.tsv"
    pd.DataFrame(
        {
            "source_contig": ["c1"],
            "source_start": [0],
            "length_bp": [CONTIG_BP],
            "p_context": [0.31],
        }
    ).to_csv(cellular, sep="\t", index=False)

    run(root, digest, fasta, tmp_path / "out.tsv", "--cellular", str(cellular))

    assert "cellular calls: 6000 bp" in capsys.readouterr().out


def test_without_any_cellular_threshold_the_union_is_not_computed(
    screen: tuple[Path, str, Path],
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    root, digest, fasta = screen
    cellular = tmp_path / "cellular.tsv"
    cellular.write_text("source_contig\tsource_start\tlength_bp\tp_context\n")

    run(root, digest, fasta, tmp_path / "out.tsv", "--cellular", str(cellular))

    assert "union not computed" in capsys.readouterr().out


def test_main_scores_inside_the_thread_limit(
    screen: tuple[Path, str, Path],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root, digest, fasta = screen
    seen: list[int] = []

    @contextlib.contextmanager
    def spy(threads: int) -> Iterator[None]:
        seen.append(threads)
        yield

    monkeypatch.setattr(giant_virus.torch_stack, "limit_threads", spy)

    code = run(root, digest, fasta, tmp_path / "out.tsv", "--threads", "2")

    assert code == 0
    assert seen == [2]
