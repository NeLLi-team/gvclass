"""Tests of the set model scoring (needs torch)."""

from __future__ import annotations

import io
import json
import pickle
from pathlib import Path

import numpy as np
import pytest

from src.contig_screen import hit_tokens as ht, layout
from tests.contig_screen.conftest import write_bundle

torch = pytest.importorskip("torch")
sm = pytest.importorskip("src.contig_screen.set_model")

N_FEATURES = len(ht.GENE_FEATURES)
N_SCALARS = 2 + len(sm.STRUCT_COLUMNS)
HIDDEN = 8


def make_net(seed: int) -> sm.SetNet:
    torch.manual_seed(seed)
    return sm.SetNet(N_FEATURES, N_SCALARS, HIDDEN).eval()


def toy_tokens() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(0)
    genes = rng.random((5, N_FEATURES), dtype=np.float32)
    offsets = np.array([0, 3, 3, 5])
    scalars = rng.random((3, N_SCALARS), dtype=np.float32)
    return genes, offsets, scalars


def state_bytes(net: sm.SetNet) -> bytes:
    buffer = io.BytesIO()
    torch.save(net.state_dict(), buffer)
    return buffer.getvalue()


def scaler_json(n_models: int) -> bytes:
    meta = {
        "fill": [0.0] * N_SCALARS,
        "mean": [0.0] * N_SCALARS,
        "std": [1.0] * N_SCALARS,
        "n_scalars": N_SCALARS,
        "n_features": N_FEATURES,
        "ablation": "full",
        "scalar_set": "hits_struct",
        "hidden": HIDDEN,
        "n_models": n_models,
    }
    return json.dumps(meta).encode()


def test_batch_genes_collects_the_tokens_of_the_chosen_contigs_in_order() -> None:
    genes = np.arange(10, dtype=np.float32).reshape(5, 2)
    offsets = np.array([0, 2, 3, 5])

    tokens, piece = sm.batch_genes(genes, offsets, np.array([2, 0]))

    assert tokens[:, 0].tolist() == [6.0, 8.0, 0.0, 2.0]
    assert piece.tolist() == [0, 0, 1, 1]


def test_a_window_without_genes_is_scored_and_does_not_disturb_the_others() -> None:
    genes, offsets, scalars = toy_tokens()
    net = make_net(1)

    together = sm.predict(net, genes, offsets, scalars, np.array([0, 1, 2]))
    alone = sm.predict(net, genes, offsets, scalars, np.array([2]))

    assert np.isfinite(together).all()
    assert together[2] == pytest.approx(alone[0], abs=1e-6)


def test_the_score_of_a_window_does_not_depend_on_the_order_of_its_genes() -> None:
    genes, offsets, scalars = toy_tokens()
    net = make_net(2)
    shuffled = genes.copy()
    shuffled[:3] = genes[[2, 0, 1]]

    first = sm.predict(net, genes, offsets, scalars, np.array([0]))
    second = sm.predict(net, shuffled, offsets, scalars, np.array([0]))

    assert first[0] == pytest.approx(second[0], abs=1e-6)


def test_scalar_matrix_uses_log_length_and_genes_per_kb() -> None:
    raw = np.zeros((1, len(ht.SCALARS)), dtype=np.float32)
    raw[0, :2] = [10.0, 5.0]

    plain = sm.scalar_matrix(raw, "hits")
    struct = sm.scalar_matrix(raw, "hits_struct")

    assert plain.shape == (1, 2)
    assert plain[0].tolist() == pytest.approx([1.0, 0.5])
    assert struct.shape == (1, N_SCALARS)


def test_the_scaler_fills_missing_values_before_it_standardizes() -> None:
    scaler = sm.Scaler(
        np.array([4.0, 0.0], dtype=np.float32),
        np.array([2.0, 1.0], dtype=np.float32),
        np.array([2.0, 1.0], dtype=np.float32),
    )

    out = scaler.transform(np.array([[np.nan, 3.0]], dtype=np.float32))

    assert out.tolist() == [[1.0, 2.0]]


def test_embedding_scale_standardizes_present_genes_and_zeroes_the_others() -> None:
    scale = sm.EmbeddingScale(np.array([1.0, 1.0]), np.array([2.0, 2.0]))
    width = sm.EMBEDDING_START
    genes = np.zeros((2, width + 3), dtype=np.float32)
    genes[0, width : width + 2] = [3.0, 5.0]
    genes[0, -1] = 1.0
    genes[1, width : width + 2] = [9.0, 9.0]

    out = scale.transform(genes)

    assert out[0, width : width + 2].tolist() == [1.0, 2.0]
    assert not out[1, width:].any()


def test_scoring_a_saved_variant_averages_the_logits_of_its_models(
    tmp_path: Path,
) -> None:
    nets = [make_net(3), make_net(4)]
    base = layout.SET_MODEL_DIR
    files = {f"{base}/scaler_hits_struct.json": scaler_json(2)}
    files.update(
        {
            f"{base}/model_hits_struct_k{k}.pt": state_bytes(n)
            for k, n in enumerate(nets)
        }
    )
    digest = write_bundle(tmp_path, files)
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)
    genes, offsets, _ = toy_tokens()
    raw = np.random.default_rng(1).random((3, len(ht.SCALARS))).astype(np.float32)

    logits = sm.score_contigs(resources, "hits_struct", genes, offsets, raw)

    scaled = sm.Scaler(
        np.zeros(N_SCALARS, np.float32),
        np.zeros(N_SCALARS, np.float32),
        np.ones(N_SCALARS, np.float32),
    ).transform(sm.scalar_matrix(raw, "hits_struct"))
    expected = np.mean(
        [sm.predict(n, genes, offsets, scaled, np.arange(3)) for n in nets], axis=0
    )
    assert logits.tolist() == pytest.approx(expected.tolist(), abs=1e-6)


def test_a_model_file_that_differs_from_the_manifest_is_refused(
    tmp_path: Path,
) -> None:
    base = layout.SET_MODEL_DIR
    files = {
        f"{base}/scaler_hits_struct.json": scaler_json(1),
        f"{base}/model_hits_struct_k0.pt": state_bytes(make_net(5)),
    }
    digest = write_bundle(tmp_path, files)
    (tmp_path / f"{base}/model_hits_struct_k0.pt").write_bytes(state_bytes(make_net(6)))
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    with pytest.raises(layout.IntegrityError, match="SHA-256 mismatch"):
        sm.load_fitted(resources, "hits_struct")


class Exploit:
    """An object whose unpickling writes a marker file."""

    def __init__(self, marker: Path) -> None:
        self.marker = marker

    def __reduce__(self) -> tuple[object, tuple[str, str]]:
        return (Path.write_text, (str(self.marker), "ran"))


def test_a_listed_model_file_is_loaded_with_weights_only(tmp_path: Path) -> None:
    marker = tmp_path / "ran"
    base = layout.SET_MODEL_DIR
    buffer = io.BytesIO()
    torch.save({"w": Exploit(marker)}, buffer)
    files = {
        f"{base}/scaler_hits_struct.json": scaler_json(1),
        f"{base}/model_hits_struct_k0.pt": buffer.getvalue(),
    }
    digest = write_bundle(tmp_path, files)
    resources = layout.ScreenResources.open(tmp_path, manifest_sha256=digest)

    with pytest.raises(pickle.UnpicklingError):
        sm.load_fitted(resources, "hits_struct")
    assert not marker.exists()
