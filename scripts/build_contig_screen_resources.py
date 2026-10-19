"""Assemble the ``contig_screen/`` resource tree and its manifest.

Run from the repository root:

    python -m scripts.build_contig_screen_resources --out-root <bundle root> \
        --ref-dir <ref7 dir> --boosting <model.joblib> --set-dir <set model dir> \
        --stage2-dir <stage 2 dir> --threshold1 <T1> --threshold2 <T2> \
        --esm-dir <ESM-2 files> --label-domain <label_domain.parquet> \
        --refs-dir <marker references> [--hmm combined.hmm] [--labels labels.parquet] \
        [--gv-dir <giant-virus model dir>]

The tree follows ``src/contig_screen/layout.py``. ``manifest.json`` lists the relative
path, size and SHA-256 of every file. The script prints the digest of the manifest,
which belongs in ``layout.MANIFEST_SHA256``. ``--link`` links the inputs instead of
copying them (for a scratch bundle). The HMM catalog and the label table are the
bundle's own files: pass them to place them at ``hmm/combined.hmm`` and
``parquet/labels/labels.parquet``, or leave them out when the bundle root has them.
``--gv-dir`` (a folder with ``model.joblib`` and ``thresholds.json``) adds the files
of the giant-virus screen; its other inputs are the catalog files of ``--ref-dir`` and
the HMM catalog, which the tree holds once.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from src.contig_screen import layout

logger = logging.getLogger(__name__)

REF_FILES = {
    "ref7.dmnd": layout.REF_DMND,
    "ref7x.dmnd": layout.REFX_DMND,
    "ref_genomes.tsv": layout.REF_GENOMES_FILE,
    "ref_genomes_x.tsv": layout.REFX_GENOMES_FILE,
}
DEFAULT_SET_VARIANT = "hits_struct"


@dataclass(frozen=True, slots=True)
class Inputs:
    """Where the files of the screen come from."""

    ref_dir: Path
    boosting: Path
    set_dir: Path
    stage2_dir: Path
    esm_dir: Path
    label_domain: Path
    refs_dir: Path
    threshold1: float
    threshold2: float
    set_variant: str = DEFAULT_SET_VARIANT
    hmm: Path | None = None
    labels: Path | None = None
    gv_dir: Path | None = None


def build_resources(inputs: Inputs, out_root: Path, *, link: bool = False) -> str:
    """Write the ``contig_screen/`` tree under ``out_root``; return the manifest digest.

    Raises:
        FileNotFoundError: An input file is missing, or the bundle root lacks the
            HMM catalog or the label table.
    """
    placed: list[str] = []
    for source, relative in planned_files(inputs):
        place(source, out_root / relative, link=link)
        placed.append(relative)
    thresholds = {"stage1": inputs.threshold1, "stage2": inputs.threshold2}
    thresholds_path = out_root / layout.THRESHOLDS_FILE
    thresholds_path.parent.mkdir(parents=True, exist_ok=True)
    thresholds_path.write_text(
        json.dumps(thresholds, indent=1) + "\n", encoding="utf-8"
    )
    placed.append(layout.THRESHOLDS_FILE)
    placed.extend(shared_files(out_root))
    manifest = layout.build_manifest(out_root, placed)
    (out_root / layout.MANIFEST_FILE).write_bytes(manifest)
    digest = hashlib.sha256(manifest).hexdigest()
    check_complete(out_root, digest)
    return digest


def planned_files(inputs: Inputs) -> list[tuple[Path, str]]:
    """Return (source, relative destination) for every file the tree copies."""
    files = [(inputs.ref_dir / name, rel) for name, rel in REF_FILES.items()]
    files.append((inputs.boosting, layout.BOOSTING_FILE))
    files.append((inputs.stage2_dir / "stage2.joblib", layout.STAGE2_FILE))
    files.append((inputs.stage2_dir / "norm.npz", layout.STAGE2_NORM_FILE))
    files.append((inputs.label_domain, layout.LABEL_DOMAIN_FILE))
    files.append((inputs.refs_dir / "model_to_group.json", layout.MODEL_TO_GROUP_FILE))
    files.extend(
        (path, f"{layout.MARKER_REFS_DIR}/{path.name}")
        for path in sorted(inputs.refs_dir.glob("*.faa"))
    )
    variant = inputs.set_variant
    set_files = [
        *sorted(inputs.set_dir.glob(f"model_{variant}_k*.pt")),
        inputs.set_dir / f"scaler_{variant}.json",
    ]
    files.extend((path, f"{layout.SET_MODEL_DIR}/{path.name}") for path in set_files)
    files.extend(
        (path, f"{layout.ESM_DIR}/{path.name}")
        for path in sorted(inputs.esm_dir.iterdir())
        if path.is_file()
    )
    if inputs.gv_dir is not None:
        files.append((inputs.gv_dir / "model.joblib", layout.GV_MODEL_FILE))
        files.append((inputs.gv_dir / "thresholds.json", layout.GV_THRESHOLDS_FILE))
    if inputs.hmm is not None:
        files.append((inputs.hmm, layout.HMM_FILE))
    if inputs.labels is not None:
        files.append((inputs.labels, layout.LABELS_FILE))
    return files


def shared_files(out_root: Path) -> list[str]:
    """Return the bundle's own files that the manifest pins; they must exist."""
    shared = [layout.HMM_FILE, layout.LABELS_FILE]
    for relative in shared:
        if not (out_root / relative).is_file():
            msg = f"{out_root / relative} is missing: pass --hmm and --labels"
            raise FileNotFoundError(msg)
    return shared


def place(source: Path, destination: Path, *, link: bool) -> None:
    """Copy or link ``source`` to ``destination``, replacing what is there."""
    if not source.is_file():
        msg = f"input file not found: {source}"
        raise FileNotFoundError(msg)
    if destination.exists() and destination.resolve() == source.resolve():
        msg = f"{destination} is the same file as its source"
        raise ValueError(msg)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.exists():
        destination.unlink()
    if link:
        destination.symlink_to(source.resolve())
    else:
        shutil.copy2(source, destination)


def check_complete(out_root: Path, digest: str) -> None:
    """Open the new tree as the screen does and check that every file is listed."""
    resources = layout.ScreenResources.open(out_root, manifest_sha256=digest)
    absent = [rel for rel in layout.FIXED_FILES if not resources.has(rel)]
    directories = (layout.MARKER_REFS_DIR, layout.SET_MODEL_DIR, layout.ESM_DIR)
    absent.extend(f"{d}/" for d in directories if not resources.names(d))
    if resources.has(layout.GV_MODEL_FILE):
        absent.extend(rel for rel in layout.GV_FILES if not resources.has(rel))
    if absent:
        msg = f"the manifest lacks {absent}"
        raise FileNotFoundError(msg)


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    """Return the parsed command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--ref-dir", type=Path, required=True)
    parser.add_argument("--boosting", type=Path, required=True)
    parser.add_argument("--set-dir", type=Path, required=True)
    parser.add_argument("--set-variant", default=DEFAULT_SET_VARIANT)
    parser.add_argument("--stage2-dir", type=Path, required=True)
    parser.add_argument("--threshold1", type=float, required=True)
    parser.add_argument("--threshold2", type=float, required=True)
    parser.add_argument("--esm-dir", type=Path, required=True)
    parser.add_argument("--label-domain", type=Path, required=True)
    parser.add_argument("--refs-dir", type=Path, required=True)
    parser.add_argument("--hmm", type=Path, default=None)
    parser.add_argument("--labels", type=Path, default=None)
    parser.add_argument(
        "--gv-dir",
        type=Path,
        default=None,
        help="folder with model.joblib and thresholds.json of the giant-virus screen",
    )
    parser.add_argument("--link", action="store_true", help="link, do not copy")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Build the tree and print the manifest digest."""
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    inputs = Inputs(
        ref_dir=args.ref_dir,
        boosting=args.boosting,
        set_dir=args.set_dir,
        stage2_dir=args.stage2_dir,
        esm_dir=args.esm_dir,
        label_domain=args.label_domain,
        refs_dir=args.refs_dir,
        threshold1=args.threshold1,
        threshold2=args.threshold2,
        set_variant=args.set_variant,
        hmm=args.hmm,
        labels=args.labels,
        gv_dir=args.gv_dir,
    )
    digest = build_resources(inputs, args.out_root, link=args.link)
    print(f"MANIFEST_SHA256 = {digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
