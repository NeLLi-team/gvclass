"""Layout and integrity checks of the resources of the cellular contig screen.

The screen reads every input from one GVClass resource bundle root. The files it
owns sit under ``contig_screen/``; the HMM catalog and the label table are the
bundle's own files and are shared with the rest of GVClass::

    <root>/hmm/combined.hmm
    <root>/parquet/labels/labels.parquet
    <root>/contig_screen/
        manifest.json                  relative path, size and SHA-256 of every file
        thresholds.json                {"stage1": T1, "stage2": T2}
        label_domain.parquet           label id to domain
        ref/ref7.dmnd, ref7x.dmnd      DIAMOND databases (proteomes; giant viruses)
        ref/ref_genomes.tsv, ref_genomes_x.tsv
        refs/<marker>.faa, model_to_group.json   marker references (pyswrd search)
        models/boosting.joblib         stage 1 boosting model
        models/set/model_<variant>_k<k>.pt, scaler_<variant>.json   set model
        models/stage2/stage2.joblib, norm.npz                       bin-context model
        esm2_t6_8M_UR50D/              local copy of the Hugging Face files

``manifest.json`` lists every file above, with paths relative to the bundle root.
Its own SHA-256 is pinned in :data:`MANIFEST_SHA256`. A pickle, joblib, torch or
JSON file is read only after its bytes match the manifest; a file of 256 MiB or
more, such as a DIAMOND database, is checked by size unless ``verify_all`` is set.
"""

from __future__ import annotations

import hashlib
import io
import json
import logging
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any

import joblib
import numpy as np

from src.utils.common import sha256_file

logger = logging.getLogger(__name__)

#: SHA-256 of ``contig_screen/manifest.json``. The release build sets it; a code
#: constant (not a sidecar file) means that a swapped bundle needs a reviewed code
#: change, as for ``CONTAMINATION_MODEL_SHA256``. ``None`` makes the screen refuse
#: to run unless the caller passes the digest explicitly. The value is that of the
#: 10% reference bundle of release 2.1.0 (``scripts/build_contig_screen_resources.py``).
MANIFEST_SHA256: str | None = (
    "51be1dfd061676e27f42576b850ea9fcd6c58a515581417b18818158be36f08c"
)

SCREEN_DIR = "contig_screen"
MANIFEST_FILE = f"{SCREEN_DIR}/manifest.json"
THRESHOLDS_FILE = f"{SCREEN_DIR}/thresholds.json"
LABEL_DOMAIN_FILE = f"{SCREEN_DIR}/label_domain.parquet"
REF_DMND = f"{SCREEN_DIR}/ref/ref7.dmnd"
REFX_DMND = f"{SCREEN_DIR}/ref/ref7x.dmnd"
REF_GENOMES_FILE = f"{SCREEN_DIR}/ref/ref_genomes.tsv"
REFX_GENOMES_FILE = f"{SCREEN_DIR}/ref/ref_genomes_x.tsv"
MARKER_REFS_DIR = f"{SCREEN_DIR}/refs"
MODEL_TO_GROUP_FILE = f"{MARKER_REFS_DIR}/model_to_group.json"
BOOSTING_FILE = f"{SCREEN_DIR}/models/boosting.joblib"
SET_MODEL_DIR = f"{SCREEN_DIR}/models/set"
STAGE2_FILE = f"{SCREEN_DIR}/models/stage2/stage2.joblib"
STAGE2_NORM_FILE = f"{SCREEN_DIR}/models/stage2/norm.npz"
ESM_DIR = f"{SCREEN_DIR}/esm2_t6_8M_UR50D"
HMM_FILE = "hmm/combined.hmm"
LABELS_FILE = "parquet/labels/labels.parquet"
#: Files of the giant-virus screen (``giant_virus.py``). It also reads the catalog
#: DIAMOND database and table (``REFX_DMND``, ``REFX_GENOMES_FILE``) and the HMM
#: catalog (``HMM_FILE``) of the cellular screen; nothing is copied.
GV_DIR = f"{SCREEN_DIR}/giant_virus"
GV_MODEL_FILE = f"{GV_DIR}/model.joblib"
GV_THRESHOLDS_FILE = f"{GV_DIR}/thresholds.json"
GV_FILES = (
    GV_MODEL_FILE,
    GV_THRESHOLDS_FILE,
    REFX_DMND,
    REFX_GENOMES_FILE,
    HMM_FILE,
)

#: Files that the screen needs and the builder must list, apart from the
#: directories (marker references, set model, ESM-2 files).
FIXED_FILES = (
    THRESHOLDS_FILE,
    LABEL_DOMAIN_FILE,
    REF_DMND,
    REFX_DMND,
    REF_GENOMES_FILE,
    REFX_GENOMES_FILE,
    MODEL_TO_GROUP_FILE,
    BOOSTING_FILE,
    STAGE2_FILE,
    STAGE2_NORM_FILE,
    HMM_FILE,
    LABELS_FILE,
)
MANIFEST_VERSION = 1
HASH_SIZE_LIMIT = 256 * 1024 * 1024
#: Deserialized by Python (pickle-based or executable on load): always hashed.
ALWAYS_HASHED_SUFFIXES = frozenset({".joblib", ".pt", ".npz", ".pkl", ".bin", ".json"})


class IntegrityError(RuntimeError):
    """A resource file is missing, unlisted or does not match the manifest."""


@dataclass(frozen=True, slots=True)
class ManifestEntry:
    """Size and SHA-256 of one resource file."""

    size: int
    sha256: str


class ScreenResources:
    """The files of a resource bundle that the manifest lists, checked on access."""

    def __init__(
        self,
        root: Path,
        entries: Mapping[str, ManifestEntry],
        *,
        verify_all: bool = False,
    ) -> None:
        """Hold the bundle root and the manifest entries keyed by relative path."""
        self.root = root
        self.entries = dict(entries)
        self.verify_all = verify_all

    @classmethod
    def open(
        cls,
        root: Path,
        *,
        manifest_sha256: str | None = None,
        verify_all: bool = False,
    ) -> ScreenResources:
        """Read the manifest of ``root`` after checking its digest.

        ``manifest_sha256`` replaces the pinned :data:`MANIFEST_SHA256`; it is for
        bundles that are not the release bundle. The digest must still match.

        Raises:
            IntegrityError: No digest is pinned or given, the manifest is missing,
                or its digest differs.
        """
        expected = manifest_sha256 or MANIFEST_SHA256
        if expected is None:
            msg = (
                "no manifest digest is pinned (layout.MANIFEST_SHA256 is unset): "
                "pass manifest_sha256 / --manifest-sha256 for this bundle"
            )
            raise IntegrityError(msg)
        manifest_path = root / MANIFEST_FILE
        if not manifest_path.is_file():
            msg = f"contig screen manifest not found: {manifest_path}"
            raise IntegrityError(msg)
        raw = manifest_path.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != expected.lower():
            msg = (
                f"manifest SHA-256 mismatch at {manifest_path}: expected "
                f"{expected}, got {digest}; refusing to use the resources"
            )
            raise IntegrityError(msg)
        return cls(root, parse_manifest(raw), verify_all=verify_all)

    def has(self, relative: str) -> bool:
        """Return True when the manifest lists ``relative``."""
        return relative in self.entries

    def names(self, directory: str) -> list[str]:
        """Return the listed paths directly inside ``directory``, sorted."""
        return sorted(
            rel for rel in self.entries if str(PurePosixPath(rel).parent) == directory
        )

    def path(self, relative: str) -> Path:
        """Return the path of a listed file after the size (and hash) check."""
        entry = self._entry(relative)
        path = self.root / relative
        size = self._size(path, relative)
        if size != entry.size:
            msg = f"{path} has {size} bytes, the manifest lists {entry.size}"
            raise IntegrityError(msg)
        if self._needs_hash(relative, entry):
            self._check_digest(path, sha256_file(path), entry)
        return path

    def read_bytes(self, relative: str) -> bytes:
        """Return the bytes of a listed file, read once and checked against the hash.

        Deserialize these bytes, not the path, so that a swap of the file between
        the check and the load cannot reach the loader.
        """
        entry = self._entry(relative)
        path = self.root / relative
        raw = path.read_bytes()
        self._check_digest(path, hashlib.sha256(raw).hexdigest(), entry)
        return raw

    def directory(self, relative: str) -> Path:
        """Return a directory whose listed files all pass the checks."""
        names = self.names(relative)
        if not names:
            msg = f"the manifest lists no file under {relative}"
            raise IntegrityError(msg)
        for name in names:
            self.path(name)
        folder = self.root / relative
        listed = set(names)
        extra = sorted(
            path.relative_to(self.root).as_posix()
            for path in folder.rglob("*")
            if path.is_file() and path.relative_to(self.root).as_posix() not in listed
        )
        if extra:
            msg = f"{relative} holds files that the manifest does not list: {extra}"
            raise IntegrityError(msg)
        return folder

    def _entry(self, relative: str) -> ManifestEntry:
        if relative not in self.entries:
            msg = f"{relative} is not listed in the contig screen manifest"
            raise IntegrityError(msg)
        return self.entries[relative]

    @staticmethod
    def _size(path: Path, relative: str) -> int:
        if not path.is_file():
            msg = f"contig screen resource missing: {path} ({relative})"
            raise IntegrityError(msg)
        return path.stat().st_size

    def _needs_hash(self, relative: str, entry: ManifestEntry) -> bool:
        suffix = PurePosixPath(relative).suffix
        return (
            self.verify_all
            or suffix in ALWAYS_HASHED_SUFFIXES
            or entry.size < HASH_SIZE_LIMIT
        )

    @staticmethod
    def _check_digest(path: Path, digest: str, entry: ManifestEntry) -> None:
        if digest != entry.sha256:
            msg = (
                f"SHA-256 mismatch at {path}: the manifest lists {entry.sha256}, "
                f"the file has {digest}; refusing to load it"
            )
            raise IntegrityError(msg)


def parse_manifest(raw: bytes) -> dict[str, ManifestEntry]:
    """Return the entries of a manifest, refusing unsafe paths."""
    document = json.loads(raw)
    if document.get("version") != MANIFEST_VERSION:
        msg = f"unsupported contig screen manifest version {document.get('version')!r}"
        raise IntegrityError(msg)
    entries: dict[str, ManifestEntry] = {}
    for relative, item in document["files"].items():
        parts = PurePosixPath(relative).parts
        if PurePosixPath(relative).is_absolute() or ".." in parts or not parts:
            msg = f"unsafe path in the contig screen manifest: {relative!r}"
            raise IntegrityError(msg)
        entries[relative] = ManifestEntry(int(item["size"]), str(item["sha256"]))
    return entries


def build_manifest(root: Path, relative_paths: Iterable[str]) -> bytes:
    """Return the manifest bytes for the listed files under ``root``."""
    files = {}
    for relative in sorted(set(relative_paths)):
        path = root / relative
        files[relative] = {"size": path.stat().st_size, "sha256": sha256_file(path)}
    document = {"version": MANIFEST_VERSION, "files": files}
    return (json.dumps(document, indent=1, sort_keys=True) + "\n").encode("utf-8")


def load_joblib(resources: ScreenResources, relative: str) -> Any:  # noqa: ANN401  # joblib objects are untyped
    """Return the object of a verified joblib file."""
    return joblib.load(io.BytesIO(resources.read_bytes(relative)))


def load_npz(resources: ScreenResources, relative: str) -> dict[str, np.ndarray]:
    """Return the arrays of a verified ``.npz`` file (no pickled objects)."""
    with np.load(io.BytesIO(resources.read_bytes(relative)), allow_pickle=False) as f:
        return {name: f[name] for name in f.files}


def load_thresholds(resources: ScreenResources) -> dict[str, float]:
    """Return the stage 1 and stage 2 thresholds from ``thresholds.json``."""
    document = json.loads(resources.read_bytes(THRESHOLDS_FILE))
    return {"stage1": float(document["stage1"]), "stage2": float(document["stage2"])}
