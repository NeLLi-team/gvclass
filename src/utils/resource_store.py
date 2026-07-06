"""Runtime accessors for GVClass resource bundles.

The historical bundle stores labels as TSV files and marker references as
``database/faa/<marker>.faa``. Newer compact bundles may store those same tables
in Parquet and materialize legacy views on demand.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import tempfile
from pathlib import Path
from typing import Iterable, Sequence


class ResourceStore:
    """Resolve labels and marker FASTA paths across legacy and Parquet bundles."""

    LABEL_TABLES = {
        "labels.tsv": {
            "parquet": "labels.parquet",
            "columns": ("label_id", "taxonomy"),
            "header": False,
        },
        "inactive_labels.tsv": {
            "parquet": "inactive_labels.parquet",
            "columns": ("label_id", "taxonomy", "reason"),
            "header": True,
        },
        "aliases.tsv": {
            "parquet": "aliases.parquet",
            "columns": (
                "alias_id",
                "canonical_id",
                "reason",
                "alias_taxonomy",
                "canonical_taxonomy",
            ),
            "header": True,
        },
        "label_context.tsv": {
            "parquet": "label_context.parquet",
            "columns": ("label_id", "context_type", "context_value", "source"),
            "header": True,
        },
    }

    def __init__(self, database_path: Path | str):
        self.database_path = Path(database_path)

    @property
    def parquet_root(self) -> Path:
        return self.database_path / "parquet"

    @property
    def faa_parquet(self) -> Path:
        return self.parquet_root / "faa.parquet"

    @property
    def label_parquet_dir(self) -> Path:
        return self.parquet_root / "labels"

    def has_parquet_labels(self, filename: str = "labels.tsv") -> bool:
        table = self.LABEL_TABLES.get(filename)
        if table is None:
            return False
        return (self.label_parquet_dir / str(table["parquet"])).exists()

    def has_parquet_faa(self) -> bool:
        return self.faa_parquet.exists()

    def label_path(self, filename: str = "labels.tsv") -> Path:
        """Return a TSV label path, materializing from Parquet when necessary."""
        legacy_path = self.database_path / filename
        if legacy_path.exists():
            return legacy_path

        table = self.LABEL_TABLES.get(filename)
        if table is None:
            return legacy_path

        parquet_path = self.label_parquet_dir / str(table["parquet"])
        if not parquet_path.exists():
            return legacy_path

        return self._materialize_label_table(filename, parquet_path, table)

    def marker_faa_path(self, marker: str) -> Path:
        """Return a marker FASTA path, materializing from Parquet when necessary."""
        for candidate in self._legacy_marker_candidates(marker):
            if candidate.exists():
                return candidate

        if not self.faa_parquet.exists():
            raise FileNotFoundError(
                f"Database file not found for marker {marker}: "
                f"{self.database_path / 'database' / 'faa' / f'{marker}.faa'}"
            )

        return self._materialize_marker_faa(marker)

    def materialized_faa_dir(
        self, markers: Iterable[str], *, strict: bool = True
    ) -> Path:
        """Return a directory containing ``<marker>.faa`` for every marker.

        Legacy marker files are linked into the cache directory. Parquet markers
        are materialized there. This gives downstream code one ordinary directory
        without forcing the whole Parquet FAA table to expand.

        When ``strict`` is false, missing markers are left absent from the returned
        directory so downstream per-marker code can skip them individually.
        """
        cache_dir = self._faa_cache_dir()
        for marker in markers:
            try:
                source = self.marker_faa_path(marker)
            except FileNotFoundError:
                if strict:
                    raise
                continue
            target = cache_dir / self._marker_filename(marker)
            if self._same_existing_path(source, target):
                continue
            self._link_or_copy(source, target)
        return cache_dir

    @classmethod
    def labels_satisfied(cls, database_path: Path | str) -> bool:
        store = cls(database_path)
        return (
            store.database_path / "labels.tsv"
        ).exists() or store.has_parquet_labels()

    def _legacy_marker_candidates(self, marker: str) -> tuple[Path, Path]:
        filename = self._marker_filename(marker)
        return (
            self.database_path / "database" / "faa" / filename,
            self.database_path / "faa" / filename,
        )

    @staticmethod
    def _marker_filename(marker: str) -> str:
        return marker if marker.endswith(".faa") else f"{marker}.faa"

    def _materialize_label_table(
        self, filename: str, parquet_path: Path, table: dict
    ) -> Path:
        target = self._label_cache_dir() / filename
        if self._target_current(target, parquet_path, allow_empty=True):
            return target

        pq, _ds = self._pyarrow_modules()
        columns: Sequence[str] = table["columns"]
        data = pq.read_table(parquet_path, columns=list(columns)).to_pylist()
        tmp = target.with_name(f"{target.name}.{os.getpid()}.tmp")
        with tmp.open("w") as handle:
            if table["header"]:
                handle.write("\t".join(columns) + "\n")
            for row in data:
                handle.write(
                    "\t".join(
                        "" if row.get(column) is None else str(row.get(column))
                        for column in columns
                    )
                    + "\n"
                )
        os.replace(tmp, target)
        return target

    def _materialize_marker_faa(self, marker: str) -> Path:
        filename = self._marker_filename(marker)
        target = self._faa_cache_dir() / filename
        if self._target_current(target, self.faa_parquet, allow_empty=False):
            return target

        _pq, ds = self._pyarrow_modules()
        dataset = ds.dataset(self.faa_parquet, format="parquet")
        marker_filter = ds.field("marker") == marker
        if marker != filename:
            marker_filter = marker_filter | (ds.field("marker") == filename)
        table = dataset.to_table(
            columns=["seq_id", "sequence"],
            filter=marker_filter,
        )
        rows = table.to_pylist()
        if not rows:
            raise FileNotFoundError(
                f"Marker {marker} was not found in Parquet FAA table {self.faa_parquet}"
            )

        tmp = target.with_name(f"{target.name}.{os.getpid()}.tmp")
        with tmp.open("w") as handle:
            for row in rows:
                seq_id = row.get("seq_id")
                sequence = row.get("sequence")
                if not seq_id or not sequence:
                    continue
                handle.write(f">{seq_id}\n")
                seq = str(sequence)
                for start in range(0, len(seq), 80):
                    handle.write(seq[start : start + 80] + "\n")
        if tmp.stat().st_size == 0:
            tmp.unlink(missing_ok=True)
            raise FileNotFoundError(
                f"Marker {marker} has no usable sequences in {self.faa_parquet}"
            )
        os.replace(tmp, target)
        return target

    def _label_cache_dir(self) -> Path:
        path = self._cache_root() / "labels"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _faa_cache_dir(self) -> Path:
        path = self._cache_root() / "faa"
        path.mkdir(parents=True, exist_ok=True)
        return path

    def _cache_root(self) -> Path:
        env_root = os.environ.get("GVCLASS_RESOURCE_CACHE", "").strip()
        if env_root:
            root = Path(env_root).expanduser()
        else:
            root = self.database_path / ".gvclass_cache"

        try:
            root.mkdir(parents=True, exist_ok=True)
            probe = root / ".write_test"
            probe.write_text("")
            probe.unlink(missing_ok=True)
            return root
        except OSError:
            key = hashlib.sha256(
                str(self.database_path.resolve()).encode()
            ).hexdigest()[:16]
            fallback = Path(tempfile.gettempdir()) / "gvclass_resource_cache" / key
            fallback.mkdir(parents=True, exist_ok=True)
            return fallback

    @staticmethod
    def _target_current(target: Path, source: Path, allow_empty: bool) -> bool:
        if not target.exists():
            return False
        if not allow_empty and target.stat().st_size == 0:
            return False
        return target.stat().st_mtime >= source.stat().st_mtime

    @staticmethod
    def _same_existing_path(source: Path, target: Path) -> bool:
        if not target.exists():
            return False
        try:
            return source.resolve() == target.resolve()
        except OSError:
            return False

    @staticmethod
    def _link_or_copy(source: Path, target: Path) -> None:
        tmp = target.with_name(f"{target.name}.{os.getpid()}.tmp")
        tmp.unlink(missing_ok=True)
        try:
            tmp.symlink_to(source)
        except OSError:
            tmp.unlink(missing_ok=True)
            shutil.copyfile(source, tmp)
        os.replace(tmp, target)

    @staticmethod
    def _pyarrow_modules():
        try:
            import pyarrow.dataset as ds
            import pyarrow.parquet as pq
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "This GVClass database stores labels or reference FAA in Parquet. "
                "Install pyarrow or run GVClass through the pixi environment."
            ) from exc
        return pq, ds
