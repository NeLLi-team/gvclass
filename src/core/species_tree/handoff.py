"""Run-level handoff for the combined species tree (Section 6).

The per-query output directories are tarred and deleted by ``_post_process_query``
before the combined tree is built, so combined mode cannot read them back. Instead
each NCLDV query writes a small **sidecar** (its representative GVOG8 sequences +
placed neighbors) into a run-level scratch dir, and the combined step reads only
the sidecars listed in this run's manifest that carry a ``.done`` marker.

Atomicity: the manifest and each sidecar JSON are written ``tmp -> os.replace`` and
the ``.done`` marker is the LAST write, so a partially-written sidecar (or a worker
that crashed mid-write) is never trusted by the combined step, which runs in the
same process after all workers finish. Writes are rename-atomic but not
fsync-durable across an OS/power crash; that is sufficient here because combined
mode never reads sidecars across a process boundary and ``--resume`` re-runs from a
cleared scratch. Clearing the scratch dir at run start guarantees a fresh run never
ingests a previous run's queries (``--resume`` therefore covers only this run's
queries — a documented limitation).

The scratch dir (``<output_base>/.species_tree_work/``) is left in place after a run
as a dot-prefixed debugging artifact; it never enters a per-query archive and does
not affect a standard (no-``--species-tree``) run.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)

SCRATCH_DIRNAME = ".species_tree_work"
_MANIFEST = "manifest.json"
_SIDECAR = "sidecar.json"
_DONE = ".done"


@dataclass
class SpeciesTreeHandoff:
    scratch_dir: Path
    manifest_path: Path


def _atomic_write_json(path: Path, obj: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(obj))
    os.replace(tmp, path)


def init_handoff(
    output_base: Path, query_names: List[str], resume: bool = False
) -> SpeciesTreeHandoff:
    """Clear the scratch + prior species-tree outputs and write this run's manifest.

    On a fresh run the whole ``species_tree/`` output dir is cleared so a run that
    produces fewer/no placements never leaves a stale tree behind. On a ``--resume``
    run the deletion is SCOPED to this run's queries: only the per-query subdirs
    being (re)built plus the combined artifacts are removed, so resume-skipped
    queries keep their existing ``species_tree/<query>/`` deliverables.
    """
    outputs = output_base / "species_tree"
    if outputs.exists():
        if resume:
            for query_name in query_names:
                query_dir = outputs / query_name
                if query_dir.is_dir():
                    shutil.rmtree(query_dir)
            # Combined artifacts are always rebuilt from this run's sidecars, so
            # clear them regardless: every top-level file (single-panel combined.*)
            # plus the multi-panel _combined/ dir. Per-query subdirs are dirs, so
            # this never touches a skipped query's output.
            for child in outputs.iterdir():
                if child.is_file():
                    child.unlink()
            combined_dir = outputs / "_combined"
            if combined_dir.is_dir():
                shutil.rmtree(combined_dir)
        else:
            shutil.rmtree(outputs)
    scratch = output_base / SCRATCH_DIRNAME
    if scratch.exists():
        shutil.rmtree(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    manifest = scratch / _MANIFEST
    _atomic_write_json(manifest, {"queries": sorted(query_names)})
    logger.info("Species-tree handoff initialised for %d queries at %s", len(query_names), scratch)
    return SpeciesTreeHandoff(scratch_dir=scratch, manifest_path=manifest)


def scratch_dir_for(output_base: Path) -> Path:
    return output_base / SCRATCH_DIRNAME


def write_sidecar(scratch_dir: Path, query_id: str, payload: Dict) -> Path:
    """Write a query's sidecar atomically; the ``.done`` marker is written last."""
    query_dir = scratch_dir / query_id
    query_dir.mkdir(parents=True, exist_ok=True)
    _atomic_write_json(query_dir / _SIDECAR, payload)
    # .done is the completion marker: combined mode trusts a sidecar only when it
    # exists, so it must land strictly after the JSON (rename-atomic ordering).
    done = query_dir / _DONE
    done_tmp = done.with_name(_DONE + ".tmp")
    done_tmp.write_text("")
    os.replace(done_tmp, done)
    return query_dir / _SIDECAR


def read_sidecars(scratch_dir: Path) -> Dict[str, Dict]:
    """Read sidecars for manifest-listed queries that carry a ``.done`` marker.

    Queries absent from the manifest (e.g. resumed-skip from a prior run, or a
    stale leftover) are never read, so combined mode contains exactly this run's
    completed queries.
    """
    manifest_path = scratch_dir / _MANIFEST
    if not manifest_path.exists():
        return {}
    try:
        manifest = json.loads(manifest_path.read_text())
    except (json.JSONDecodeError, OSError) as exc:
        logger.warning("Unreadable species-tree manifest %s: %s", manifest_path, exc)
        return {}

    sidecars: Dict[str, Dict] = {}
    for query_id in manifest.get("queries", []):
        query_dir = scratch_dir / query_id
        if not (query_dir / _DONE).exists():
            continue
        sidecar_path = query_dir / _SIDECAR
        try:
            sidecars[query_id] = json.loads(sidecar_path.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning("Skipping unreadable sidecar %s: %s", sidecar_path, exc)
    return sidecars
