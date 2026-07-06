"""Tests for species-tree Section 6: run-level handoff + combined orchestration.

Handoff tests cover manifest init (clearing stale scratch), atomic sidecar +
``.done`` marker semantics, and that combined mode reads only manifest-listed,
completed sidecars (so resume-skipped / stale queries are excluded). The
integration tests drive the real supermatrix -> tree -> placement path on tiny
fixtures (no reference DB): ``run_combined_species_tree`` for the opt-in combined
artifact, and ``_build_and_place`` for the per-query summary-merge that owns the
placement columns.
"""

from __future__ import annotations

import json
from pathlib import Path

import src.core.species_tree.config as st_config
import src.core.species_tree.orchestration as orch
from src.core.species_tree.config import SpeciesTreePanel
from src.core.species_tree.handoff import (
    SCRATCH_DIRNAME,
    init_handoff,
    read_sidecars,
    scratch_dir_for,
    write_sidecar,
)


# ---------------------------------------------------------------------------
# Handoff
# ---------------------------------------------------------------------------


def test_init_clears_stale_and_writes_manifest(tmp_path: Path) -> None:
    scratch = tmp_path / SCRATCH_DIRNAME
    scratch.mkdir()
    (scratch / "stale.txt").write_text("old run")
    handoff = init_handoff(tmp_path, ["q2", "q1"])
    assert not (scratch / "stale.txt").exists()
    assert json.loads(handoff.manifest_path.read_text())["queries"] == ["q1", "q2"]


def test_init_fresh_clears_all_species_tree_outputs(tmp_path: Path) -> None:
    st = tmp_path / "species_tree"
    (st / "old").mkdir(parents=True)
    (st / "old" / "old.treefile").write_text("(a,b);")
    init_handoff(tmp_path, ["q1"], resume=False)
    # Fresh run wipes the whole species_tree dir (no stale trees survive).
    assert not st.exists()


def test_init_resume_scopes_deletion_to_this_run(tmp_path: Path) -> None:
    """--resume clears only this run's query subdirs + combined artifacts, keeping
    resume-skipped queries' per-query trees on disk."""
    st = tmp_path / "species_tree"
    (st / "q1").mkdir(parents=True)
    (st / "q1" / "q1.treefile").write_text("(a,b);")        # this run -> rebuilt
    (st / "kept").mkdir()
    (st / "kept" / "kept.treefile").write_text("(c,d);")    # resume-skipped -> keep
    (st / "combined.treefile").write_text("(e,f);")          # stale combined -> clear
    (st / "_combined" / "NCLDV").mkdir(parents=True)         # multi-panel combined -> clear

    init_handoff(tmp_path, ["q1"], resume=True)

    assert not (st / "q1").exists()
    assert (st / "kept" / "kept.treefile").exists()
    assert not (st / "combined.treefile").exists()
    assert not (st / "_combined").exists()


def test_sidecar_roundtrip(tmp_path: Path) -> None:
    handoff = init_handoff(tmp_path, ["q1"])
    write_sidecar(handoff.scratch_dir, "q1", {"query_id": "q1", "x": 1})
    assert read_sidecars(handoff.scratch_dir) == {"q1": {"query_id": "q1", "x": 1}}


def test_read_ignores_sidecar_without_done_marker(tmp_path: Path) -> None:
    handoff = init_handoff(tmp_path, ["q1", "q2"])
    write_sidecar(handoff.scratch_dir, "q1", {"a": 1})
    # q2 wrote its JSON but crashed before the .done marker landed.
    (handoff.scratch_dir / "q2").mkdir()
    (handoff.scratch_dir / "q2" / "sidecar.json").write_text(json.dumps({"a": 2}))
    assert set(read_sidecars(handoff.scratch_dir)) == {"q1"}


def test_read_ignores_query_not_in_manifest(tmp_path: Path) -> None:
    handoff = init_handoff(tmp_path, ["q1"])
    write_sidecar(handoff.scratch_dir, "q1", {"a": 1})
    # qX is complete but not part of this run's manifest (resume/stale).
    write_sidecar(handoff.scratch_dir, "qX", {"a": 9})
    assert set(read_sidecars(handoff.scratch_dir)) == {"q1"}


def test_fresh_run_does_not_ingest_previous_run(tmp_path: Path) -> None:
    init_handoff(tmp_path, ["q1", "qOld"])
    write_sidecar(scratch_dir_for(tmp_path), "qOld", {"a": 1})
    # A fresh run with a different query set clears the scratch.
    handoff2 = init_handoff(tmp_path, ["q1"])
    write_sidecar(handoff2.scratch_dir, "q1", {"a": 2})
    sidecars = read_sidecars(handoff2.scratch_dir)
    assert set(sidecars) == {"q1"}


# ---------------------------------------------------------------------------
# Combined orchestration (real supermatrix -> tree -> placement, tiny fixtures)
# ---------------------------------------------------------------------------

_GRP_A = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLS"
_GRP_B = "MSEQNNTSGFLGKKVDLSSLTGKKVAVDASHALYQFLIAVRQEGGQLTN"
_GRP_C = "MAQVINTNSLSLLTQNNLNKSQSALGTAIERLSSGLRINSAKDDAAGQA"
_GROUPS = ("grpA", "grpB", "grpC")
_BASE = {"grpA": _GRP_A, "grpB": _GRP_B, "grpC": _GRP_C}


def _variant(seq: str, n: int) -> str:
    # Deterministic small divergence: swap n residues to 'V'.
    chars = list(seq)
    for i in range(n):
        chars[(i * 7 + 3) % len(chars)] = "V"
    return "".join(chars)


def _setup_db(tmp_path: Path) -> Path:
    db = tmp_path / "db"
    faa_dir = db / "database" / "faa"
    faa_dir.mkdir(parents=True)
    labels = []
    for group in _GROUPS:
        with open(faa_dir / f"{group}.faa", "w") as handle:
            for ref, nmut in (("R1", 0), ("R2", 4), ("R3", 8)):
                handle.write(f">NCLDV__{ref}|{group}\n{_variant(_BASE[group], nmut)}\n")
    for ref, order in (("R1", "Imitervirales"), ("R2", "Algavirales"), ("R3", "Asfuvirales")):
        labels.append(
            f"NCLDV__{ref}\tNCLDV|Nucleocytoviricota|Megaviricetes|{order}|fam|gen|sp_{ref}"
        )
    (db / "labels.tsv").write_text("\n".join(labels) + "\n")
    return db


def _sidecar_payload(query_id: str) -> dict:
    # Query sequences close to R1; neighbors list all three refs per group.
    neighbors = []
    for group in _GROUPS:
        for ref, dist in (("R1", 0.1), ("R2", 0.2), ("R3", 0.3)):
            neighbors.append(
                {
                    "genome_id": f"NCLDV__{ref}",
                    "protein_id": f"NCLDV__{ref}|{group}",
                    "group": group,
                    "distance": dist,
                }
            )
    return {
        "query_id": query_id,
        "panel": "NCLDV",
        "query_reps": {group: _variant(_BASE[group], 1) for group in _GROUPS},
        "neighbors": neighbors,
    }


def test_run_combined_species_tree_end_to_end(tmp_path: Path, monkeypatch) -> None:
    # Use a 3-group NCLDV panel (min_markers=2) so the tiny fixtures suffice.
    test_panel = SpeciesTreePanel(
        name="NCLDV",
        domain_prefix="NCLDV__",
        gate_token="d_NCLDV",
        groups=_GROUPS,
        min_markers=2,
    )
    monkeypatch.setitem(st_config.PANELS, "NCLDV", test_panel)

    db = _setup_db(tmp_path)
    output_base = tmp_path / "out"
    output_base.mkdir()

    init_handoff(output_base, ["q1", "q2"])
    scratch = scratch_dir_for(output_base)
    write_sidecar(scratch, "q1", _sidecar_payload("q1"))
    write_sidecar(scratch, "q2", _sidecar_payload("q2"))

    results = [
        {"query": "q1", "status": "complete", "summary_data": {"taxonomy_majority": "d_NCLDV"}},
        {"query": "q2", "status": "complete", "summary_data": {"taxonomy_majority": "d_NCLDV"}},
    ]

    orch.run_combined_species_tree(output_base, db, results, threads=1)

    # Combined is an additional ARTIFACT only: the tree + taxonomy table are
    # written, covering both queries.
    st_dir = output_base / "species_tree"
    assert (st_dir / "combined.treefile").exists()
    tax_tsv = st_dir / "species_tree_taxonomy.tsv"
    assert tax_tsv.exists()
    rows = [line.split("\t") for line in tax_tsv.read_text().splitlines()]
    header = rows[0]
    assert header[0] == "query" and "species_tree_nn_taxonomy" in header
    query_rows = {r[0] for r in rows[1:]}
    assert query_rows == {"q1", "q2"}

    # The combined step does NOT touch the summary columns — the per-query
    # placement owns them (see test_per_query_merges_summary_columns).
    for res in results:
        assert "species_tree_nn_genome" not in res["summary_data"]


def test_per_query_merges_summary_columns(tmp_path: Path, monkeypatch) -> None:
    """The per-query assembly path (the default product) merges placement columns
    into the query's summary_data."""
    test_panel = SpeciesTreePanel(
        name="NCLDV", domain_prefix="NCLDV__", gate_token="d_NCLDV",
        groups=_GROUPS, min_markers=2,
    )
    monkeypatch.setitem(st_config.PANELS, "NCLDV", test_panel)
    db = _setup_db(tmp_path)
    output_base = tmp_path / "out"
    output_base.mkdir()
    analyzer = orch.TreeAnalyzer(db / "labels.tsv")

    sd = {"taxonomy_majority": "d_NCLDV"}
    orch._build_and_place(
        test_panel,
        {"q1": _sidecar_payload("q1")},
        db,
        analyzer,
        tmp_path / "work",
        output_base / "species_tree" / "q1",
        threads=1,
        summary_by_query={"q1": sd},
        trim_method="pytrimal",
        k_neighbors=20,
        basename="q1",
    )

    assert sd["species_tree_nn_genome"].startswith("NCLDV__")
    assert "p_Nucleocytoviricota" in sd["species_tree_nn_taxonomy"]
    assert "species_tree_nn_distance" in sd
    # Per-query files are basenamed by the query id (not "combined").
    assert (output_base / "species_tree" / "q1" / "q1.treefile").exists()


def test_cli_argv_forwards_species_tree_flags() -> None:
    """The subprocess seam: the launcher must put --species-tree/
    --species-tree-combined on the runner argv when set, and omit them otherwise."""
    import argparse

    from src.bin.gvclass_cli import _append_optional_pipeline_flags

    base = dict(
        cluster_queue=None,
        cluster_project=None,
        cluster_walltime=None,
        verbose=False,
        resume=False,
        allow_short=False,
    )
    on = argparse.Namespace(
        species_tree=True, species_tree_combined=True, species_tree_trim="pytrimal", **base
    )
    cmd: list = []
    _append_optional_pipeline_flags(cmd, on, True, "legacy", False)
    assert "--species-tree" in cmd and "--species-tree-combined" in cmd
    assert "--species-tree-trim" in cmd and "pytrimal" in cmd

    off = argparse.Namespace(
        species_tree=False, species_tree_combined=False, species_tree_trim="witchi", **base
    )
    cmd_off: list = []
    _append_optional_pipeline_flags(cmd_off, off, True, "legacy", False)
    assert "--species-tree" not in cmd_off and "--species-tree-combined" not in cmd_off
    # Default trim is not forwarded (keeps argv minimal).
    assert "--species-tree-trim" not in cmd_off


def test_combined_excludes_failed_query(tmp_path: Path, monkeypatch) -> None:
    """A query whose sidecar landed but whose pipeline later FAILED is not placed."""
    test_panel = SpeciesTreePanel(
        name="NCLDV", domain_prefix="NCLDV__", gate_token="d_NCLDV",
        groups=_GROUPS, min_markers=2,
    )
    monkeypatch.setitem(st_config.PANELS, "NCLDV", test_panel)
    db = _setup_db(tmp_path)
    output_base = tmp_path / "out"
    output_base.mkdir()
    init_handoff(output_base, ["q1", "q2"])
    scratch = scratch_dir_for(output_base)
    write_sidecar(scratch, "q1", _sidecar_payload("q1"))
    write_sidecar(scratch, "q2", _sidecar_payload("q2"))  # sidecar landed...
    results = [
        {"query": "q1", "status": "complete", "summary_data": {"taxonomy_majority": "d_NCLDV"}},
        {"query": "q2", "status": "failed", "error": "post-process boom"},  # ...but failed
    ]
    orch.run_combined_species_tree(output_base, db, results, threads=1)

    tax = (output_base / "species_tree" / "species_tree_taxonomy.tsv").read_text()
    placed = {line.split("\t")[0] for line in tax.splitlines()[1:]}
    assert placed == {"q1"}  # q2 excluded (failed)
    assert "q2" not in tax


def test_combined_no_sidecars_is_noop(tmp_path: Path) -> None:
    db = _setup_db(tmp_path)
    output_base = tmp_path / "out"
    output_base.mkdir()
    init_handoff(output_base, [])  # no queries
    results = []
    # Must not raise and must not create a species_tree dir.
    orch.run_combined_species_tree(output_base, db, results, threads=1)
    assert not (output_base / "species_tree" / "combined.treefile").exists()
