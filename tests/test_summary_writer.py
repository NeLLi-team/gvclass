from __future__ import annotations

import csv
import tarfile
from pathlib import Path


def test_failed_queries_are_reported_separately_from_final_summary(
    tmp_path: Path,
) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "completed_query",
                "status": "complete",
                "summary_data": {
                    "query": "completed_query",
                    "taxonomy_majority": "d_NCLDV;p_Nucleocytoviricota",
                    "taxonomy_confidence": "high",
                },
            },
            {
                "query": "failed_query",
                "status": "failed",
                "error": "simulated HMM search failure",
            },
        ],
        tmp_path,
    )

    with open(summary_file, newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert [row["query"] for row in summary_rows] == ["completed_query"]

    with open(tmp_path / "gvclass_failed_queries.tsv", newline="") as handle:
        failed_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert failed_rows == [
        {
            "query": "failed_query",
            "status": "failed",
            "error": "simulated HMM search failure",
        }
    ]
    assert (tmp_path / "gvclass_failed_queries.csv").exists()


def test_species_tree_sentinel_when_no_tree(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [{"query": "q1", "status": "complete", "summary_data": {"query": "q1"}}],
        tmp_path,
    )
    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["species_tree_nn_taxonomy"] == "nd"
    assert row["species_tree_nn_genome"] == "nd"
    assert row["species_tree_nn_distance"] == "nd"
    assert row["species_tree_clade_id"] == "nd"


def test_legacy_species_tree_sentinel_migrated_to_nd(tmp_path: Path) -> None:
    """A --resume over per-query files from an older gvclass version may carry the
    retired ``no-species-tree-calculated`` literal; it must migrate to ``nd``."""
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "species_tree_nn_taxonomy": "no-species-tree-calculated",
                },
            }
        ],
        tmp_path,
    )
    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["species_tree_nn_taxonomy"] == "nd"


def test_species_tree_taxonomy_value_passthrough(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_files

    summary_file = write_final_summary_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "species_tree_nn_taxonomy": "d_NCLDV;p_Nucleocytoviricota",
                },
            }
        ],
        tmp_path,
    )
    with open(summary_file, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["species_tree_nn_taxonomy"] == "d_NCLDV;p_Nucleocytoviricota"


def test_extended_summary_written(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import write_final_summary_extended_files

    extended = write_final_summary_extended_files(
        [
            {
                "query": "q1",
                "status": "complete",
                "summary_data": {
                    "query": "q1",
                    "cellular_coherent_contig_count": 0,
                    "contig_attribution_mode": "fna_gene_calling",
                },
            }
        ],
        tmp_path,
    )
    with open(extended, newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["contig_attribution_mode"] == "fna_gene_calling"
    assert row["cellular_coherent_contig_count"] == "0"


def test_extended_summary_archive_removes_loose_files(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import (
        archive_final_summary_extended_files,
        write_final_summary_extended_files,
    )

    write_final_summary_extended_files(
        [{"query": "q1", "status": "complete", "summary_data": {"query": "q1"}}],
        tmp_path,
    )

    archive = archive_final_summary_extended_files(tmp_path)

    assert archive == tmp_path / "gvclass_summary.extended.tar.gz"
    assert archive.exists()
    assert not (tmp_path / "gvclass_summary.extended.tsv").exists()
    assert not (tmp_path / "gvclass_summary.extended.csv").exists()
    with tarfile.open(archive, "r:gz") as tar_handle:
        assert {
            "gvclass_summary.extended.tsv",
            "gvclass_summary.extended.csv",
        }.issubset(set(tar_handle.getnames()))


def test_read_individual_summary_results_from_query_archive(tmp_path: Path) -> None:
    from src.pipeline.summary_writer import read_individual_summary_results

    work = tmp_path / "q1"
    work.mkdir()
    (work / "q1.summary.tab").write_text("query\ttaxonomy_majority\nq1\tlegacy\n")
    (work / "q1.final_summary.tsv").write_text("query\ttaxonomy_majority\nq1\tfinal\n")
    with tarfile.open(tmp_path / "q1.tar.gz", "w:gz") as tar_handle:
        tar_handle.add(work, arcname="q1")

    results = read_individual_summary_results(tmp_path, ["q1"])

    assert len(results) == 1
    assert results[0]["query"] == "q1"
    assert results[0]["summary_data"]["taxonomy_majority"] == "final"
