"""Validate the pinned example screen outputs before releasing the SIF."""

import csv
import math
import sys
from pathlib import Path

source, output = map(Path, sys.argv[1:])
bin_name = source.stem
with source.open(encoding="utf-8") as handle:
    input_contigs = [line[1:].split()[0] for line in handle if line.startswith(">")]
assert len(input_contigs) == len(set(input_contigs)) == 10


def read_table(name: str) -> list[dict[str, str]]:
    """Read a TSV from the mounted screen output directory."""
    with (output / name).open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


status = read_table("contig_screen_summary.tsv")
assert [
    {key: row[key] for key in ("bin", "cellular_status", "giant_virus_status")}
    for row in status
] == [{"bin": bin_name, "cellular_status": "ok", "giant_virus_status": "ok"}], status
cellular = read_table(f"{bin_name}.cellular.tsv")
giant_virus = read_table(f"{bin_name}.giant_virus.tsv")
assert len(cellular) == 33, len(cellular)
assert len(giant_virus) == 10, len(giant_virus)
assert len({row["contig_id"] for row in cellular}) == 33, "Duplicate cellular window"
assert len({(row["source_contig"], row["source_start"]) for row in cellular}) == 33, (
    "Duplicate window position"
)
assert {row["source_contig"] for row in cellular} == set(input_contigs), (
    "Cellular contig coverage differs"
)
assert {row["contig"] for row in giant_virus} == set(input_contigs), (
    "Giant-virus contig coverage differs"
)
for rows, columns in [
    (cellular, ("p_cellular", "p_gbm", "p_context")),
    (giant_virus, ("p_foreign",)),
]:
    for row in rows:
        for column in columns:
            probability = float(row[column])
            assert math.isfinite(probability), (column, row)
            assert 0 <= probability <= 1, (column, row)
print(
    "Verified both screens: 33 cellular windows, 10 giant-virus contigs, "
    "complete coverage and valid probabilities"
)
