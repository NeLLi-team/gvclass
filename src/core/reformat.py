import os
from pathlib import Path

import pandas as pd
from Bio.SeqRecord import SeqRecord

from src.utils.error_handling import ValidationError


def reformat_sequences(records: list[SeqRecord], input_file: str) -> list[SeqRecord]:
    """Rename every record to ``<input basename>|<original ID>``.

    The original ID is kept whole so that distinct IDs stay distinct. A leading
    ``<basename>|`` is dropped first, so GVClass output used as input keeps its
    IDs. Empty fields at either end are dropped, every other ``|`` becomes
    ``__``, and ``:`` becomes ``__``. Descriptions are cleared.

    Raises:
        ValidationError: If an ID has no content, or if two records end up with
            the same new ID. The error names both original IDs.
    """
    basename = Path(input_file).stem
    originals_by_new_id: dict[str, str] = {}
    for record in records:
        new_id = f"{basename}|{_flatten_id(record.id, basename)}"
        if new_id in originals_by_new_id:
            raise ValidationError(
                f"Sequence IDs {originals_by_new_id[new_id]!r} and {record.id!r} "
                f"in {input_file} both become {new_id!r}. Rename one of them.",
                field="sequence_id",
                value=record.id,
            )
        originals_by_new_id[new_id] = record.id
    for record, new_id in zip(records, originals_by_new_id, strict=True):
        record.id = new_id
        record.description = ""
    return records


def _flatten_id(sequence_id: str, basename: str) -> str:
    """Drop the basename prefix and empty end fields; replace ``|`` and ``:``."""
    prefix = f"{basename}|"
    unprefixed = sequence_id.removeprefix(prefix)
    fields = unprefixed.split("|")
    while fields and not fields[0]:
        fields.pop(0)
    while fields and not fields[-1]:
        fields.pop()
    if not fields:
        raise ValidationError(
            f"Sequence ID {sequence_id!r} has no content",
            field="sequence_id",
            value=sequence_id,
        )
    return "__".join(fields).replace(":", "__")


def calculate_stats(records: list[SeqRecord], input_file: str) -> pd.DataFrame:
    """Calculate statistics for the input FASTA file.

    Args:
        records (list[SeqRecord]): List of sequence records.
        input_file (str): Path to the input FASTA file.

    Returns:
        pd.DataFrame: DataFrame containing the calculated statistics.
    """
    total_length = sum(len(record.seq) for record in records)
    gc_count = sum(
        str(record.seq).upper().count("G") + str(record.seq).upper().count("C")
        for record in records
    )
    gc_percentage = (gc_count / total_length * 100) if total_length > 0 else 0
    stats_dict = {
        "query": os.path.splitext(os.path.basename(input_file))[0],
        "contigs": len(records),
        "LENbp": total_length,
        "GCperc": gc_percentage,
        "genecount": len(records),
        "CODINGperc": "no_fna",
        "ttable": "no_fna",
    }
    return pd.DataFrame([stats_dict])
