"""Read the tables of the contig screen without changing their identifiers."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import pandas as pd


def read_table(path: Path, id_columns: Sequence[str]) -> pd.DataFrame:
    """Read a tab-separated table; the identifier columns stay literal text.

    pandas turns ``007`` into 7 and ``NA`` into a missing value, so a contig of that
    name would fail the lookup of its length or drop out of a group, which reads as
    a clean bin. The other columns keep the usual parsing.

    Raises:
        ValueError: A column of ``id_columns`` is not in the table.
    """
    table = pd.read_csv(path, sep="\t")
    names = pd.read_csv(
        path, sep="\t", usecols=list(id_columns), dtype=str, keep_default_na=False
    )
    return table.assign(**{column: names[column] for column in id_columns})
