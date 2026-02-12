"""CLI entry point for HMM search utilities."""

import os

import click
import pandas as pd

from src.core.hmm_search import (
    get_models,
    load_cutoffs,
    process_hmmout,
    run_pyhmmer_search,
)
from src.utils.common import GVClassError, validate_file_path


def _run_search_if_needed(
    queryfaa: str, modelscombined: str, hmmout: str, threads: int, sensitive: bool
) -> None:
    query_file = validate_file_path(queryfaa, must_exist=True)
    models_file = validate_file_path(modelscombined, must_exist=True)
    output_file = validate_file_path(hmmout, must_exist=False)
    if not isinstance(threads, int) or threads < 1:
        raise ValueError("Threads must be a positive integer")

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print(f"Skipping hmmsearch as the output file {output_file} already exists.")
        return

    run_pyhmmer_search(
        str(models_file),
        str(query_file),
        str(output_file),
        threads,
        sensitive_mode=sensitive,
    )


def _write_empty_outputs(models: list[str], countout: str, scoreout: str) -> None:
    print("No hits found in the hmmsearch output. Generating empty output files.")
    pd.DataFrame(columns=["genomeid"] + models).set_index("genomeid").to_csv(
        countout, sep="\t"
    )
    pd.DataFrame(columns=["genomeid"] + models).set_index("genomeid").to_csv(
        scoreout, sep="\t"
    )


def _write_populated_outputs(count_matrix, score_matrix, countout: str, scoreout: str) -> None:
    count_matrix.reset_index().rename(columns={"index": "genomeid"}).to_csv(
        countout, sep="\t", index=False
    )
    score_matrix.reset_index().rename(columns={"index": "genomeid"}).to_csv(
        scoreout, sep="\t", index=False
    )


@click.command()
@click.option(
    "--queryfaa",
    "-q",
    type=click.Path(exists=True),
    required=True,
    help="Input query FASTA file",
)
@click.option(
    "--modelscombined",
    "-m",
    type=click.Path(exists=True),
    required=True,
    help="Combined HMM model file",
)
@click.option(
    "--hmmout", "-h", type=click.Path(), required=True, help="Output file for HMM hits"
)
@click.option(
    "--filtered_hmmout",
    "-hf",
    type=click.Path(),
    required=True,
    help="Output file for HMM hits filtered",
)
@click.option(
    "--countout",
    "-c",
    type=click.Path(),
    required=True,
    help="Output file for marker gene counts",
)
@click.option(
    "--scoreout",
    "-s",
    type=click.Path(),
    required=True,
    help="Output file for highest bitscore matrix",
)
@click.option(
    "--cutoff_file",
    "-f",
    type=click.Path(exists=True),
    help="Cutoff file for filtering HMM hits",
)
@click.option("--threads", "-t", default=4, help="Number of CPU threads to use")
@click.option(
    "--sensitive",
    is_flag=True,
    help="Use sensitive search mode (E-value 1e-5) instead of GA cutoffs",
)
def main(
    queryfaa: str,
    modelscombined: str,
    hmmout: str,
    filtered_hmmout: str,
    countout: str,
    scoreout: str,
    cutoff_file: str = None,
    threads: int = 4,
    sensitive: bool = False,
) -> None:
    """Run HMM search and post-process output matrices."""
    try:
        _run_search_if_needed(queryfaa, modelscombined, hmmout, threads, sensitive)
    except (ValueError, FileNotFoundError, GVClassError) as exc:
        print(f"Error: {exc}")
        raise click.ClickException(str(exc))

    models = get_models(modelscombined)
    cutoffs = (
        load_cutoffs(cutoff_file) if cutoff_file else {model: (0.0, 0.0) for model in models}
    )
    count_matrix, score_matrix = process_hmmout(models, hmmout, filtered_hmmout, cutoffs)

    if count_matrix.empty or score_matrix.empty:
        _write_empty_outputs(models, countout, scoreout)
        return

    _write_populated_outputs(count_matrix, score_matrix, countout, scoreout)


if __name__ == "__main__":
    main()
