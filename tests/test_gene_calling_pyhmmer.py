from pathlib import Path

import pyhmmer
import pytest

from src.core.gene_calling import run_hmmsearch


MODEL_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQANLQDKPEAIRRLVDELAQLGVE"


@pytest.fixture
def tiny_amino_hmm(tmp_path: Path) -> Path:
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequence = pyhmmer.easel.TextSequence(
        name=b"tiny_model",
        sequence=MODEL_PROTEIN,
    ).digitize(alphabet)
    hmm, _, _ = pyhmmer.plan7.Builder(alphabet).build(
        sequence,
        pyhmmer.plan7.Background(alphabet),
    )

    hmm_path = tmp_path / "tiny.hmm"
    with hmm_path.open("wb") as handle:
        hmm.write(handle)
    return hmm_path


def test_run_hmmsearch_handles_ambiguous_first_protein(
    tmp_path: Path,
    tiny_amino_hmm: Path,
) -> None:
    faa_path = tmp_path / "ambiguous_first.faa"
    faa_path.write_text(
        ">ambiguous_first\n"
        "MGKKNRKKNNNDNNNDNKNDNKNDNKN\n"
        ">matching_protein\n"
        f"{MODEL_PROTEIN}\n"
    )

    hit_count, *_ = run_hmmsearch(str(faa_path), str(tiny_amino_hmm))

    assert hit_count >= 1
    assert "matching_protein" in faa_path.with_suffix(".hmmout").read_text()


def test_run_hmmsearch_forces_amino_for_acgt_only_protein(
    tmp_path: Path,
    tiny_amino_hmm: Path,
) -> None:
    faa_path = tmp_path / "acgt_only.faa"
    faa_path.write_text(">acgt_only\n" + "ACGT" * 20 + "\n")

    run_hmmsearch(str(faa_path), str(tiny_amino_hmm))

    header = faa_path.with_suffix(".hmmout").read_text().splitlines()[0]
    assert header.startswith("# target name\taccession\t")
