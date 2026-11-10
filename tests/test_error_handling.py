"""Tests for src/utils/error_handling.py."""

from src.utils.error_handling import ValidationError


def test_validation_error_str_names_code_and_context():
    error = ValidationError("bad input", field="filepath", value="x.fna")

    assert str(error) == (
        "[VALIDATION_ERROR] bad input "
        "(Context: {'field': 'filepath', 'value': 'x.fna'})"
    )


def test_library_loggers_attach_no_handler_of_their_own():
    from src.core import marker_processing
    from src.utils import input_validation

    assert marker_processing.logger.handlers == []
    assert input_validation.logger.handlers == []


def test_sha256_file_reports_every_chunk_it_hashes(tmp_path):
    import hashlib

    from src.utils.common import sha256_file

    path = tmp_path / "blob.bin"
    payload = b"x" * 2_500_000
    path.write_bytes(payload)
    sizes: list[int] = []

    digest = sha256_file(path, on_chunk=sizes.append)

    assert digest == hashlib.sha256(payload).hexdigest()
    assert sum(sizes) == len(payload)
    assert len(sizes) == 3
