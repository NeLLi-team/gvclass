"""Tests for src/utils/error_handling.py."""

import logging

from src.utils.error_handling import ErrorHandler, ValidationError


def test_log_warning_appends_context(caplog):
    handler = ErrorHandler(logger_name="gvclass_test_log_warning")

    with caplog.at_level(logging.WARNING, logger="gvclass_test_log_warning"):
        handler.log_warning("careful", context={"step": "unit"})

    assert "careful | Context: {'step': 'unit'}" in caplog.text


def test_validation_error_str_names_code_and_context():
    error = ValidationError("bad input", field="filepath", value="x.fna")

    assert str(error) == (
        "[VALIDATION_ERROR] bad input "
        "(Context: {'field': 'filepath', 'value': 'x.fna'})"
    )
