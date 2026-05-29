"""Tests for src/utils/error_handling.py (H4: ErrorHandler.log_error)."""

import logging

import pytest

from src.utils.error_handling import ErrorHandler, error_handler, retry_on_failure


def test_log_error_exists_and_logs(caplog):
    handler = ErrorHandler(logger_name="gvclass_test_log_error")
    with caplog.at_level(logging.ERROR, logger="gvclass_test_log_error"):
        handler.log_error("boom", context={"step": "unit"})
    assert "boom" in caplog.text
    assert "step" in caplog.text


def test_retry_on_failure_reraises_original(caplog):
    calls = {"n": 0}

    @retry_on_failure(max_retries=2, delay=0)
    def always_fails():
        calls["n"] += 1
        raise ValueError("original failure")

    # Must raise the ORIGINAL exception after retries, not AttributeError from a
    # missing log_error method.
    with caplog.at_level(logging.ERROR, logger="gvclass"):
        with pytest.raises(ValueError, match="original failure"):
            always_fails()

    assert calls["n"] == 3  # initial + 2 retries
    assert "All 3 attempts failed" in caplog.text


def test_global_error_handler_has_log_error():
    # The live call sites in hmm_search.py / hmm_search_multi.py use the module
    # singleton; it must expose log_error.
    assert hasattr(error_handler, "log_error")
