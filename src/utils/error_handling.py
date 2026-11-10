"""Exception hierarchy and warning logger for the GVClass pipeline."""

import logging
import sys
from typing import Any


class GVClassError(Exception):
    """Base exception for GVClass errors."""

    def __init__(
        self,
        message: str,
        error_code: str = "GVCLASS_ERROR",
        context: dict[str, Any] | None = None,
    ):
        super().__init__(message)
        self.message = message
        self.error_code = error_code
        self.context = context or {}

    def __str__(self) -> str:
        context_str = f" (Context: {self.context})" if self.context else ""
        return f"[{self.error_code}] {self.message}{context_str}"


class ValidationError(GVClassError):
    """Input validation error."""

    def __init__(
        self, message: str, field: str | None = None, value: object | None = None
    ):
        context = {}
        if field:
            context["field"] = field
        if value is not None:
            context["value"] = str(value)
        super().__init__(message, "VALIDATION_ERROR", context)


class ProcessingError(GVClassError):
    """Processing pipeline error."""

    def __init__(
        self, message: str, step: str | None = None, input_file: str | None = None
    ):
        context = {}
        if step:
            context["step"] = step
        if input_file:
            context["input_file"] = input_file
        super().__init__(message, "PROCESSING_ERROR", context)


class FileError(GVClassError):
    """File operation error."""

    def __init__(
        self,
        message: str,
        file_path: str | None = None,
        operation: str | None = None,
    ):
        context = {}
        if file_path:
            context["file_path"] = file_path
        if operation:
            context["operation"] = operation
        super().__init__(message, "FILE_ERROR", context)


class ErrorHandler:
    """Logs warnings with structured context to a named logger."""

    def __init__(self, logger_name: str = "gvclass", log_level: str = "INFO"):
        self.logger = logging.getLogger(logger_name)
        self.logger.setLevel(getattr(logging, log_level.upper()))

        if not self.logger.handlers:
            console_handler = logging.StreamHandler(sys.stderr)
            console_handler.setLevel(logging.WARNING)
            console_handler.setFormatter(
                logging.Formatter(
                    "%(asctime)s - %(name)s - %(levelname)s - "
                    "%(filename)s:%(lineno)d - %(message)s"
                )
            )
            self.logger.addHandler(console_handler)

    def log_warning(self, message: str, context: dict[str, Any] | None = None) -> None:
        """Log a warning, appending ``context`` when given."""
        if context:
            message += f" | Context: {context}"
        self.logger.warning(message)


error_handler = ErrorHandler()
