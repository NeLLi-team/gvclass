"""
Error handling framework for GVClass workflow scripts.

This module provides structured error handling, logging, and recovery mechanisms
for the GVClass pipeline.
"""

import logging
import logging.handlers
import traceback
import sys
from typing import Optional, Dict, Any, Callable, TypeVar, Union
from pathlib import Path
from functools import wraps
from contextlib import contextmanager

# Type variable for decorated functions
F = TypeVar("F", bound=Callable[..., Any])


class GVClassError(Exception):
    """Base exception for GVClass errors."""

    def __init__(
        self,
        message: str,
        error_code: str = "GVCLASS_ERROR",
        context: Optional[Dict[str, Any]] = None,
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
        self, message: str, field: Optional[str] = None, value: Optional[Any] = None
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
        self, message: str, step: Optional[str] = None, input_file: Optional[str] = None
    ):
        context = {}
        if step:
            context["step"] = step
        if input_file:
            context["input_file"] = input_file
        super().__init__(message, "PROCESSING_ERROR", context)


class CommandError(GVClassError):
    """External command execution error."""

    def __init__(
        self,
        message: str,
        command: Optional[str] = None,
        return_code: Optional[int] = None,
    ):
        context = {}
        if command:
            context["command"] = command
        if return_code is not None:
            context["return_code"] = return_code
        super().__init__(message, "COMMAND_ERROR", context)


class FileError(GVClassError):
    """File operation error."""

    def __init__(
        self,
        message: str,
        file_path: Optional[str] = None,
        operation: Optional[str] = None,
    ):
        context = {}
        if file_path:
            context["file_path"] = file_path
        if operation:
            context["operation"] = operation
        super().__init__(message, "FILE_ERROR", context)


class ErrorHandler:
    """Centralized error handling and logging."""

    def __init__(self, logger_name: str = "gvclass", log_level: str = "INFO"):
        self.logger = self._setup_logger(logger_name, log_level)
        self.error_counts: Dict[str, int] = {}
        self.max_errors_per_type = 10

    def _setup_logger(self, name: str, level: str) -> logging.Logger:
        """Set up logger with structured formatting."""
        logger = logging.getLogger(name)
        logger.setLevel(getattr(logging, level.upper()))

        if not logger.handlers:
            # Console handler
            console_handler = logging.StreamHandler(sys.stderr)
            console_handler.setLevel(logging.WARNING)

            # Detailed formatter for console
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s"
            )
            console_handler.setFormatter(formatter)
            logger.addHandler(console_handler)

        return logger

    def handle_error(
        self, error: Exception, context: Optional[Dict[str, Any]] = None
    ) -> None:
        """Handle an error with logging and context."""
        error_type = type(error).__name__

        # Track error counts
        self.error_counts[error_type] = self.error_counts.get(error_type, 0) + 1

        # Log the error
        if isinstance(error, GVClassError):
            error_msg = str(error)
            if error.context:
                error_msg += f" | Context: {error.context}"
        else:
            error_msg = f"{error_type}: {str(error)}"

        if context:
            error_msg += f" | Additional context: {context}"

        self.logger.error(error_msg)

        # Log traceback for debugging
        if self.logger.level <= logging.DEBUG:
            self.logger.debug(traceback.format_exc())

        # Check if we've exceeded error limits
        if self.error_counts[error_type] >= self.max_errors_per_type:
            raise ProcessingError(
                f"Too many errors of type {error_type} ({self.error_counts[error_type]}). "
                f"Stopping processing to prevent cascading failures."
            )

    def log_info(self, message: str, context: Optional[Dict[str, Any]] = None) -> None:
        """Log informational message."""
        msg = message
        if context:
            msg += f" | Context: {context}"
        self.logger.info(msg)

    def log_warning(
        self, message: str, context: Optional[Dict[str, Any]] = None
    ) -> None:
        """Log warning message."""
        msg = message
        if context:
            msg += f" | Context: {context}"
        self.logger.warning(msg)

    def log_debug(self, message: str, context: Optional[Dict[str, Any]] = None) -> None:
        """Log debug message."""
        msg = message
        if context:
            msg += f" | Context: {context}"
        self.logger.debug(msg)


# Global error handler instance
error_handler = ErrorHandler()


def handle_errors(
    error_types: Optional[Union[type, tuple]] = None,
    reraise: bool = True,
    context: Optional[Dict[str, Any]] = None,
) -> Callable[[F], F]:
    """
    Decorator to handle errors in functions.

    Args:
        error_types: Exception types to handle (default: all)
        reraise: Whether to reraise the exception after handling
        context: Additional context to include in error logging

    Returns:
        Decorated function
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                if error_types is None or isinstance(e, error_types):
                    error_handler.handle_error(e, context)
                    if reraise:
                        raise
                else:
                    raise

        return wrapper

    return decorator


@contextmanager
def error_context(step: str, input_file: Optional[str] = None):
    """
    Context manager for error handling in processing steps.

    Args:
        step: Name of the processing step
        input_file: Input file being processed
    """
    context = {"step": step}
    if input_file:
        context["input_file"] = input_file

    try:
        error_handler.log_info(f"Starting step: {step}", context)
        yield
        error_handler.log_info(f"Completed step: {step}", context)
    except Exception as e:
        error_handler.handle_error(e, context)
        raise


def validate_input(
    value: Any, field: str, validator: Callable[[Any], bool], error_msg: str
) -> Any:
    """
    Validate input value with custom validator.

    Args:
        value: Value to validate
        field: Field name for error reporting
        validator: Validation function
        error_msg: Error message if validation fails

    Returns:
        Validated value

    Raises:
        ValidationError: If validation fails
    """
    if not validator(value):
        raise ValidationError(error_msg, field=field, value=value)
    return value


def safe_file_operation(operation: str, file_path: Union[str, Path]) -> Callable:
    """
    Decorator for safe file operations.

    Args:
        operation: Name of the file operation
        file_path: Path to the file

    Returns:
        Decorator function
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except IOError as e:
                raise FileError(
                    f"File operation '{operation}' failed: {str(e)}",
                    file_path=str(file_path),
                    operation=operation,
                )
            except Exception as e:
                if "Permission denied" in str(e):
                    raise FileError(
                        f"Permission denied for file operation '{operation}'",
                        file_path=str(file_path),
                        operation=operation,
                    )
                raise

        return wrapper

    return decorator


def require_files(*file_paths: Union[str, Path]) -> Callable:
    """
    Decorator to ensure required files exist before function execution.

    Args:
        file_paths: Paths to required files

    Returns:
        Decorator function
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            for file_path in file_paths:
                path = Path(file_path)
                if not path.exists():
                    raise FileError(
                        f"Required file not found: {file_path}",
                        file_path=str(file_path),
                        operation="require_files",
                    )
                if not path.is_file():
                    raise FileError(
                        f"Path is not a file: {file_path}",
                        file_path=str(file_path),
                        operation="require_files",
                    )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def retry_on_failure(max_retries: int = 3, delay: float = 1.0) -> Callable:
    """
    Decorator to retry function on failure.

    Args:
        max_retries: Maximum number of retry attempts
        delay: Delay between retries in seconds

    Returns:
        Decorator function
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            import time

            last_exception = None
            for attempt in range(max_retries + 1):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    last_exception = e
                    if attempt < max_retries:
                        error_handler.log_warning(
                            f"Attempt {attempt + 1} failed, retrying in {delay}s",
                            context={"function": func.__name__, "error": str(e)},
                        )
                        time.sleep(delay)
                    else:
                        error_handler.log_error(
                            f"All {max_retries + 1} attempts failed",
                            context={"function": func.__name__},
                        )

            if last_exception:
                raise last_exception

        return wrapper

    return decorator
