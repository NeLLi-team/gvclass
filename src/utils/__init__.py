"""
Utility modules for GVClass pipeline.
"""

from .error_handling import ErrorHandler, GVClassError, ValidationError, ProcessingError, error_handler
from .input_validation import InputValidator
from .common import (
    run_command_safe,
    run_command_safe_legacy,
    validate_file_path,
    setup_logging
)
from .database_manager import DatabaseManager

# Apply pyrodigal patches on module import for genetic codes 106 and 129
try:
    from .pyrodigal_patcher import patch_pyrodigal_genetic_codes
    patch_pyrodigal_genetic_codes()
except Exception:
    pass  # Fail silently if pyrodigal not installed

__all__ = [
    # Error handling
    'ErrorHandler',
    'error_handler',
    'GVClassError', 
    'ValidationError',
    'ProcessingError',
    # Input validation
    'InputValidator',
    # Common utilities
    'run_command_safe',
    'run_command_safe_legacy',
    'validate_file_path',
    'setup_logging',
    # Database management
    'DatabaseManager'
]