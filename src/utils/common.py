"""
Utility functions for GVClass workflow scripts.

This module provides common utilities for secure subprocess execution,
logging, and error handling.
"""

import logging
import subprocess
import shlex
from typing import List, Tuple, Optional, Union
from pathlib import Path

# Import error handling framework
from .error_handling import (
    error_handler, 
    CommandError, 
    ValidationError, 
    FileError,
    handle_errors
)


def setup_logging(name: str, level: str = "INFO") -> logging.Logger:
    """
    Set up logging with consistent formatting.
    
    Args:
        name: Logger name
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        
    Returns:
        Configured logger instance
    """
    return setup_logger(name, level)


def setup_logger(name: str, level: str = "INFO") -> logging.Logger:
    """
    Set up a logger with consistent formatting.
    
    Args:
        name: Logger name
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


@handle_errors(error_types=(subprocess.CalledProcessError, subprocess.TimeoutExpired))
def run_command_safe(
    cmd: Union[str, List[str]], 
    cwd: Optional[Union[str, Path]] = None,
    timeout: Optional[int] = None
) -> Tuple[str, str]:
    """
    Execute a command safely without shell injection vulnerability.
    
    Args:
        cmd: Command to execute (string will be split safely)
        cwd: Working directory for command execution
        timeout: Command timeout in seconds
        
    Returns:
        Tuple of (stdout, stderr)
        
    Raises:
        CommandError: If command fails
        subprocess.TimeoutExpired: If command times out
    """
    # Convert string command to list safely
    if isinstance(cmd, str):
        cmd_list = shlex.split(cmd)
    else:
        cmd_list = cmd
    
    # Validate command is not empty
    if not cmd_list:
        raise CommandError("Command cannot be empty")
    
    command_str = ' '.join(cmd_list)
    error_handler.log_debug(f"Executing command: {command_str}")
    
    try:
        result = subprocess.run(
            cmd_list,
            cwd=cwd,
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout
        )
        
        if result.stderr:
            error_handler.log_debug(f"Command stderr: {result.stderr}")
        if result.stdout:
            error_handler.log_debug(f"Command stdout: {result.stdout}")
            
        return result.stdout, result.stderr
        
    except subprocess.CalledProcessError as e:
        raise CommandError(
            f"Command failed: {str(e)}",
            command=command_str,
            return_code=e.returncode
        ) from e
    except subprocess.TimeoutExpired as e:
        raise CommandError(
            f"Command timed out after {timeout} seconds",
            command=command_str
        ) from e


def run_command_safe_legacy(cmd: Union[str, List[str]]) -> None:
    """
    Legacy wrapper for run_command_safe that prints output.
    
    This function maintains compatibility with existing code that expects
    output to be printed rather than returned.
    
    Args:
        cmd: Command to execute
    """
    try:
        stdout, stderr = run_command_safe(cmd)
        if stderr:
            print(f'std_err: {stderr}')
        if stdout:
            print(f'std_out: {stdout}')
    except CommandError as e:
        print(f'Command failed: {e}')
        raise


@handle_errors(error_types=(ValueError, FileNotFoundError))
def validate_file_path(file_path: Union[str, Path], must_exist: bool = True) -> Path:
    """
    Validate and normalize a file path.
    
    Args:
        file_path: Path to validate
        must_exist: Whether the file must exist
        
    Returns:
        Normalized Path object
        
    Raises:
        ValidationError: If path is invalid
        FileError: If file doesn't exist and must_exist is True
    """
    if not file_path:
        raise ValidationError("File path cannot be empty", field="file_path")
    
    path = Path(file_path).resolve()
    
    # Check for path traversal attempts
    if '..' in str(file_path):
        raise ValidationError(
            f"Path traversal detected in: {file_path}",
            field="file_path",
            value=file_path
        )
    
    if must_exist and not path.exists():
        raise FileError(
            f"File not found: {file_path}",
            file_path=str(file_path),
            operation="validate_path"
        )
    
    return path


# Legacy compatibility
class GVClassError(Exception):
    """Base exception for GVClass errors."""
    pass