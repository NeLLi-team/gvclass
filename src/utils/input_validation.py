"""
Input validation module for GVClass workflow scripts.

This module provides comprehensive input validation for files, parameters,
and other user inputs to ensure data integrity and security.
"""

import os
import re
from pathlib import Path
from typing import Union, Optional, Dict, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .error_handling import ValidationError, FileError, error_handler


class InputValidator:
    """Comprehensive input validation for GVClass pipeline."""
    
    # Constants for validation
    VALID_EXTENSIONS = {'.fna', '.faa', '.fas', '.fasta'}
    MIN_SEQUENCE_LENGTH = 20000  # 20kb minimum
    RECOMMENDED_LENGTH = 50000   # 50kb recommended
    MAX_SEQUENCE_LENGTH = 50000000  # 50MB maximum
    
    # Filename validation patterns
    INVALID_FILENAME_CHARS = re.compile(r'[<>:"|?*\x00-\x1f]')
    INVALID_FILENAME_PATTERNS = re.compile(r'^\.|^CON$|^PRN$|^AUX$|^NUL$|^COM[1-9]$|^LPT[1-9]$', re.IGNORECASE)
    
    # Sequence validation patterns
    VALID_DNA_CHARS = set('ATCGRYSWKMBDHVN-')
    VALID_PROTEIN_CHARS = set('ACDEFGHIKLMNPQRSTVWYXZBJOU*-')
    
    # Parameter validation ranges
    THREAD_LIMITS = (1, 64)
    PROCESS_LIMITS = (1, 128)
    EVALUE_LIMITS = (1e-20, 1.0)
    COVERAGE_LIMITS = (0.0, 100.0)
    
    @classmethod
    def validate_filename(cls, filename: str) -> str:
        """
        Validate and sanitize filename.
        
        Args:
            filename: Filename to validate
            
        Returns:
            Sanitized filename
            
        Raises:
            ValidationError: If filename is invalid
        """
        if not filename:
            raise ValidationError("Filename cannot be empty", field="filename")
        
        # Remove directory components for security
        filename = os.path.basename(filename)
        
        # Check for invalid characters
        if cls.INVALID_FILENAME_CHARS.search(filename):
            raise ValidationError(
                f"Filename contains invalid characters: {filename}",
                field="filename",
                value=filename
            )
        
        # Check for reserved names
        if cls.INVALID_FILENAME_PATTERNS.match(filename):
            raise ValidationError(
                f"Filename uses reserved name: {filename}",
                field="filename",
                value=filename
            )
        
        # Check for path traversal
        if '..' in filename or filename.startswith('/'):
            raise ValidationError(
                f"Filename contains path traversal: {filename}",
                field="filename",
                value=filename
            )
        
        # Check length
        if len(filename) > 255:
            raise ValidationError(
                f"Filename too long (max 255 characters): {len(filename)}",
                field="filename",
                value=filename
            )
        
        return filename
    
    @classmethod
    def validate_sequence_file(cls, filepath: Union[str, Path], file_type: Optional[str] = None) -> Path:
        """
        Validate sequence file format and content.
        
        Args:
            filepath: Path to sequence file
            file_type: Expected file type ('fna' or 'faa')
            
        Returns:
            Validated Path object
            
        Raises:
            ValidationError: If file is invalid
            FileError: If file doesn't exist or can't be read
        """
        if not filepath:
            raise ValidationError("File path cannot be empty", field="filepath")
        
        path = Path(filepath).resolve()
        
        # Check if file exists
        if not path.exists():
            raise FileError(
                f"File not found: {filepath}",
                file_path=str(filepath),
                operation="validate_sequence_file"
            )
        
        # Check if it's a file (not directory)
        if not path.is_file():
            raise FileError(
                f"Path is not a file: {filepath}",
                file_path=str(filepath),
                operation="validate_sequence_file"
            )
        
        # Check file extension
        if path.suffix.lower() not in cls.VALID_EXTENSIONS:
            raise ValidationError(
                f"Invalid file extension: {path.suffix}. Must be one of: {cls.VALID_EXTENSIONS}",
                field="file_extension",
                value=path.suffix
            )
        
        # Check file size
        file_size = path.stat().st_size
        if file_size == 0:
            raise ValidationError(
                f"File is empty: {filepath}",
                field="file_size",
                value=file_size
            )
        
        if file_size > cls.MAX_SEQUENCE_LENGTH:
            raise ValidationError(
                f"File too large: {file_size} bytes (max: {cls.MAX_SEQUENCE_LENGTH})",
                field="file_size",
                value=file_size
            )
        
        # Validate file content
        try:
            sequences = list(SeqIO.parse(str(path), "fasta"))
        except Exception as e:
            raise ValidationError(
                f"Error parsing FASTA file: {str(e)}",
                field="file_content",
                value=str(filepath)
            )
        
        if not sequences:
            raise ValidationError(
                f"No sequences found in file: {filepath}",
                field="sequence_count",
                value=0
            )
        
        # Determine file type if not provided
        if file_type is None:
            file_type = path.suffix.lower().replace('.', '')
            if file_type == 'fas' or file_type == 'fasta':
                # Need to guess based on content
                file_type = None
        
        # Validate individual sequences
        for i, record in enumerate(sequences):
            cls._validate_sequence_record(record, i + 1, file_type or path.suffix.lower().replace('.', ''))
        
        # Check total sequence length (only for nucleotide sequences)
        if path.suffix.lower() == '.fna':
            total_length = sum(len(seq.seq) for seq in sequences)
            if total_length < cls.MIN_SEQUENCE_LENGTH:
                error_handler.log_warning(
                    f"Total sequence length ({total_length} bp) is below recommended minimum ({cls.MIN_SEQUENCE_LENGTH} bp)",
                    context={"file": str(filepath), "total_length": total_length}
                )
        
        return path
    
    @classmethod
    def _validate_sequence_record(cls, record: SeqRecord, seq_num: int, file_type: Optional[str] = None) -> None:
        """
        Validate individual sequence record.
        
        Args:
            record: SeqRecord to validate
            seq_num: Sequence number for error reporting
            file_type: Expected file type
            
        Raises:
            ValidationError: If sequence is invalid
        """
        # Validate sequence ID
        if not record.id:
            raise ValidationError(
                f"Sequence {seq_num} has no ID",
                field="sequence_id",
                value=seq_num
            )
        
        # Check for invalid characters in sequence ID (allow pipes for formatting)
        # Note: We allow pipes (|) as they're commonly used in sequence IDs
        invalid_chars_pattern = re.compile(r'[<>:"?*\x00-\x1f]')
        if invalid_chars_pattern.search(record.id):
            raise ValidationError(
                f"Sequence {seq_num} ID contains invalid characters: {record.id}",
                field="sequence_id",
                value=record.id
            )
        
        # Validate sequence length
        seq_length = len(record.seq)
        if seq_length == 0:
            raise ValidationError(
                f"Sequence {seq_num} ({record.id}) is empty",
                field="sequence_length",
                value=seq_length
            )
        
        # Only check minimum length for nucleotide sequences
        if file_type == 'fna':
            if seq_length < cls.MIN_SEQUENCE_LENGTH:
                error_handler.log_warning(
                    f"Sequence {seq_num} ({record.id}) is shorter than recommended minimum ({seq_length} bp < {cls.MIN_SEQUENCE_LENGTH} bp)",
                    context={"sequence_id": record.id, "length": seq_length}
                )
        elif file_type == 'faa':
            # For proteins, check a different minimum (e.g., 100 amino acids ~ 300 bp equivalent)
            min_protein_length = 100
            if seq_length < min_protein_length:
                error_handler.log_warning(
                    f"Protein sequence {seq_num} ({record.id}) is shorter than recommended minimum ({seq_length} aa < {min_protein_length} aa)",
                    context={"sequence_id": record.id, "length": seq_length}
                )
        
        # Validate sequence content based on type
        seq_str = str(record.seq).upper()
        
        if file_type == 'fna' or record.id.endswith('.fna'):
            # DNA sequence validation
            invalid_chars = set(seq_str) - cls.VALID_DNA_CHARS
            if invalid_chars:
                raise ValidationError(
                    f"Sequence {seq_num} ({record.id}) contains invalid DNA characters: {invalid_chars}",
                    field="sequence_content",
                    value=sorted(invalid_chars)
                )
        
        elif file_type == 'faa' or record.id.endswith('.faa'):
            # Protein sequence validation
            invalid_chars = set(seq_str) - cls.VALID_PROTEIN_CHARS
            if invalid_chars:
                raise ValidationError(
                    f"Sequence {seq_num} ({record.id}) contains invalid protein characters: {invalid_chars}",
                    field="sequence_content",
                    value=sorted(invalid_chars)
                )
    
    @classmethod
    def validate_directory(cls, dirpath: Union[str, Path], must_exist: bool = True, must_be_writable: bool = False) -> Path:
        """
        Validate directory path.
        
        Args:
            dirpath: Directory path to validate
            must_exist: Whether directory must exist
            must_be_writable: Whether directory must be writable
            
        Returns:
            Validated Path object
            
        Raises:
            ValidationError: If directory is invalid
            FileError: If directory doesn't exist or isn't accessible
        """
        if not dirpath:
            raise ValidationError("Directory path cannot be empty", field="dirpath")
        
        path = Path(dirpath).resolve()
        
        # Check for path traversal
        if '..' in str(dirpath):
            raise ValidationError(
                f"Path traversal detected in directory: {dirpath}",
                field="dirpath",
                value=dirpath
            )
        
        if must_exist:
            if not path.exists():
                raise FileError(
                    f"Directory not found: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory"
                )
            
            if not path.is_dir():
                raise FileError(
                    f"Path is not a directory: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory"
                )
        
        if must_be_writable and path.exists():
            if not os.access(path, os.W_OK):
                raise FileError(
                    f"Directory is not writable: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory"
                )
        
        return path
    
    @classmethod
    def validate_threads(cls, threads: int) -> int:
        """
        Validate thread count parameter.
        
        Args:
            threads: Number of threads
            
        Returns:
            Validated thread count
            
        Raises:
            ValidationError: If thread count is invalid
        """
        if not isinstance(threads, int):
            raise ValidationError(
                f"Thread count must be an integer, got {type(threads).__name__}",
                field="threads",
                value=threads
            )
        
        if not (cls.THREAD_LIMITS[0] <= threads <= cls.THREAD_LIMITS[1]):
            raise ValidationError(
                f"Thread count must be between {cls.THREAD_LIMITS[0]} and {cls.THREAD_LIMITS[1]}, got {threads}",
                field="threads",
                value=threads
            )
        
        return threads
    
    @classmethod
    def validate_processes(cls, processes: int) -> int:
        """
        Validate process count parameter.
        
        Args:
            processes: Number of processes
            
        Returns:
            Validated process count
            
        Raises:
            ValidationError: If process count is invalid
        """
        if not isinstance(processes, int):
            raise ValidationError(
                f"Process count must be an integer, got {type(processes).__name__}",
                field="processes",
                value=processes
            )
        
        if not (cls.PROCESS_LIMITS[0] <= processes <= cls.PROCESS_LIMITS[1]):
            raise ValidationError(
                f"Process count must be between {cls.PROCESS_LIMITS[0]} and {cls.PROCESS_LIMITS[1]}, got {processes}",
                field="processes",
                value=processes
            )
        
        return processes
    
    @classmethod
    def validate_evalue(cls, evalue: float) -> float:
        """
        Validate E-value parameter.
        
        Args:
            evalue: E-value threshold
            
        Returns:
            Validated E-value
            
        Raises:
            ValidationError: If E-value is invalid
        """
        if not isinstance(evalue, (int, float)):
            raise ValidationError(
                f"E-value must be a number, got {type(evalue).__name__}",
                field="evalue",
                value=evalue
            )
        
        if not (cls.EVALUE_LIMITS[0] <= evalue <= cls.EVALUE_LIMITS[1]):
            raise ValidationError(
                f"E-value must be between {cls.EVALUE_LIMITS[0]} and {cls.EVALUE_LIMITS[1]}, got {evalue}",
                field="evalue",
                value=evalue
            )
        
        return float(evalue)
    
    @classmethod
    def validate_coverage(cls, coverage: float) -> float:
        """
        Validate coverage percentage parameter.
        
        Args:
            coverage: Coverage percentage
            
        Returns:
            Validated coverage percentage
            
        Raises:
            ValidationError: If coverage is invalid
        """
        if not isinstance(coverage, (int, float)):
            raise ValidationError(
                f"Coverage must be a number, got {type(coverage).__name__}",
                field="coverage",
                value=coverage
            )
        
        if not (cls.COVERAGE_LIMITS[0] <= coverage <= cls.COVERAGE_LIMITS[1]):
            raise ValidationError(
                f"Coverage must be between {cls.COVERAGE_LIMITS[0]}% and {cls.COVERAGE_LIMITS[1]}%, got {coverage}%",
                field="coverage",
                value=coverage
            )
        
        return float(coverage)
    
    @classmethod
    def validate_genetic_code(cls, code: int) -> int:
        """
        Validate genetic code parameter.
        
        Args:
            code: Genetic code number
            
        Returns:
            Validated genetic code
            
        Raises:
            ValidationError: If genetic code is invalid
        """
        # Valid genetic codes from NCBI
        valid_codes = {0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 106, 129}
        
        if not isinstance(code, int):
            raise ValidationError(
                f"Genetic code must be an integer, got {type(code).__name__}",
                field="genetic_code",
                value=code
            )
        
        if code not in valid_codes:
            raise ValidationError(
                f"Invalid genetic code: {code}. Must be one of: {sorted(valid_codes)}",
                field="genetic_code",
                value=code
            )
        
        return code
    
    @classmethod
    def validate_mafft_option(cls, option: str) -> str:
        """
        Validate MAFFT alignment option.
        
        Args:
            option: MAFFT option
            
        Returns:
            Validated option
            
        Raises:
            ValidationError: If option is invalid
        """
        valid_options = {"auto", "linsi", "ginsi", "einsi", "fftns", "fftnsi", "nwns", "nwnsi"}
        
        if not isinstance(option, str):
            raise ValidationError(
                f"MAFFT option must be a string, got {type(option).__name__}",
                field="mafft_option",
                value=option
            )
        
        if option not in valid_options:
            raise ValidationError(
                f"Invalid MAFFT option: {option}. Must be one of: {sorted(valid_options)}",
                field="mafft_option",
                value=option
            )
        
        return option
    
    @classmethod
    def validate_tree_method(cls, method: str) -> str:
        """
        Validate tree building method.
        
        Args:
            method: Tree building method
            
        Returns:
            Validated method
            
        Raises:
            ValidationError: If method is invalid
        """
        valid_methods = {"iqtree", "fasttree"}
        
        if not isinstance(method, str):
            raise ValidationError(
                f"Tree method must be a string, got {type(method).__name__}",
                field="tree_method",
                value=method
            )
        
        if method not in valid_methods:
            raise ValidationError(
                f"Invalid tree method: {method}. Must be one of: {sorted(valid_methods)}",
                field="tree_method",
                value=method
            )
        
        return method
    
    @classmethod
    def validate_query_directory(cls, query_dir: Union[str, Path]) -> Path:
        """
        Validate query directory and its contents.
        
        Args:
            query_dir: Path to query directory
            
        Returns:
            Validated Path object
            
        Raises:
            ValidationError: If directory or contents are invalid
        """
        dir_path = cls.validate_directory(query_dir, must_exist=True)
        
        # Check for valid sequence files
        valid_files = []
        for ext in cls.VALID_EXTENSIONS:
            valid_files.extend(dir_path.glob(f"*{ext}"))
        
        if not valid_files:
            raise ValidationError(
                f"No valid sequence files found in query directory: {query_dir}",
                field="query_directory",
                value=str(query_dir)
            )
        
        # Validate each file
        for file_path in valid_files:
            try:
                cls.validate_sequence_file(file_path)
            except ValidationError as e:
                error_handler.log_warning(
                    f"Invalid sequence file found: {file_path}",
                    context={"file": str(file_path), "error": str(e)}
                )
        
        return dir_path
    
    @classmethod
    def validate_config_parameters(cls, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate configuration parameters.
        
        Args:
            config: Configuration dictionary
            
        Returns:
            Validated configuration dictionary
            
        Raises:
            ValidationError: If any parameter is invalid
        """
        validated_config = {}
        
        # Validate known parameters
        if 'threads' in config:
            validated_config['threads'] = cls.validate_threads(config['threads'])
        
        if 'processes' in config:
            validated_config['processes'] = cls.validate_processes(config['processes'])
        
        if 'evalue' in config:
            validated_config['evalue'] = cls.validate_evalue(config['evalue'])
        
        if 'query_coverage' in config:
            validated_config['query_coverage'] = cls.validate_coverage(config['query_coverage'])
        
        if 'subject_coverage' in config:
            validated_config['subject_coverage'] = cls.validate_coverage(config['subject_coverage'])
        
        if 'mafft_option' in config:
            validated_config['mafft_option'] = cls.validate_mafft_option(config['mafft_option'])
        
        if 'tree_method' in config:
            validated_config['tree_method'] = cls.validate_tree_method(config['tree_method'])
        
        if 'genetic_codes' in config:
            if isinstance(config['genetic_codes'], list):
                validated_config['genetic_codes'] = [cls.validate_genetic_code(code) for code in config['genetic_codes']]
            else:
                raise ValidationError(
                    "genetic_codes must be a list",
                    field="genetic_codes",
                    value=config['genetic_codes']
                )
        
        # Copy other parameters without validation (with warning)
        for key, value in config.items():
            if key not in validated_config:
                validated_config[key] = value
                error_handler.log_warning(
                    f"Unknown configuration parameter: {key}",
                    context={"parameter": key, "value": value}
                )
        
        return validated_config