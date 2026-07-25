"""
Input validation module for GVClass workflow scripts.

This module provides comprehensive input validation for files, parameters,
and other user inputs to ensure data integrity and security.
"""

import os
import re
from pathlib import Path
from typing import Union, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .error_handling import ValidationError, FileError, error_handler


class InputValidator:
    """Comprehensive input validation for GVClass pipeline."""

    # Constants for validation
    VALID_EXTENSIONS = {".fna", ".faa", ".fas", ".fasta"}
    MIN_SEQUENCE_LENGTH = 20000  # 20kb minimum
    # Maximum input file size, in bytes. Was 50 MB under a constant named for
    # sequence length, which rejected the assembly FASTA files --contigs exists
    # to process. The ceiling is bounded by memory, not by policy: validation
    # and gene calling each hold a whole file via list(SeqIO.parse(...)), at
    # roughly 2x the file size in RSS, and the per-query stages run
    # concurrently. 500 MB keeps one parse near 1 GB. Raising it further needs
    # streaming parsers first.
    MAX_FILE_SIZE_BYTES = 500_000_000

    # Sequence validation patterns
    VALID_DNA_CHARS = set("ATCGRYSWKMBDHVN-")
    VALID_PROTEIN_CHARS = set("ACDEFGHIKLMNPQRSTVWYXZBJOU*-")

    # Full IUPAC nucleotide alphabet. Any character outside this set
    # (E, F, I, L, P, Q, ...) can only be an amino acid. Deliberately matches
    # VALID_DNA_CHARS: inferring a type this set admits but the alphabet
    # validator rejects would turn a working input into a hard error. U is
    # therefore absent here too, so RNA keeps its pre-2.0.2 handling.
    NUCLEOTIDE_ALPHABET = VALID_DNA_CHARS
    # The four unambiguous bases. Real nucleotide sequence is dominated by
    # these; protein is not, because A/C/G/T are only 4 of 20 residues.
    UNAMBIGUOUS_BASES = set("ACGT")
    # Characters carrying no alphabet signal either way.
    UNINFORMATIVE_CHARS = set("N-")
    # Below this many informative characters, base diversity is not evidence.
    MIN_LENGTH_FOR_BASE_DIVERSITY = 12

    @classmethod
    def _infer_fasta_type_from_content(cls, sequences) -> str:
        """Return ``'fna'`` or ``'faa'`` from the first informative sequence.

        Used for extension-agnostic ``.fasta``/``.fas`` files, both to pick the
        alphabet for :meth:`_validate_sequence_record` and to route the query to
        the nucleotide or protein branch in
        :func:`~src.pipeline.query_processing_engine._is_nucleotide_input`.

        Three steps, because no one of them is sufficient alone:

        1. Any character outside :data:`NUCLEOTIDE_ALPHABET` means protein.
           Ordinary proteins nearly always contain E, F, I, L, P or Q.
        2. A sequence long enough to judge that uses fewer than three of the
           four bases is protein. Real nucleotide sequence uses all four; a
           low-complexity protein repeat such as the silk-fibroin ``GAGAGS``
           motif uses two and would otherwise pass step 3, since ``S`` is also
           an IUPAC ambiguity code.
        3. Otherwise nucleotide when the unambiguous bases are at least half of
           the informative characters (ignoring ``N`` and gaps). This is what
           separates the peptide ``MKWVNDHRYSGACD``, whose every letter is also
           an ambiguity code, from a contig carrying many ``R``/``Y``/``S``/``W``
           positions.

        Only the first non-empty record is examined, so it decides the whole
        file. Defaults to ``'fna'`` when every sequence is empty.
        """
        for record in sequences:
            seq = str(record.seq).upper()
            if not seq:
                continue
            if not set(seq) <= cls.NUCLEOTIDE_ALPHABET:
                return "faa"
            informative = [c for c in seq if c not in cls.UNINFORMATIVE_CHARS]
            if not informative:
                return "fna"
            distinct_bases = {c for c in informative if c in cls.UNAMBIGUOUS_BASES}
            if len(informative) >= cls.MIN_LENGTH_FOR_BASE_DIVERSITY:
                if len(distinct_bases - {"U"}) < 3 and len(distinct_bases - {"T"}) < 3:
                    return "faa"
            unambiguous = sum(1 for c in informative if c in cls.UNAMBIGUOUS_BASES)
            return "fna" if unambiguous / len(informative) >= 0.5 else "faa"
        return "fna"

    @classmethod
    def validate_sequence_file(
        cls,
        filepath: Union[str, Path],
        file_type: Optional[str] = None,
        allow_short: bool = False,
        min_sequence_length: Optional[int] = None,
    ) -> Path:
        """
        Validate sequence file format and content.

        Args:
            filepath: Path to sequence file
            file_type: Expected file type ('fna' or 'faa')
            allow_short: When False (default) aggregate nucleotide files
                shorter than ``MIN_SEQUENCE_LENGTH`` raise ``ValidationError``
                instead of emitting a warning.
            min_sequence_length: Minimum aggregate nucleotide length in bp.
                Defaults to ``MIN_SEQUENCE_LENGTH``. ``--contigs`` mode passes
                its split filter here so per-contig inputs are governed by
                ``pipeline.contigs_min_length``.

        Returns:
            Validated Path object

        Raises:
            ValidationError: If file is invalid
            FileError: If file doesn't exist or can't be read
        """
        if not filepath:
            raise ValidationError("File path cannot be empty", field="filepath")

        if min_sequence_length is None:
            min_sequence_length = cls.MIN_SEQUENCE_LENGTH
        try:
            min_sequence_length = int(min_sequence_length)
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                f"Minimum sequence length must be an integer: {min_sequence_length!r}",
                field="min_sequence_length",
                value=min_sequence_length,
            ) from exc
        if min_sequence_length < 0:
            raise ValidationError(
                f"Minimum sequence length must be >= 0: {min_sequence_length}",
                field="min_sequence_length",
                value=min_sequence_length,
            )

        path = Path(filepath).resolve()

        # Check if file exists
        if not path.exists():
            raise FileError(
                f"File not found: {filepath}",
                file_path=str(filepath),
                operation="validate_sequence_file",
            )

        # Check if it's a file (not directory)
        if not path.is_file():
            raise FileError(
                f"Path is not a file: {filepath}",
                file_path=str(filepath),
                operation="validate_sequence_file",
            )

        # Check file extension
        if path.suffix.lower() not in cls.VALID_EXTENSIONS:
            raise ValidationError(
                f"Invalid file extension: {path.suffix}. Must be one of: {cls.VALID_EXTENSIONS}",
                field="file_extension",
                value=path.suffix,
            )

        # Check file size
        file_size = path.stat().st_size
        if file_size == 0:
            raise ValidationError(
                f"File is empty: {filepath}", field="file_size", value=file_size
            )

        if file_size > cls.MAX_FILE_SIZE_BYTES:
            raise ValidationError(
                f"File too large: {file_size} bytes exceeds the "
                f"{cls.MAX_FILE_SIZE_BYTES // 1_000_000} MB limit. "
                "Split the file or pass fewer sequences per file.",
                field="file_size",
                value=file_size,
            )

        # Validate file content
        try:
            sequences = list(SeqIO.parse(str(path), "fasta"))
        except Exception as e:
            raise ValidationError(
                f"Error parsing FASTA file: {str(e)}",
                field="file_content",
                value=str(filepath),
            )

        if not sequences:
            raise ValidationError(
                f"No sequences found in file: {filepath}",
                field="sequence_count",
                value=0,
            )

        # Determine file type if not provided. ``.fasta``/``.fas`` files
        # previously fell through to an extension-based branch that skipped
        # DNA/protein validation entirely (Codex-audit finding). Now we
        # always resolve to a concrete ``fna``/``faa`` type by peeking at
        # the first informative sequence.
        if file_type is None:
            suffix_type = path.suffix.lower().replace(".", "")
            if suffix_type in ("fas", "fasta"):
                file_type = cls._infer_fasta_type_from_content(sequences)
            else:
                file_type = suffix_type

        # Validate individual sequences
        for i, record in enumerate(sequences):
            cls._validate_sequence_record(
                record,
                i + 1,
                file_type,
                allow_short=allow_short,
                min_sequence_length=min_sequence_length,
            )

        # Check total sequence length for nucleotide inputs. Hard-fail unless
        # the caller opted into ``allow_short``. The floor is caller-provided
        # so ``--contigs`` can use ``pipeline.contigs_min_length``.
        if file_type == "fna":
            total_length = sum(len(seq.seq) for seq in sequences)
            if total_length < min_sequence_length:
                message = (
                    f"Total sequence length ({total_length} bp) is below the "
                    f"minimum ({min_sequence_length} bp) for "
                    f"{filepath}. Pass --allow-short to override this gate."
                )
                if allow_short:
                    error_handler.log_warning(
                        message,
                        context={"file": str(filepath), "total_length": total_length},
                    )
                else:
                    raise ValidationError(
                        message,
                        field="sequence_length",
                        value=total_length,
                    )

        return path

    @classmethod
    def _validate_sequence_record(
        cls,
        record: SeqRecord,
        seq_num: int,
        file_type: Optional[str] = None,
        allow_short: bool = False,
        min_sequence_length: Optional[int] = None,
    ) -> None:
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
                f"Sequence {seq_num} has no ID", field="sequence_id", value=seq_num
            )

        # Check for invalid characters in sequence ID (allow pipes for formatting)
        # Note: We allow pipes (|) as they're commonly used in sequence IDs
        invalid_chars_pattern = re.compile(r'[<>:"?*\x00-\x1f]')
        if invalid_chars_pattern.search(record.id):
            raise ValidationError(
                f"Sequence {seq_num} ID contains invalid characters: {record.id}",
                field="sequence_id",
                value=record.id,
            )

        # Validate sequence length
        seq_length = len(record.seq)
        if seq_length == 0:
            raise ValidationError(
                f"Sequence {seq_num} ({record.id}) is empty",
                field="sequence_length",
                value=seq_length,
            )

        # Only check minimum length for nucleotide sequences. Per-record
        # length is still a soft warning (a single short contig in a
        # multi-contig FNA is normal); the aggregate hard-fail lives in
        # ``validate_sequence_file`` where it can also honour allow_short.
        if file_type == "fna":
            length_floor = (
                cls.MIN_SEQUENCE_LENGTH
                if min_sequence_length is None
                else int(min_sequence_length)
            )
            if seq_length < length_floor:
                error_handler.log_warning(
                    f"Sequence {seq_num} ({record.id}) is shorter than recommended minimum ({seq_length} bp < {length_floor} bp)",
                    context={"sequence_id": record.id, "length": seq_length},
                )
        elif file_type == "faa":
            # For proteins, check a different minimum (e.g., 100 amino acids ~ 300 bp equivalent)
            min_protein_length = 100
            if seq_length < min_protein_length:
                error_handler.log_warning(
                    f"Protein sequence {seq_num} ({record.id}) is shorter than recommended minimum ({seq_length} aa < {min_protein_length} aa)",
                    context={"sequence_id": record.id, "length": seq_length},
                )
        _ = allow_short

        # Validate sequence content based on type
        seq_str = str(record.seq).upper()

        if file_type == "fna" or record.id.endswith(".fna"):
            # DNA sequence validation
            invalid_chars = set(seq_str) - cls.VALID_DNA_CHARS
            if invalid_chars:
                raise ValidationError(
                    f"Sequence {seq_num} ({record.id}) contains invalid DNA characters: {invalid_chars}",
                    field="sequence_content",
                    value=sorted(invalid_chars),
                )

        elif file_type == "faa" or record.id.endswith(".faa"):
            # Protein sequence validation
            invalid_chars = set(seq_str) - cls.VALID_PROTEIN_CHARS
            if invalid_chars:
                raise ValidationError(
                    f"Sequence {seq_num} ({record.id}) contains invalid protein characters: {invalid_chars}",
                    field="sequence_content",
                    value=sorted(invalid_chars),
                )
        # ``.fasta``/``.fas`` files resolved to a concrete file_type above;
        # any remaining unknown type here means the caller invoked us with
        # an unrecognized extension, which earlier validation already
        # rejected. Guard against future regressions.

    @classmethod
    def validate_directory(
        cls,
        dirpath: Union[str, Path],
        must_exist: bool = True,
        must_be_writable: bool = False,
    ) -> Path:
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
        if ".." in str(dirpath):
            raise ValidationError(
                f"Path traversal detected in directory: {dirpath}",
                field="dirpath",
                value=dirpath,
            )

        if must_exist:
            if not path.exists():
                raise FileError(
                    f"Directory not found: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory",
                )

            if not path.is_dir():
                raise FileError(
                    f"Path is not a directory: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory",
                )

        if must_be_writable and path.exists():
            if not os.access(path, os.W_OK):
                raise FileError(
                    f"Directory is not writable: {dirpath}",
                    file_path=str(dirpath),
                    operation="validate_directory",
                )

        return path

    @classmethod
    def validate_query_directory(
        cls,
        query_dir: Union[str, Path],
        allow_short: bool = False,
        min_sequence_length: Optional[int] = None,
    ) -> Path:
        """
        Validate query directory and its contents.

        Args:
            query_dir: Path to query directory
            allow_short: Forwarded to :meth:`validate_sequence_file` so
                ``--contigs`` and ``--allow-short`` runs can opt out of the
                configured aggregate nucleotide minimum.
            min_sequence_length: Minimum aggregate nucleotide length in bp.
                Forwarded to :meth:`validate_sequence_file`.

        Returns:
            Validated Path object

        Raises:
            ValidationError: If directory or contents are invalid. The
                previous implementation swallowed per-file
                ``ValidationError`` into ``log_warning`` and silently
                proceeded, which allowed malformed FASTA to reach the
                downstream pipeline. The Codex audit required this to be
                fail-closed.
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
                value=str(query_dir),
            )

        # Validate each file — fail-closed on the first invalid file so
        # bad inputs do not reach the pipeline silently.
        for file_path in sorted(valid_files):
            cls.validate_sequence_file(
                file_path,
                allow_short=allow_short,
                min_sequence_length=min_sequence_length,
            )

        return dir_path
