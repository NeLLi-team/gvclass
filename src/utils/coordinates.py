"""
Coordinate system conversion utilities.

This module provides consistent conversion between 0-based and 1-based coordinate systems.

Convention:
- Internal processing: 0-based (Python standard)
- Output formats: 1-based (bioinformatics standard for BLAST, HMM, etc.)

0-based: First position is 0, ranges are [start, end)
1-based: First position is 1, ranges are [start, end]
"""


def to_1based(start: int, end: int) -> tuple[int, int]:
    """
    Convert 0-based coordinates to 1-based coordinates.

    Args:
        start: 0-based start position (inclusive)
        end: 0-based end position (exclusive in 0-based, becomes inclusive in 1-based)

    Returns:
        Tuple of (1-based start, 1-based end) both inclusive

    Example:
        >>> to_1based(0, 10)  # First 10 positions in 0-based
        (1, 10)  # Positions 1-10 in 1-based
    """
    return start + 1, end


def to_1based_inclusive(start: int, end: int) -> tuple[int, int]:
    """
    Convert 0-based coordinates to 1-based when both systems use inclusive end.

    Args:
        start: 0-based start position (inclusive)
        end: 0-based end position (inclusive)

    Returns:
        Tuple of (1-based start, 1-based end) both inclusive

    Example:
        >>> to_1based_inclusive(0, 9)  # Positions 0-9 inclusive in 0-based
        (1, 10)  # Positions 1-10 inclusive in 1-based
    """
    return start + 1, end + 1


def to_0based(start: int, end: int) -> tuple[int, int]:
    """
    Convert 1-based coordinates to 0-based coordinates.

    Args:
        start: 1-based start position (inclusive)
        end: 1-based end position (inclusive)

    Returns:
        Tuple of (0-based start, 0-based end) where end is exclusive

    Example:
        >>> to_0based(1, 10)  # Positions 1-10 in 1-based
        (0, 10)  # Positions [0, 10) in 0-based
    """
    return start - 1, end


def to_0based_inclusive(start: int, end: int) -> tuple[int, int]:
    """
    Convert 1-based coordinates to 0-based when both systems use inclusive end.

    Args:
        start: 1-based start position (inclusive)
        end: 1-based end position (inclusive)

    Returns:
        Tuple of (0-based start, 0-based end) both inclusive

    Example:
        >>> to_0based_inclusive(1, 10)  # Positions 1-10 in 1-based
        (0, 9)  # Positions 0-9 inclusive in 0-based
    """
    return start - 1, end - 1


def validate_coordinates(start: int, end: int, system: str = "0-based") -> bool:
    """
    Validate that coordinates are sensible.

    Args:
        start: Start position
        end: End position
        system: Either "0-based" or "1-based"

    Returns:
        True if coordinates are valid

    Raises:
        ValueError: If coordinates are invalid
    """
    min_value = 0 if system == "0-based" else 1

    if start < min_value:
        raise ValueError(
            f"Start position {start} cannot be less than {min_value} in {system} system"
        )

    if end < start:
        raise ValueError(
            f"End position {end} cannot be less than start position {start}"
        )

    return True
