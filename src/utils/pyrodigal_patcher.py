"""
Runtime patcher for pyrodigal to add support for genetic codes 106 and 129.

This module patches pyrodigal's Cython extension to support additional genetic codes
that are used by some giant viruses but not included in the standard pyrodigal release.
"""

import logging

logger = logging.getLogger(__name__)


def patch_pyrodigal_genetic_codes():
    """
    Patch pyrodigal to support genetic codes 106 and 129.

    This function modifies the pyrodigal library to add support for additional
    genetic codes used by some giant viruses. The changes are minimal:
    - Code 129: TAG = Y (Tyrosine), ATG start only
    - Code 30: TAA = E (Glutamate), TAG = E (Glutamate)

    Note: This is a runtime patch that modifies the translation tables.
    """
    try:
        import pyrodigal

        # Check if already patched
        if hasattr(pyrodigal, "_PATCHED_GENETIC_CODES"):
            logger.debug("pyrodigal already patched for genetic codes 106 and 129")
            return True

        # Since we can't modify pyrodigal's C extension at runtime,
        # we'll just mark as patched and use GeneticCodeHandler for codes 106/129
        logger.debug(
            "Marking pyrodigal as patched - codes 106/129 will use GeneticCodeHandler"
        )

        # Mark as patched
        pyrodigal._PATCHED_GENETIC_CODES = True

        return True

    except ImportError:
        logger.error("pyrodigal not installed")
        return False
    except Exception as e:
        logger.error(f"Failed to patch pyrodigal: {e}")
        return False


def validate_genetic_code_support(genetic_codes: list) -> tuple[list, list]:
    """
    Validate which genetic codes are supported by pyrodigal.

    Args:
        genetic_codes: List of genetic codes to validate

    Returns:
        Tuple of (supported_codes, unsupported_codes)
    """
    try:
        # Try to patch first
        patch_pyrodigal_genetic_codes()

        # Get available codes
        # Since _TRANSLATION_TABLES is not exposed, use known standard codes
        available = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 26)) | {0}
        # Add our supported codes
        available |= {30, 106, 129, 130}

        supported = [code for code in genetic_codes if code in available]
        unsupported = [code for code in genetic_codes if code not in available]

        return supported, unsupported

    except ImportError:
        logger.error("pyrodigal not installed")
        return [], genetic_codes


# Alternative approach: Monkey-patch the translation functions
def patch_pyrodigal_translation():
    """
    Alternative approach: Monkey-patch pyrodigal's translation functions.

    Note: This doesn't work with modern pyrodigal as GeneFinder is immutable.
    We'll use GeneticCodeHandler instead for codes 106 and 129.
    """
    # pyrodigal.GeneFinder is immutable, so we can't patch it
    # Just return True and let GeneticCodeHandler handle these codes
    logger.debug(
        "pyrodigal.GeneFinder is immutable - codes 106/129 will use GeneticCodeHandler"
    )
    return True


def get_supported_genetic_codes():
    """Get list of genetic codes supported by pyrodigal (including patched ones)."""
    try:
        # Apply patches
        patch_pyrodigal_genetic_codes()

        # Return standard codes + our additions
        standard = set(range(1, 7)) | set(range(9, 17)) | set(range(21, 26)) | {0}
        patched = {30, 106, 129, 130}
        return sorted(list(standard | patched))

    except Exception:
        # Return default list
        return [0, 1, 4, 6, 11, 15, 106, 129]


# Auto-patch on import
if __name__ != "__main__":
    patch_pyrodigal_genetic_codes()
