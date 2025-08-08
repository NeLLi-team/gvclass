#!/usr/bin/env python3
"""
Test script to verify pyrodigal genetic code support.
Tests codes 1 (standard), 106, 129 (should work with fork), and 135 (should fail).
"""

import sys
import pyrodigal

# Test sequence - a small viral genome fragment with multiple ORFs
TEST_SEQUENCE = """
>test_sequence
ATGAAAATTGAAGAGATTTTAGATCATGAGAAGGAAATTGGCAAAACTGAAGGAGGAATTAGAGAGAGACTTACAAAAATTG
GCAGAAGCTATGAGAGAAATTTTTAAAGATCATGAGATGGAGCAGAACTAAAGATAGTAAGAGATGATATTGAATTAGATGAA
AAATTAAGACAGATTTTAGAAGAAGAGAAGGAAATGGCAAAAGAAGTAAAAGAAGTAAGAGATGAGATTGAGAGAGAAACTT
ACAAGAATTGGCAGAAGCTATGAGAGAAATTTTTAAAGATCATTAGATGGAGCAGAACTAAAGATAGTAAGAGATGATATT
GAATTAGATGAAAAATTAAGACAGATTTTAGAAGAAGAGAAGGAAATGGCAAAAGAAGTAAAAGAAGTAAGAGATGAGATT
GAGAGAGAAACTTACAAGAATTGGCAGAAGCTATGAGAGAAATTTTTAAAGATCATGAGATGGAGCAGAACTAAAGATAGT
AAGAGATGATATTGAATTAGATGAAAAATTAAGACAGATTTTAGAAGAAGAGAAGGAAATGGCAAAAGAAGTAAAAGAAGT
AAGAGGATGAGATTGAGAGAGAAACTTACAAGAATTGGCAGAAGCTTGATAGAAATTTTTAAAGATCATGAGATGGAGCAG
AACTAAAGATAGTAAGAGATGATATTGAATTAGATGAAAAATTAAGACAGATTTTAGAAGAAGAGAAGGAAATGGCAAAAG
AAGTAAAAGAAGTAAGAGATGAGATTGAGAGAGAAACTTACAAGAATTGGCAGAAGCTATGAGAGAAATTTTTAAAGATCAT
"""


def test_genetic_code(code: int) -> dict:
    """
    Test gene calling with a specific genetic code.

    Returns:
        dict with status, gene_count, and any error message
    """
    result = {
        "code": code,
        "status": "unknown",
        "gene_count": 0,
        "error": None,
        "genes": [],
    }

    try:
        # Clean sequence (remove header and newlines)
        seq_lines = TEST_SEQUENCE.strip().split("\n")
        sequence = "".join(seq_lines[1:])  # Skip the header line

        print(f"\n{'='*60}")
        print(f"Testing genetic code {code}")
        print(f"{'='*60}")

        # For testing, we'll use meta mode since it doesn't require training
        # and our test sequence is short
        print("  Using meta mode (no training required)")

        # Create gene finder in meta mode
        gene_finder = pyrodigal.GeneFinder(meta=True)

        # Try to use specific translation table in meta mode
        # First attempt: direct approach (works with newer versions)
        try:
            # This is how it would work if meta mode supported translation_table
            genes = gene_finder.find_genes(
                bytes(sequence, "utf-8"), translation_table=code
            )
            print("  ✓ Meta mode with translation_table parameter worked")
        except TypeError:
            # Fallback: Try training mode with the code
            print(
                "  Meta mode doesn't support translation_table, trying training mode..."
            )
            gene_finder = pyrodigal.GeneFinder(meta=False)

            # For testing, repeat the sequence to meet minimum length requirement
            extended_sequence = sequence * 30  # Make it >20000 chars
            print(f"  Extended sequence to {len(extended_sequence)} chars for training")

            # Train with the specific genetic code
            gene_finder.train(bytes(extended_sequence, "utf-8"), translation_table=code)

            # Find genes in original sequence
            genes = gene_finder.find_genes(bytes(sequence, "utf-8"))

        result["gene_count"] = len(genes)
        result["status"] = "success"

        print(f"✅ SUCCESS: Found {len(genes)} genes with code {code}")

        # Show first few genes (safely)
        try:
            for i, gene in enumerate(genes[:3]):
                print(
                    f"  Gene {i+1}: {gene.start}-{gene.end} ({'+' if gene.strand == 1 else '-'})"
                )
                result["genes"].append(
                    {
                        "start": gene.start,
                        "end": gene.end,
                        "strand": gene.strand,
                    }
                )
            if len(genes) > 3:
                print(f"  ... and {len(genes) - 3} more genes")
        except Exception:
            # Gene details might not be accessible, but that's ok
            pass

    except ValueError as e:
        result["status"] = "failed"
        result["error"] = str(e)
        print(f"❌ FAILED: {e}")

    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)
        print(f"⚠️  ERROR: {e}")

    return result


def main():
    """Run tests for different genetic codes."""

    print("Testing pyrodigal genetic code support")
    print(f"pyrodigal version: {pyrodigal.__version__}")
    print()

    # Test codes
    test_codes = [
        (1, "Standard genetic code - should work"),
        (106, "Code 106 (giant virus) - should work with fork"),
        (129, "Code 129 (giant virus) - should work with fork"),
        (135, "Code 135 - should fail (not available)"),
    ]

    results = []

    for code, description in test_codes:
        print(f"\n{description}")
        result = test_genetic_code(code)
        results.append(result)

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")

    success_count = sum(1 for r in results if r["status"] == "success")
    failed_count = sum(1 for r in results if r["status"] == "failed")

    print(f"\n✅ Successful: {success_count}")
    print(f"❌ Failed: {failed_count}")

    print("\nDetailed Results:")
    print(f"{'Code':<8} {'Status':<10} {'Genes':<10} {'Notes':<40}")
    print("-" * 68)

    for r in results:
        notes = ""
        if r["status"] == "success":
            notes = f"Found {r['gene_count']} genes"
        elif r["error"]:
            # Truncate error message if too long
            notes = r["error"][:37] + "..." if len(r["error"]) > 40 else r["error"]

        print(f"{r['code']:<8} {r['status']:<10} {r['gene_count']:<10} {notes:<40}")

    print("\n" + "=" * 60)

    # Verify expected behavior
    print("\nVERIFICATION:")

    checks = [
        (results[0]["status"] == "success", "✅ Code 1 (standard) works"),
        (
            results[1]["status"] == "success",
            (
                "✅ Code 106 works with fork"
                if results[1]["status"] == "success"
                else "❌ Code 106 FAILED - fork may not be installed correctly"
            ),
        ),
        (
            results[2]["status"] == "success",
            (
                "✅ Code 129 works with fork"
                if results[2]["status"] == "success"
                else "❌ Code 129 FAILED - fork may not be installed correctly"
            ),
        ),
        (
            results[3]["status"] == "failed",
            (
                "✅ Code 135 correctly fails (not available)"
                if results[3]["status"] == "failed"
                else "❌ Code 135 should have failed but didn't"
            ),
        ),
    ]

    for check, message in checks:
        print(message)

    # Overall test result
    all_correct = all(check for check, _ in checks)

    print("\n" + "=" * 60)
    if all_correct:
        print("🎉 ALL TESTS PASSED! The forked pyrodigal is working correctly.")
        print("   Codes 106 and 129 are supported.")
    else:
        print("⚠️  SOME TESTS FAILED!")
        if results[1]["status"] != "success" or results[2]["status"] != "success":
            print("   The forked pyrodigal may not be installed correctly.")
            print(
                "   Please ensure you're using: https://github.com/tomasbruna/pyrodigal.git"
            )

    return 0 if all_correct else 1


if __name__ == "__main__":
    sys.exit(main())
