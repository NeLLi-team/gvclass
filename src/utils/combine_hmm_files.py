#!/usr/bin/env python
"""
Combine HMM files using pyhmmer's native reading and writing capabilities.
"""
import pyhmmer.plan7
import pyhmmer.easel


def main():
    print("=" * 60)
    print("Combining HMM Files with pyhmmer")
    print("=" * 60)

    all_models = []

    # Read combined_old.hmm
    print("\n1. Reading combined_old.hmm...")
    with pyhmmer.plan7.HMMFile("resources/models/combined_old.hmm") as f:
        models_old = list(f)
    print(f"   ✅ Read {len(models_old)} models")
    all_models.extend(models_old)

    # Read OGv2_order.hmm
    print("\n2. Reading OGv2_order.hmm...")
    with pyhmmer.plan7.HMMFile("../models/OGv2_order.hmm") as f:
        models_ogv2 = list(f)
    print(f"   ✅ Read {len(models_ogv2)} models")
    all_models.extend(models_ogv2)

    print(f"\n3. Total models to write: {len(all_models)}")

    # Write combined file
    output_file = "resources/models/combined.hmm"
    print(f"\n4. Writing to {output_file}...")

    with open(output_file, "wb") as f:
        for model in all_models:
            model.write(f, binary=True)
    print(f"   ✅ Written {len(all_models)} models")

    # Validate
    print(f"\n5. Validating {output_file}...")
    with pyhmmer.plan7.HMMFile(output_file) as f:
        validated = list(f)
    print(f"   ✅ Successfully validated {len(validated)} models")

    if len(validated) == len(all_models):
        # Show first and last models
        print("\n   First 5 models:")
        for i in range(min(5, len(validated))):
            name = (
                validated[i].name.decode()
                if isinstance(validated[i].name, bytes)
                else validated[i].name
            )
            print(f"     {i+1}. {name}")

        print("\n   Last 5 models:")
        for i in range(max(0, len(validated) - 5), len(validated)):
            name = (
                validated[i].name.decode()
                if isinstance(validated[i].name, bytes)
                else validated[i].name
            )
            print(f"     {i+1}. {name}")

        # Quick test
        print("\n6. Running quick test...")

        with pyhmmer.plan7.HMMFile(output_file) as f:
            test_hmms = list(f)[:5]

        with pyhmmer.easel.SequenceFile("example/PkV-RF01.faa", digital=True) as f:
            test_seqs = list(f)[:5]

        hits = sum(
            1
            for top_hits in pyhmmer.hmmsearch(test_hmms, test_seqs, E=10.0)
            for hit in top_hits
        )
        print(f"   ✅ Test search with 5 models found {hits} hits")

        print("\n" + "=" * 60)
        print("✅ SUCCESS! Combined HMM file created")
        print(f"📁 File: {output_file}")
        print(f"📊 Total models: {len(validated)}")
        print("=" * 60)

        return True
    else:
        print(
            f"\n❌ Validation failed: expected {len(all_models)}, got {len(validated)}"
        )
        return False


if __name__ == "__main__":
    success = main()
    if not success:
        exit(1)
