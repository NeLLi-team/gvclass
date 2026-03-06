# Benchmarking Runbook

This runbook records the local layout and commands used for reproducible benchmarking and release packaging.

## Current Layout

- Contamination benchmarking:
  - `benchmarking/contamination/contamination_benchmark`
  - `benchmarking/contamination/contamination_benchmark_faa`
- Completeness/reference workspace:
  - `benchmarking/completeness/refs-Feb-2026-fna`
  - `benchmarking/completeness/refs-Feb-2026-fna-pox10`
  - `benchmarking/completeness/refs-Feb-2026-fna-pox10-genus10`

## Contamination Model Export

Export the trained contamination bundle into runtime resources:

```bash
pixi run python -m src.bin.benchmark_contamination_methods_faa \
  --export-model resources/contamination_model.joblib
```

Expected outcome:

- `resources/contamination_model.joblib` exists
- bundle metadata reports `model_name = hist_gbm`

## Runtime Archive

Build a runtime archive containing only production assets:

```bash
tar -czf resources_v1_5_0.tar.gz resources
sha256sum resources_v1_5_0.tar.gz > resources_v1_5_0.tar.gz.sha256
```

## Private Benchmark Archive

Build a separate private benchmarking archive:

```bash
tar -czf benchmarking_v1_5_0.tar.gz benchmarking
sha256sum benchmarking_v1_5_0.tar.gz > benchmarking_v1_5_0.tar.gz.sha256
```

## Fresh-Start Local Validation

Validate the runtime archive without relying on the in-place `resources/` directory:

```bash
mkdir -p /tmp/gvclass_v1_5_0_test
tar -xzf resources_v1_5_0.tar.gz -C /tmp/gvclass_v1_5_0_test
pixi run gvclass example \
  -d /tmp/gvclass_v1_5_0_test/resources \
  -o /tmp/gvclass_v1_5_0_test/example_results \
  --threads 4 \
  --plain-output
```

## Private Transfer

Do not publish the benchmarking archive yet. If internal sharing is needed, copy it to private storage explicitly, for example:

```bash
scp benchmarking_v1_5_0.tar.gz user@private-host:/path/to/staging/
scp benchmarking_v1_5_0.tar.gz.sha256 user@private-host:/path/to/staging/
```

## Public Release Follow-Up

After the runtime archive is uploaded publicly:

1. Update the published download URL in:
   - `src/utils/database_manager.py`
   - `config/gvclass_config.yaml`
   - `config/gvclass_config_test.yaml`
2. Re-test a clean install using the published archive.
3. Push the repo release changes.
4. Publish the updated Apptainer image to Sylabs.
