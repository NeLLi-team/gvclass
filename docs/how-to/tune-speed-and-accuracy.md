# Choose analysis settings

The defaults use fast mode, VeryFastTree, and sensitive HMM searches. Run the
commands below from the GVClass repository after database setup. They use the
bundled `example/` inputs; replace that path with your genome directory.

## Include order-level marker trees

Fast mode searches the order-level HMMs but skips their trees. To include those
trees, use `--extended`:

```bash
pixi run gvclass example -o example_extended -t 8 --extended
```

Inspect `taxonomy_majority`, `taxonomy_confidence`, and the per-rank vote columns
in `example_extended/gvclass_summary.tsv`. `reduced_fastmode` indicates that an
order assignment lacks the extended order-tree evidence.

## Use IQ-TREE

Add `--tree-method iqtree` to use IQ-TREE with the `Q.pfam+R10+F` model.
Per-marker trees use IQ-TREE's `--fast` search. Compare any change in taxonomy
with the per-marker votes in the summary.

## Species-tree support values

With `--tree-method iqtree`, select the species-tree search using:

| Option | Search |
| --- | --- |
| `--iqtree-mode fast` | Fast search without bootstrap replicates; the default. |
| `--iqtree-mode ufboot` | Ultrafast bootstrap with 1,000 replicates and `-bnni`. |

These settings affect species trees only. See [Build a species tree](build-a-species-tree.md)
for the run command and output paths.

## Threads and workers

`-t` sets the total thread budget. `-j` sets the number of queries processed at
once. For eight threads split across two workers, add `-t 8 -j 2` to a run; each
worker gets four threads. Alternatively, set `--threads-per-worker` and let
GVClass derive the worker count. Leave both options unset to use automatic
allocation.

Use only the CPUs available to the job. See [Run on HPC](run-on-hpc.md) for batch
submission.

## Marker sensitivity

Sensitive search is enabled by default and uses sequence and domain E-value
thresholds of `1e-5`. With `pipeline.sensitive_mode: false` in the configuration,
GVClass uses HMM gathering (GA) cutoffs where available. The `--sensitive` flag
overrides that setting for one run.

Use separate output directories when comparing settings. `--resume` skips
completed queries, so it does not recompute them with new options. See the
[CLI reference](../reference/cli.md) for all flags and the
[configuration reference](../reference/configuration.md) for persistent settings.
