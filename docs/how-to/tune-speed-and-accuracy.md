# Tune speed and accuracy

GVClass ships with speed-first defaults: fast mode is on, VeryFastTree builds every tree, and sensitive marker search is enabled. Tuning happens through a few knobs: fast mode, tree method (including the species-tree search), thread layout, and marker sensitivity. Change one at a time, read the effect in the summary, then combine.

## Fast mode vs extended trees

Fast mode is the default. It skips the 576-marker order-level panel, which cuts runtime 2-3x at the cost of coarser order resolution.

```bash
# Default: fast mode on, order-level marker trees skipped
pixi run gvclass my_genomes -t 32
```

To build trees for every marker, pass `-e`/`--extended`. This turns fast mode off, sharpens the order assignment, and runs slower.

```bash
# Build all marker trees for a finer order call (slower)
pixi run gvclass my_genomes -t 32 --extended
```

`--mode-fast` forces fast mode back on, which only matters when a config file sets `mode_fast: false`.

!!! tip
    When `taxonomy_confidence` reports `reduced_fastmode`, that order was called without the order panel. Rerun the query with `--extended` if you need a firmer order. See the [output reference](../reference/output.md).

## Tree method

VeryFastTree builds all gene trees by default. For a more accurate placement, switch to IQ-TREE (model `Q.pfam+R10+F`), which is much slower.

```bash
# Default, fast
pixi run gvclass my_genomes -t 32 --tree-method veryfasttree

# IQ-TREE, more accurate, much slower
pixi run gvclass my_genomes -t 32 --tree-method iqtree
```

`fasttree` is an alias for `veryfasttree`.

## Species-tree support values

The species tree is opt-in (`--species-tree`) and built with VeryFastTree unless you also pass `--tree-method iqtree`. When IQ-TREE builds it, `--iqtree-mode` sets the search:

- `fast` (default): FastTree-like search, no bootstrap.
- `ufboot`: ultrafast bootstrap (`-B 1000 -bnni`), adds branch-support values, writes a `.contree` consensus, much slower.

```bash
# Final, well-supported species tree
pixi run gvclass my_genomes -t 32 --species-tree --tree-method iqtree --iqtree-mode ufboot
```

`--iqtree-mode` affects only the species tree. Per-marker gene trees always use fast search, since they are internal nearest-neighbor scaffolds where bootstrap adds cost and no signal. Full recipe: [build a species tree](build-a-species-tree.md).

## Threads and workers

`-t` is the total thread budget. `-j` is the number of queries run in parallel. `--threads-per-worker` sets the threads each worker gets. Set `-t` and one of the other two; GVClass derives the rest.

```bash
# 32 threads total, 4 queries at once, 8 threads each
pixi run gvclass my_genomes -t 32 -j 4
```

Use more workers when you have many small genomes; use more threads per worker for a few large ones. On a cluster, see [run on HPC](run-on-hpc.md).

## Sensitivity

Sensitive marker search is on by default. It uses `E=1e-5` (`domE=1e-5`) instead of the GA model cutoffs, recovering weaker marker hits at some risk of false positives. If a config turned it off (`sensitive_mode: false`), force it back on:

```bash
pixi run gvclass my_genomes -t 32 --sensitive
```

## Decision table

| Goal | Setting |
| --- | --- |
| Fastest triage | Defaults (fast mode + `--tree-method veryfasttree`) |
| Most accurate placement | `--extended --tree-method iqtree` |
| Well-supported final tree | `--species-tree --tree-method iqtree --iqtree-mode ufboot` |
| Maximise marker recovery | `--sensitive` (on by default, `E=1e-5`) |

For every flag and default, see the [CLI reference](../reference/cli.md). For why fast mode trades order resolution for speed, see [how GVClass works](../explanation/how-it-works.md).
