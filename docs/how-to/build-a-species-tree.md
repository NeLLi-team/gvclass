# Build a species tree

The default run places each marker in its own tree and takes a majority vote across them. A species tree instead concatenates the core markers into one supermatrix and places the whole genome at once, which gives a stronger genome-level signal. This guide turns that on and reads the placement back out. For why the two methods differ, see [the species tree explanation](../explanation/species-tree.md).

## Build a per-query tree

Add `--species-tree` to a normal [classification run](classify-bins.md):

```bash
pixi run gvclass my_bins -o my_results -t 32 --species-tree
```

For each query, GVClass concatenates the classified domain's panel markers into a per-query supermatrix and builds that query's tree. For an NCLDV query, that panel is the eight GVOG8 markers. The placement is the nearest reference genome in the resulting tree.

The artifacts land in `my_results/species_tree/<query>/`:

| File | Contents |
| --- | --- |
| `<query>.treefile` | the per-query species tree, Newick |
| `<query>.partitions.txt` | the supermatrix partition map, one block per marker |
| `species_tree_taxonomy.tsv` | the placement result: nearest references and their taxonomy |

The run also fills four columns in `gvclass_summary.tsv`:

- `species_tree_nn_taxonomy`: taxonomy of the nearest reference in the tree.
- `species_tree_nn_genome`: that reference's genome id.
- `species_tree_nn_distance`: tree distance to it, six decimals.
- `species_tree_clade_id`: the clade it falls in.

!!! note

    Without `--species-tree`, all four `species_tree_*` columns read `nd`. The per-query placement owns these columns. See [Output](../reference/output.md) for the full 44-column schema.

## Also build one combined tree

To place every query in a single shared tree as well, add `--species-tree-combined`:

```bash
pixi run gvclass my_bins -o my_results -t 32 --species-tree-combined
```

`--species-tree-combined` implies `--species-tree`, so you get the per-query trees plus one combined tree over all queries, written to `species_tree/combined.*`. For a batch that spans more than one domain, the combined trees are split by panel under `species_tree/_combined/<panel>/`. The combined tree is an extra artifact for seeing how your queries relate to each other and to the references. It does not change the summary columns, which always come from the per-query placement.

## Choose a column trimmer

`--species-tree-trim` controls how supermatrix columns are pruned before tree building:

| Value | What it does |
| --- | --- |
| `witchi` | default, rigorous chi-squared column pruner |
| `pytrimal` | faster, lighter trimming |
| `none` | no column trimming |

```bash
pixi run gvclass my_bins -o my_results -t 32 --species-tree --species-tree-trim pytrimal
```

Keep `witchi` for a careful placement. Switch to `pytrimal` when you are trimming many large supermatrices and want the time back. Every flag is listed in the [CLI reference](../reference/cli.md).

## Cost and reproducibility

!!! warning

    Tree inference dominates wall-clock. That is why the species tree is opt-in, and why a plain classification run skips it. Budget extra time, more so with `--species-tree-combined`. See [Tune speed and accuracy](tune-speed-and-accuracy.md) for the levers.

!!! note

    Tree inference is multithreaded, so the topology, the selected reference set, and the six-decimal distances vary from run to run. The nearest-reference taxonomy, the placement itself, is stable across runs. For byte-identical trees, run with `--threads 1`.

## Resume and rebuild

On `--resume`, queries that already finished keep their existing `species_tree/<query>/` artifacts and are not re-placed. To rebuild a query's species tree, re-run it without `--resume`.
