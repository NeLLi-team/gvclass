# Build a species tree

Use `--species-tree-combined` to build a shared tree of query genomes and selected references. Run these commands from the GVClass repository after [installation and database setup](../tutorials/getting-started.md).

## 1. Prepare the input directory

Put one genome assembly in each `.fna` file, or one predicted proteome in each `.faa` file. To try the bundled nucleotide examples:

```bash
mkdir -p query_genomes
cp example/*.fna query_genomes/
ls query_genomes/
```

For your own data, replace the example files with your genome files. Do not combine separate genomes into one FASTA file.

## 2. Build the shared tree

```bash
pixi run gvclass query_genomes -o species_tree_results -t 8 --species-tree-combined
```

GVClass builds each query's tree, then a combined tree for each viral panel:

| Panel | Markers required |
| --- | --- |
| NCLDV | At least 3 of 8 |
| PPV | At least 2 of 4 |
| MIRUS | At least 3 of 6 |

Queries outside these panels or below the marker threshold are excluded.

## 3. Find the tree and placements

For a run with one eligible panel:

```bash
ls species_tree_results/species_tree/combined.*
head -n 6 species_tree_results/species_tree/species_tree_taxonomy.tsv
```

| File in `species_tree_results/species_tree/` | Contents |
| --- | --- |
| `combined.treefile` | Shared tree in Newick format |
| `combined.supermatrix.faa` | Concatenated protein alignment |
| `combined.partitions.txt` | Marker column ranges in the alignment |
| `species_tree_taxonomy.tsv` | Query placements from the shared tree |
| `<query>/<query>.treefile` | Tree for one query and its references |

For multiple panels, the combined files are under `species_tree/_combined/<panel>/`.

The `species_tree_*` columns in `gvclass_summary.tsv` contain per-query placements. Combined placements are in `species_tree/species_tree_taxonomy.tsv` for one panel, or `species_tree/_combined/<panel>/species_tree_taxonomy.tsv` for multiple panels.

## Build only per-query trees

```bash
pixi run gvclass query_genomes -o per_query_results -t 8 --species-tree
```

## Select the tree method or trimmer

The default uses VeryFastTree and `witchi` trimming.

| Option | Effect |
| --- | --- |
| `--tree-method iqtree` | Use IQ-TREE for tree inference |
| `--iqtree-mode ufboot` | With IQ-TREE, run 1,000 ultrafast bootstrap replicates for species trees |

See [tree settings](tune-speed-and-accuracy.md#species-tree-support-values) for output details.

To use `pytrimal` instead of `witchi`:

```bash
pixi run gvclass query_genomes -o pytrimal_results -t 8 \
  --species-tree-combined --species-tree-trim pytrimal
```

Use `--species-tree-trim none` to retain all alignment columns. Keep the alignment, tree, command, software version and database version with any reported analysis.

## Rebuild a combined tree

Run all queries into a new output directory, without `--resume`:

```bash
pixi run gvclass query_genomes -o rebuilt_species_tree_results -t 8 --species-tree-combined
```

Queries skipped by `--resume` are excluded from the combined tree for that run. If all queries are skipped, the combined tree is not rebuilt.

See [the species-tree method](../explanation/species-tree.md) for marker selection and reference selection.
