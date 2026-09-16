# Species trees

GVClass can combine conserved marker alignments into a species tree containing query genomes and reference genomes. `--species-tree` builds a separate tree for each eligible query. `--species-tree-combined` also builds a shared tree of the eligible queries in each viral panel.

For commands and output paths, see [Build a species tree](../how-to/build-a-species-tree.md).

## Eligible queries and references

The domain in `taxonomy_majority` selects the marker panel. Queries and references must meet its minimum marker count.

| Panel | Marker groups | Minimum present |
| --- | --- | --- |
| NCLDV | 8 GVOG8 markers | 3 |
| PPV | 4 capsid and packaging groups | 2 |
| MIRUS | 6 groups covering MCP, ATPase, portal, and triplex proteins | 3 |

A query without a supported panel assignment or enough markers is excluded from species-tree inference. Mixed-panel input produces separate combined trees for NCLDV, PPV, and MIRUS.

References are selected from nearby genomes in each marker tree: up to 30 per marker for an individual query tree, and 20 per query and marker for a combined tree. The pooled references are then filtered by marker count.

Database bundles with putative endogenous viral element (pEVE) references can include `EUK-pEVE__...` genomes in these trees. They must meet the same marker threshold and retain their eukaryotic source labels. Ordinary `EUK__...` references are excluded.

## Alignment and inference

Each marker group is aligned and trimmed separately before concatenation. Missing markers are represented by gaps. The partition file records each marker's columns in the combined alignment.

`--species-tree-trim` selects the trimming method:

- `witchi` (default) prunes compositionally heterogeneous sites. If pruning fails, GVClass tries pytrimal, then the untrimmed alignment.
- `pytrimal` uses the `automated1` method, with the untrimmed alignment as fallback.
- `none` retains all aligned columns.

VeryFastTree infers the tree by default. `--tree-method iqtree` uses IQ-TREE with `Q.pfam+R10+F`. IQ-TREE species trees support either a fast search or 1,000 ultrafast bootstrap replicates; see [Tree settings](../how-to/tune-speed-and-accuracy.md).

## Interpreting the placement

`taxonomy_majority` summarises the separate marker-tree votes. The `species_tree_nn_*` columns give the nearest reference's taxonomy, genome identifier and tree distance in each query's concatenated tree. The combined tree has its own placement table and does not replace these values.

`species_tree_clade_id` identifies a clade containing multiple queries and no references. It is unset for per-query trees; the main summary reports `nd`.

Nearest-reference labels describe relationships to the sampled genomes. They do not establish an ICTV species assignment. Disagreement between the marker vote and the species tree warrants inspection of marker coverage, reference sampling, and individual gene trees.

See [Taxonomy and classification](taxonomy.md) for confidence flags and the [output reference](../reference/output.md) for column definitions.
