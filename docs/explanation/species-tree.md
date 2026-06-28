# The species tree

Every GVClass run ends with a per-marker majority vote. For an NCLDV genome, each of the eight GVOG8 markers builds its own gene tree, the nearest reference leaf in that tree casts one vote for a lineage, and the votes are tallied rank by rank into `taxonomy_majority`. The vote is conservative and stable, because no single marker can carry a rank on its own and a contaminating or misplaced gene is outvoted. It also throws away information. Each gene tree sees only one marker's worth of aligned sites, and the votes are counted as if the markers were independent, so the joint signal across the panel never gets used.

The opt-in `--species-tree` route recovers that joint signal. For a genome GVClass has already routed to a registered viral domain (NCLDV, PPV, or MIRUS), it concatenates that domain's marker panel into one alignment with a separate partition per marker, then builds a single tree for the query against its reference neighbors. For NCLDV the panel is the eight [GVOG8 core markers](../reference/markers.md), and a genome needs at least three of them to enter the supermatrix; PPV and MIRUS use their own structural panels of capsid and packaging genes. The product is a supermatrix tree where the query sits among references on combined evidence rather than on eight separate scaffolds. GVClass reads the nearest reference leaf off that tree and fills four [summary columns](../reference/output.md): `species_tree_nn_taxonomy` (the lineage of that neighbor), `species_tree_nn_genome` (its genome id), `species_tree_nn_distance` (the tree distance to it, reported to six decimal places), and `species_tree_clade_id` (a clade identifier resolved from the tree). Without the flag, all four hold `nd`.

## Which references go in

The placement depends on which references enter the supermatrix. GVClass reuses the nearest neighbors it already computed during classification: from each marker's gene tree it keeps the top reference genomes, then pools them across the panel. Two knobs in `src/core/species_tree/config.py` set the breadth.

- `NEIGHBORS_PER_QUERY_TREE` (default 30): top references kept per marker tree for the per-query supermatrix, the larger set chosen for resolution.
- `NEIGHBORS_PER_COMBINED_TREE` (default 20): the lighter set used only when `--species-tree-combined` adds one tree across all queries at once (an extra artifact that does not change the summary columns).

Both values are read at run time, so editing one changes behavior with no other code change.

## Trimming the supermatrix

Before inference the concatenated columns are trimmed, and the trimmer is a deliberate trade of runtime against how aggressively non-phylogenetic signal is stripped. It is set by `--species-tree-trim`.

- `witchi` (default): a chi-squared column pruner that drops compositionally heterogeneous sites. The slowest and most aggressive option.
- `pytrimal`: a faster gap-and-similarity trimmer.
- `none`: keep every column.

For routine placement the default is the safe choice. On large batches where wall-clock matters, `pytrimal` trims the same alignment in a fraction of the time. See [Tune speed and accuracy](../how-to/tune-speed-and-accuracy.md) for the broader speed picture.

## Inference and the substitution model

The trimmed supermatrix goes to VeryFastTree by default, the same speed-first inference used for the per-marker gene trees. For a more accurate, fully model-based placement, pass `--tree-method iqtree`: GVClass then builds the species tree with IQ-TREE under `Q.pfam+R10+F`, the empirical Q.pfam amino-acid exchange matrix paired with a ten-category FreeRate model of among-site rate variation (`+R10`) and empirical amino-acid frequencies (`+F`). This model was not arbitrary. Across the GVMAGs V2 framework several substitution models were tested for Nucleocytoviricota species-tree inference, and `Q.pfam+R10+F` gave the most robust topologies across different datasets and marker-gene sets ([Vasquez et al. 2025](https://doi.org/10.1101/2025.09.26.678796)). The `--iqtree-mode` flag then chooses a quick `--fast` search or an ultrafast-bootstrap search; see [Tune speed and accuracy](../how-to/tune-speed-and-accuracy.md).

## What is reproducible, and what is not

The numbers around a placement are not deterministic, and it helps to be precise about which parts move. Tree inference runs multithreaded, and with more than one thread the topology, the exact set of reference genomes that end up in the supermatrix, and the six-decimal distances can all differ from one run to the next on identical input. The placement itself is stable to this. The nearest-reference taxonomy in `species_tree_nn_taxonomy` does not drift, even when the distance in `species_tree_nn_distance` changes in its trailing digits. When you need byte-identical output, for a regression test or a figure you have to reproduce exactly, run with `--threads 1`.

!!! note

    `--threads 1` is the only setting that guarantees identical species-tree output across runs. It is also the slowest, so reserve it for cases where exact reproducibility matters more than wall-clock.

## Why it is opt-in

The route is off by default for one reason: tree inference dominates the wall-clock. Building the per-group neighbor trees and then the supermatrix tree for each query is far slower than the HMM search and gene-tree votes that produce the standard call, so GVClass keeps it behind `--species-tree` and runs it only when you ask for a genome-level placement. A standard run, or a genome in no registered domain, never pays that cost and is unchanged.

Together the two products answer different questions. `taxonomy_majority` is the conservative consensus to trust for the domain-to-family backbone. `species_tree_nn_taxonomy` is a sharper, genome-level pointer at the single closest reference, useful when you want to pin a new GVMAG to the specific described genome it sits next to, finer than the family. Read them side by side. Agreement raises confidence, and a species-tree neighbor that refines the majority call to a finer clade is the common and informative case.

For the commands and the output layout, see [Build a species tree](../how-to/build-a-species-tree.md). For how the per-marker vote becomes `taxonomy_majority`, and why genus and species are a nearest-reference label rather than an ICTV assignment, see [Taxonomy](../explanation/taxonomy.md).
