# Explanation

These pages give the background and the reasoning behind GVClass: why the pipeline is shaped the way it is, how it turns marker placements into a taxonomic call, and what its quality numbers mean for giant viruses. They make for good reading away from the keyboard. None of it is needed to get a classification done, but it helps you read a result with the right amount of trust.

- [How it works](how-it-works.md): the pipeline end to end, from gene calling across nine genetic codes through HMM marker detection, per-marker trees, and the nearest-neighbor vote that produces a call.
- [Taxonomy](taxonomy.md): how the per-marker majority vote becomes a lineage, why genus and species are a nearest-reference label rather than an ICTV assignment, and how `taxonomy_confidence` is derived.
- [Quality metrics](quality-metrics.md): what estimated completeness and contamination mean for giant viruses, and why the many eukaryote-like genes acquired by horizontal transfer are not counted as contamination.
- [Species tree](species-tree.md): the opt-in concatenated-marker tree that adds a genome-level placement and fills the four `species_tree_*` columns.

!!! note

    For step-by-step recipes see the [how-to guides](../how-to/index.md), and for exact flags, columns, and panels see the [reference](../reference/index.md).
