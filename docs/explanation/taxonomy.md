# Taxonomy and classification

`taxonomy_majority` is a consensus of the nearest-reference assignments from marker gene trees. Its resolution depends on marker support and the reference genomes available for the query's lineage.

## From marker trees to a lineage

Each marker contributes one taxonomic vote, even when several copies occur in a query. If a protein matches overlapping marker models, its closest tree placement is counted once. At each rank, GVClass selects the taxon supported by the most markers and restricts the next rank to that parent lineage.

Assignments require two distinct supporting markers at most ranks. Order assignments require three in extended mode or two in fast mode. If a rank lacks support, that rank and all lower ranks are left unassigned.

The detailed `domain` through `species` columns contain taxon counts from individual protein placements. These counts can differ from the distinct-marker votes used for `taxonomy_majority`. `avgdist` is the mean tree distance from query proteins to their nearest references.

## Confidence flags

| `taxonomy_confidence` | Meaning |
| --- | --- |
| `high` | No support warning was triggered. |
| `low_support` | Evidence was present, but too few distinct markers supported a rank. |
| `reduced_fastmode` | An order assignment relied on the relaxed fast-mode threshold of two markers. |
| `no_support` | No usable taxonomic evidence was available. |

Multiple flags are comma-separated. Fast mode alone does not cause `reduced_fastmode`: an order supported by three or more markers clears the extended-mode threshold too. `high` describes marker support, not a measured probability that the classification is correct.

## Reference labels

Genus and species fields can contain database identifiers such as `g1787` or `singleton`. They are reference labels rather than formal ICTV assignments. `s_singleton` is inherited from the reference taxonomy; it does not establish that the query is a new species.

The example PkV-RF01 lineage includes `o_Imitervirales;f_IM_19;g_g1787;s_singleton`. Its family assignment is IM_19, and the lower-rank identifiers describe its reference group. A short `avgdist` indicates a nearby reference in the marker trees; it is not a universal species boundary.

The [species-tree analysis](species-tree.md) adds a nearest-reference placement based on concatenated markers. It uses overlapping marker evidence and should be interpreted alongside `taxonomy_majority`.

## Viral panels and pEVE references

GVClass has species-tree panels for NCLDV (Nucleocytoviricota), PPV (Preplasmiviricota), and MIRUS (Mirusviricota). Classification also reports cellular and other reference lineages when supported; the species-tree panels do not define all possible classification results.

Some bundles include putative endogenous viral element (pEVE) references. Their identifiers start with `EUK-pEVE__`, and their taxonomy retains the eukaryotic source lineage with `-pEVE` labels. A source genome can supply both ordinary `EUK__` and pEVE proteins because this distinction is made per marker hit.

A pEVE match reports a relationship to one of these reference proteins. It does not by itself establish integration or assign a viral lineage. Eligible pEVE references can also appear in viral species trees while retaining their source taxonomy.

## Capsid typing

`capsid_group` reports MCP group counts, such as `Nucleocytoviricota:4,Gossevirus:1`. `plv` counts proteins matching the A32 ATPase marker `PLV_PC_054` whose tree placement supports PPV. A32 is shared with NCLDV, so an HMM match alone does not contribute to this count.

Use these counts with `vp_completeness`, `vp_mcp`, and the taxonomic assignment when examining virophages and Polinton-like viruses. See the [marker reference](../reference/markers.md), [output columns](../reference/output.md), and [genome-quality guide](../how-to/assess-genome-quality.md).
