# Taxonomy and classification

`taxonomy_majority` is a consensus of the nearest-reference assignments from marker gene trees. Its resolution depends on marker support and the reference genomes available for the query's lineage.

## From marker trees to a lineage

Each marker contributes at most one taxonomic vote per rank, even when several copies occur in a query. If a protein matches overlapping marker models, its closest tree placement is counted once. At each rank, a marker votes for the taxon with the most placements among those under the parent lineage already selected. A marker whose two best taxa have the same number of placements casts no vote. It still counts in the total for that rank. The total at a rank is the number of markers with at least one placement at that rank under the selected parent.

GVClass then selects the taxon with the most marker votes and restricts the next rank to that parent lineage. If two taxa have the most votes, that rank and all lower ranks are left unassigned and the confidence includes `tie`, whether or not the tied taxa reach the support floor below. The result does not depend on the order of the markers or the files.

Assignments require two distinct supporting markers at most ranks. Order assignments require three in extended mode or two in fast mode. If a rank lacks support, that rank and all lower ranks are left unassigned. A rank is also left unassigned, with no flag, when no marker has a placement at that rank under the selected parent. This happens when the references of that lineage carry no label at the rank.

A winner is reported even when it holds half or fewer of the markers at its rank. In that case the confidence includes `plurality`. For example, three markers for one order and nine other markers that each name a different order give that order the call, flagged `plurality`.

A eukaryotic domain needs more than a plurality. Its markers must be at least 80% of the markers that vote for a cellular domain (bacterial, archaeal, or eukaryotic). Distant homologs and contaminated reference sequences give a prokaryotic genome a minority of eukaryotic votes, and a plurality of them can win. When eukaryotic markers win without reaching 80%, the domain and all lower ranks are left unassigned and the confidence is `mixed_cellular`. Bacterial, archaeal, and viral domains keep the plurality rule.

The detailed `domain` through `species` columns contain taxon counts from individual protein placements. These counts can differ from the distinct-marker votes used for `taxonomy_majority`. `avgdist` is the mean tree distance from query proteins to their nearest references.

## Confidence flags

| `taxonomy_confidence` | Meaning |
| --- | --- |
| `high` | No flag was set: there was no tie, enough markers supported every assigned rank, and more than half of the markers at every assigned rank agreed. |
| `no_support` | No usable taxonomic evidence was available. |
| `low_support` | Evidence was present, but too few distinct markers supported the leading taxon at a rank, so that rank and all lower ranks are unassigned. |
| `tie` | Two taxa had the most marker votes at a rank, so that rank and all lower ranks are unassigned. |
| `mixed_cellular` | Eukaryotic markers won the domain vote without reaching 80% of the cellular-domain markers, so the domain is unassigned. |
| `plurality` | A rank is assigned, but the winning taxon holds half or fewer of the markers at that rank. |
| `reduced_fastmode` | An order assignment relied on the relaxed fast-mode threshold of two markers. |

Multiple flags are comma-separated, in the order of the table: `no_support`, `low_support`, `tie`, `mixed_cellular`, `plurality`, `reduced_fastmode`. Fast mode alone does not cause `reduced_fastmode`: an order supported by three or more markers clears the extended-mode threshold too. `high` describes marker support, not a measured probability that the classification is correct.

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
