# Taxonomy and classification

`taxonomy_majority` is a consensus of the nearest-reference assignments from marker gene trees. Its resolution depends on marker support and the reference genomes available for the query's lineage.

## From marker trees to a lineage

Each marker contributes at most one taxonomic vote per rank, even when several copies occur in a query. If a protein matches overlapping marker models, its closest tree placement is counted once. The closest placement can be a tie: several references at exactly the same tree distance. The protein then votes with the lineage those references share, and it casts no vote at the ranks below it. A protein tied between two genomes that share a family but not a genus votes for the family only. At each rank, a marker votes for the taxon with the most placements among those under the parent lineage already selected. A marker whose two best taxa have the same number of placements casts no vote. It still counts in the total for that rank. The total at a rank is the number of markers with at least one placement at that rank under the selected parent.

GVClass then selects the taxon with the most marker votes and restricts the next rank to that parent lineage. If two or more taxa share the most votes, the confidence includes `tie:<rank>`. At the domain rank, a tie between exactly one viral domain and one or more cellular domains selects the viral domain, subject to the usual support floor. Lower ranks are then judged as usual within that viral domain. All other ties leave that rank and all lower ranks unassigned, whether or not the tied taxa reach the support floor below. This includes ties between viral domains, ties between cellular domains only, and every tie below domain. The result does not depend on the order of the markers or the files.

Assignments require two distinct supporting markers at most ranks. Order assignments require three in extended mode or two in fast mode. If a rank lacks support, that rank and all lower ranks are left unassigned. A rank is also left unassigned, with no flag, when no marker has a placement at that rank under the selected parent. This happens when the references of that lineage carry no label at the rank.

A winner is reported even when it holds half or fewer of the markers at its rank. In that case the confidence includes `plurality:<rank>`. For example, three markers for one order and nine other markers that each name a different order give that order the call, flagged `plurality:order`.

A eukaryotic domain needs more than a plurality. Its markers must be at least 80% of the markers that vote for a cellular domain (bacterial, archaeal, or eukaryotic). Distant homologs and contaminated reference sequences give a prokaryotic genome a minority of eukaryotic votes, and a plurality of them can win. When eukaryotic markers win without reaching 80%, the domain and all lower ranks are left unassigned and the confidence is `mixed_cellular`. Bacterial, archaeal, and viral domains keep the plurality rule. A cellular-only domain tie stops with `tie:domain` before this share check. A viral winner of a domain tie is exempt from the share check and does not set `mixed_cellular`.

The detailed `domain` through `species` columns contain taxon counts from individual protein placements. A protein with tied nearest references adds counts for the ranks those references share. These counts can differ from the distinct-marker votes used for `taxonomy_majority`. `avgdist` is the mean tree distance from query proteins to their nearest references.

## Confidence flags

| `taxonomy_confidence` | Meaning |
| --- | --- |
| `high` | No flag was set: there was no tie, enough markers supported every assigned rank, and more than half of the markers at every assigned rank agreed. |
| `no_support` | No usable taxonomic evidence was available. |
| `low_support` | Evidence was present, but too few distinct markers supported the leading taxon at a rank, so that rank and all lower ranks are unassigned. |
| `tie:<rank>` | Two or more taxa shared the most marker votes. That rank and all lower ranks are unassigned, except that exactly one viral domain tied with cellular domains wins the domain vote if it meets the support floor. The flag remains, and lower ranks are judged within that domain. |
| `mixed_cellular` | Eukaryotic markers won the domain vote without reaching 80% of the cellular-domain markers, so the domain is unassigned. |
| `plurality:<rank>` | A rank is assigned, but the winning taxon holds half or fewer of the markers at that rank. `<rank>` is the shallowest such rank. |
| `reduced_fastmode` | An order assignment relied on the relaxed fast-mode threshold of two markers. |

`<rank>` is `domain`, `phylum`, `class`, `order`, `family`, `genus`, or `species`. When `tie` or `plurality` applies at several ranks, each flag names the shallowest one. Thus `tie:domain` can also cover an unresolved tie below the selected viral domain, and `plurality:phylum` can hide a weaker call at the genus. Multiple flags are comma-separated, in the order of the table: `no_support`, `low_support`, `tie:<rank>`, `mixed_cellular`, `plurality:<rank>`, `reduced_fastmode`. Fast mode alone does not cause `reduced_fastmode`: an order supported by three or more markers clears the extended-mode threshold too. `high` describes marker support, not a measured probability that the classification is correct.

## Reference labels

Genus and species fields can contain database identifiers such as `g1787` or `singleton`. They are reference labels rather than formal ICTV assignments. `s_singleton` is inherited from the reference taxonomy; it does not establish that the query is a new species.

The example PkV-RF01 lineage includes `o_Imitervirales;f_IM_19;g_NCLDV_unclassified;s_PkV-RF01`. Its family assignment is IM_19, and the lower-rank identifiers describe its reference group. The database holds this genome under two identifiers with different genus and species labels (`g1787` and `singleton`, against `NCLDV_unclassified` and `PkV-RF01`). A marker tree that places the query at the same distance to both copies casts no vote at the genus and species ranks, so those calls come from the trees that hold only one copy. A short `avgdist` indicates a nearby reference in the marker trees; it is not a universal species boundary.

The [species-tree analysis](species-tree.md) adds a nearest-reference placement based on concatenated markers. It uses overlapping marker evidence and should be interpreted alongside `taxonomy_majority`.

## Viral panels and pEVE references

GVClass has species-tree panels for NCLDV (Nucleocytoviricota), PPV (Preplasmiviricota), and MIRUS (Mirusviricota). Classification also reports cellular and other reference lineages when supported; the species-tree panels do not define all possible classification results.

Some bundles include putative endogenous viral element (pEVE) references. Their identifiers start with `EUK-pEVE__`, and their taxonomy retains the eukaryotic source lineage with `-pEVE` labels. A source genome can supply both ordinary `EUK__` and pEVE proteins because this distinction is made per marker hit.

A pEVE match reports a relationship to one of these reference proteins. It does not by itself establish integration or assign a viral lineage. Eligible pEVE references can also appear in viral species trees while retaining their source taxonomy.

## Capsid typing

`capsid_group` reports MCP group counts, such as `Nucleocytoviricota:4,Gossevirus:1`. `plv` counts proteins matching the A32 ATPase marker `PLV_PC_054` whose tree placement supports PPV. A32 is shared with NCLDV, so an HMM match alone does not contribute to this count.

Use these counts with `vp_completeness`, `vp_mcp`, and the taxonomic assignment when examining virophages and Polinton-like viruses. See the [marker reference](../reference/markers.md), [output columns](../reference/output.md), and [genome-quality guide](../how-to/assess-genome-quality.md).
