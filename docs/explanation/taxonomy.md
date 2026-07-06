# Taxonomy and classification

Most giant viruses recovered from metagenomes have no named relative. Placing one in the tree of life means asking which known reference sits closest, marker by marker, and trusting that answer only as far as the references reach. GVClass treats classification as a phylogenetic placement problem and keeps the call deliberately conservative.

## From marker trees to a lineage

Classification rests on the per-marker trees the pipeline builds during a run (see [How it works](how-it-works.md)). Every conserved marker that a query carries gets its own single-gene tree, aligned against the top reference hits for that marker. From each tree GVClass reads the nearest reference neighbour of the query and records that neighbour's lineage. One marker, one vote.

The votes are tallied rank by rank. For each level from domain down to species, the lineage that the most markers point to becomes the call, and the per-taxon vote counts are kept beside it. The concatenated result lands in `taxonomy_majority`, while the individual rank columns (`species` through `domain`) carry the winning taxon with its count. `avgdist` records the mean tree distance from the query to those nearest references, a first read on how far the query sits from anything known.

For the bundled example genome PkV-RF01, the markers agree on `d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_IM_19;g_g1787;s_singleton`. The lineage resolves cleanly to family IM_19, then narrows to a single reference genus (g1787) and a species placeholder (s_singleton).

!!! note

    Genus and species identifiers like `g1787` and `s_singleton` are reference-derived labels from the GVClass database, not formal binomial names. `s_singleton` means the query placed on its own with no close species-level reference.

## Why confidence varies

`taxonomy_confidence` summarises how well the vote held up. A call is `high` only when every emitted rank cleared its distinct-marker threshold, meaning enough independent markers, not duplicate copies of one, agreed at that level. When a rank is emitted on too few distinct markers, the value drops to `low_support`. `reduced_fastmode` appears when fast mode (the default) skipped the order-level marker panel, so order resolution rests on the core set alone; build trees for all markers with `-e/--extended` to lift that flag. `no_support` marks a query whose markers could not agree on a lineage at all. PkV-RF01 reports `high`. Read this column together with the [quality metrics](quality-metrics.md), since duplicated markers inflate raw vote counts without adding independent evidence.

## Where the calls are trustworthy

The deep ranks are the dependable ones. Domain, phylum, class, order, and family draw on conserved markers whose reference trees are well sampled, so the majority vote at those levels is stable and reproducible. Genus and species behave differently. The pipeline reports the name of whatever genome sat closest in the marker trees, a nearest-reference label rather than an ICTV assignment. A genus or species call carries weight only when the query is close to a reference: a low `avgdist` and a real species match instead of an `s_singleton` placeholder. For a novel GVMAG, read the family and treat the genus and species as a pointer to the nearest known thing.

For genome-level confirmation at the tips, the opt-in [species tree](species-tree.md) concatenates a domain's core markers into one supermatrix and reports an independent nearest-reference placement in the four `species_tree_*` columns (see also the [how-to guide](../how-to/build-a-species-tree.md)).

## Registered domains

GVClass classifies into three registered viral domains, each with its own marker panel and reference set (full panel sizes are in the [markers reference](../reference/markers.md)):

- NCLDV (Nucleocytoviricota), the nucleocytoplasmic large DNA viruses, scored on the GVOG4 and GVOG8 core panels.
- PPV (Preplasmiviricota), the domain holding virophages and Polinton-like viruses, scored on the virophage core (MCP, Penton, ATPase, Protease) plus the A32/plv PPV-membership marker.
- MIRUS (Mirusviricota), scored on the four-marker Mirus core (MCP, ATPase, Portal, Triplex).

A genome that falls in none of these is left unclassified rather than forced into a domain.

## Putative EVE references

Some database bundles include eukaryotic marker proteins with viral signal in their nearest-neighbor evidence as putative endogenous viral element (pEVE) references. These references are not assigned to NCLDV, PPV, or MIRUS. They keep the eukaryotic source lineage in a separate namespace so the signal remains visible without pooling with ordinary eukaryotic reference labels.

Reference FASTA and label IDs use `EUK-pEVE__<source>` for pEVE proteins. Their taxonomy strings carry the namespace at every rank, for example `EUK-pEVE|Discosea-pEVE|Flabellinia-pEVE|EUK_unclassified-pEVE|Vannellidae-pEVE|Vannella-pEVE|Vannella sp.-pEVE`.

A source genome can contribute both ordinary `EUK__...` reference proteins and `EUK-pEVE__...` reference proteins. The distinction is made per protein marker hit, not per source genome.

When marker trees point to these references, output taxa carry `-pEVE` and the domain is `EUK-pEVE`. Read that as a eukaryotic source-lineage pEVE signal, not as a formal viral lineage assignment.

With `--species-tree`, pEVE references can also enter the NCLDV, PPV, and MIRUS species-tree candidate sets as auxiliary references. They are not routed as their own species-tree panel and ordinary `EUK__...` references stay excluded. A pEVE reference is kept only if it carries enough markers for that viral panel: three GVOG8 markers for NCLDV, two of four PPV groups, or three of six MIRUS groups. If it becomes the nearest species-tree reference, the `species_tree_nn_taxonomy` field reports its `EUK-pEVE` lineage rather than converting it into an NCLDV, PPV, or MIRUS assignment.

## Capsid typing

Two columns describe the major capsid protein (MCP) signal. `capsid_group` is a `label:count` tally across the MCP panels, for example `Nucleocytoviricota:4,Gossevirus:1`, spanning the Nucleocytoviricota and Mirusviricota phyla and the Bellas & Sommaruga capsid groups.

Within PPV, `plv` tells the two members apart. It counts A32 ATPase proteins (the `PLV_PC_054` marker) that place with PPV references in the marker tree. A32 is shared with NCLDV, so only copies that branch inside PPV are counted, which keeps `plv` at 0 for an ordinary NCLDV genome. Among PPV members, virophages and Polinton-like viruses both carry the virophage MCP markers (`vp_completeness`, `vp_mcp`), and a higher `plv` count leans the call toward a Polinton-like virus. The count is graded, not a yes/no flag.

## Reading the result conservatively

Together these choices keep the output honest about its own reach. GVClass reports the deepest lineage the markers support, records how much support that was, and labels the tips as nearest references rather than formal species. Read the family with confidence, read the genus and species as the closest known neighbour, and check `taxonomy_confidence` and `avgdist` before trusting either. Every column named here is defined in the [output reference](../reference/output.md).
