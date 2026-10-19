# How the contig screen works

The contig screen answers one question about a giant-virus bin: which parts of it
are not from the giant virus. It combines two models. Neither replaces the
marker-based contamination estimate of GVClass; both are research-grade review aids.

## Cellular screen

The cellular screen scores overlapping windows of 10 kb at a 5 kb step. A window
carries three kinds of evidence:

- **Gene structure and composition:** gene density, gene length, strand switches,
  coding density and GC, codon and amino-acid usage.
- **Protein hits:** DIAMOND searches against a reference of bacterial, archaeal and
  eukaryotic proteomes, tailed phages, other viruses and the catalog of giant
  viruses. For each gene the model sees, per reference class, whether there is a hit,
  its identity and its bit score per amino acid, and which class holds the best hit.
  A gene without a hit is a zero token.
- **Protein embeddings:** the mean ESM-2 8M embedding of at most eight proteins per
  window (boosting model), and the embedding of every gene (set model).

Stage 1 averages the logits of a gradient-boosting model and a set model that pools
the gene tokens of a window. Stage 2 rescores each window against the other windows
of the same bin (distance to the bin centroid, neighbours in composition and
embedding space), so a window that differs from the bulk of its bin scores higher.
The thresholds come from the cross-validation of the training set: the stage 1
threshold passes at most 2% of giant-virus windows of 10 kb or more whose family was
held out; the stage 2 threshold flags at most 1% of clean synthetic bins.

The draw of the eight proteins is seeded by the sequence content of the bin, so the
file name and the order of the contigs do not change the draw.

## Giant-virus screen

The giant-virus screen looks for contigs of a second giant virus. It uses the copies
and the order of the catalog markers (OG and GVOGm), the identity of the copies,
the composition and gene structure of each contig against the rest of the bin, and
the best catalog hits of the genes and their agreement with the bin. A marker that
is single-copy in reference genomes but occurs on two contigs is evidence for a second
genome. Only contigs of 5 kb or more are scored.

## What the models were checked on

Numbers in this section are from the research record of the models (cross-validation
with the family of the test genome held out, and samples that no fit used). They
describe the 10% reference that ships with GVClass 2.1.0.

- **Unseen cellular genomes**, share of windows flagged at the stage 1 threshold, at
  3, 5, 10 and 20 kb: eukaryotes 0.96, 0.96, 0.98 and 0.99 (round 8 genomes); bacteria
  0.98, 0.99, 0.99 and 0.99; archaea 0.99, 0.99, 1.00 and 1.00.
- **Real bins**, contigs with a cellular origin verified by nucleotide alignment
  (2,981 contigs of 543 bins): stage 1 flags 214 of 232 contigs of 3 to 5 kb, 2,035
  of 2,130 of 5 to 10 kb and 517 of 533 of 10 to 20 kb; none of the 1,928 scored
  contigs with a verified viral origin is flagged. A stack fitted on the 17 GB full
  reference, with its own thresholds, flags 224, 2,101 and 532 of those contigs.
- **False positives on new viral genomes:** 1.2% of the windows of 3 kb of new
  giant-virus genomes and 0.1% at 5 kb. Herpesviruses and baculovirus-like viruses
  were not training negatives: 18.8% and 1.8% of their windows of 3 kb pass.

## Limits

- The labels of the real contigs are alignments to sequenced genomes, so the checks
  cover cellular contigs that have a relative in the references; 96% of them are
  eukaryotic.
- With the family of the virus held out, a lone eukaryotic window of 3 kb is found at
  0.75 when 1% of the viral windows pass (cross-validation, 10% reference). The checks
  on real bins and new genomes give higher rates for viruses of known families.
- The giant-virus screen cannot tell a second giant virus from a divergent region of
  the first one. It was checked on real bins with sequence evidence of origin, not on
  bins with a contamination rate known from read coverage.
- The 10% reference flags 3 to 4 points fewer verified cellular contigs of 3 to 20 kb
  than the full reference, and fewer cellular contigs in bins that are cellular
  throughout (stage 2: 1,543 against 2,069 of 2,981).
- The bundled completeness and contamination models of GVClass are unchanged.
