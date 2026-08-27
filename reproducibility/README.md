# Figure reproducibility

This directory is the provenance layer for the manuscript figures. Each figure directory records the exact analytical inputs, derived tables, statistical decisions, and plotting workflow needed to recreate the publication figure from a clean checkout.

## Rules

1. Taxon/sample selection must be defined before plotting and must not be changed manually in vector graphics.
2. Seasonal comparisons use fixed pooled selections unless explicitly stated otherwise.
3. Any standardized heatmap must document the dimension over which scaling was calculated.
4. Derived values plotted in a figure should be exported as CSV/TSV whenever practical.
5. Figure scripts must use repository-relative paths and deterministic random seeds for resampling/permutation analyses.
6. Cosmetic post-processing must be avoided. If unavoidable, it must be documented and incorporated into the plotting script before manuscript release.
7. Manuscript legends and Methods should be traceable to the analytical decisions documented here.

## Current figure hierarchy

- Figure 1: hydrographic template
- Figure 2: objectively detected ecological boundary and seasonal displacement
- Figure 3: taxonomic succession and trophic reorganization, with winter/summer comparison
- Figure 4: cross-domain concordance
- Figure 5: cross-domain FlashWeave associations (candidate for replacement by the integrated Syndiniales ecological-reorganization figure)
- Supplement: chemistry-expanded hydrography/PCA and fine-resolution 16S seasonal succession
- Syndiniales supplement Phase 1: vertical parasite-community structure and objective boundary
- Syndiniales supplement Phase 2: depth-dependent restructuring of parasite–protist association strength
- Syndiniales Phase 3: partner-repertoire stability and strict turnover null test
- Syndiniales Phases 4–8: host-evidence grading, sequence-cluster validation, representative-sequence provenance, community concordance, and cross-ocean replication

## Syndiniales analyses

The complete parasite workflow is under `reproducibility/syndiniales/`.

- `INPUT_MANIFEST.md` records the exact unfiltered 18S and metadata inputs and SHA-256 checksums used for the archived analysis.
- `phase1/` contains the Syndiniales-only Bray–Curtis boundary scan, permutation/bootstrap specification, PERMDISP, clade-level depth analysis, and supplementary figure script.
- `phase2/` contains the unfiltered-18S input preparation, regime-specific FlashWeave inference, cross-regime retesting, formal partner × regime interaction model, guild enrichment test, compositional sensitivity analysis, and supplementary figure script.
- `phase3/` tests whether exact partner identities turn over more than expected after preserving parasite-specific candidate space, regime-specific availability, degree, and partner-guild composition.
- `phase4/` grades candidate-host evidence using quantitative recurrence plus explicitly coded literature precedent, while retaining an ecological-association category for unresolved links.
- `phase5/` tests whether selected taxon-aggregate associations persist at the constituent **99%-OTU** level.
- `phase6/` validates the matching representative-sequence artifact, extracts priority centroid sequences, and records the exact sequence provenance.
- `phase7/` tests concordance between Syndiniales and the non-Syndiniales unicellular-eukaryote community after environmental adjustment and verifies that the broader 18S boundary remains after Syndiniales are removed.
- `phase8/` compares pre-specified East Mediterranean Syndiniales depth trends against independently reported photic/aphotic clade patterns at BATS.

## Critical feature provenance

The uploaded `rep-seqs-99.qza` shows that the 18S workflow used `q2-dada2 denoise_paired` followed by `q2-vsearch cluster_features_de_novo` with `perc_identity = 0.99`. Thus the current OTU-table features are **99% de-novo clusters represented by centroid sequences**, not unclustered ASVs. Any Phase-5 or Phase-6 manuscript wording must use OTU/99%-OTU/centroid terminology rather than ASV terminology.

The artifact contains 40,736 representative features, exactly matching the current OTU-table feature IDs. All 78 Phase-6 priority sequence targets were recovered.

## Integrated interpretation

Phase 2 explicitly distinguishes sparse-network selection differences from statistically supported ecological change. Broad partner-guild directional enrichment was not significant after FDR correction; individual associations must not be generalized into a community-wide host-guild shift.

Phase 3 further shows that the apparent near-zero overlap of independently selected FlashWeave edges is not evidence for wholesale partner replacement. Cross-regime residual retesting recovered 95 shared associations, 93 of which retained the same sign. Under the strict parasite-specific candidate-universe null, observed exact partner overlap was indistinguishable from expectation. The supported interpretation is therefore **changes in association strength within a constrained partner repertoire**, rather than wholesale host/partner switching.

Phases 4–8 add independent support without upgrading co-occurrence to demonstrated infection: candidate associations are ranked using known host biology, selected aggregate signals recur at the 99%-OTU level, Syndiniales track the broader protist vertical architecture after environmental adjustment, and major photic/aphotic clade directions replicate independently reported BATS patterns.

The input manifests distinguish raw/source inputs from intermediate products. Files not yet committed to `data/` are named explicitly so they can be deposited without changing the analytical specification.
