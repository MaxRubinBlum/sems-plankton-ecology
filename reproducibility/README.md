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
- Figure 5 candidate: integrated MALV-parasite ecological-reorganization figure; the generic cross-domain FlashWeave figure can move to the Supplement if this architecture is retained
- Supplement: chemistry-expanded hydrography/PCA and fine-resolution 16S seasonal succession
- parasite supplement Phase 1: vertical parasite-community structure and objective boundary
- parasite supplement Phase 2: depth-dependent restructuring of parasite–protist association strength with taxonomy-corrected partner guilds
- Phase 3: partner-repertoire stability and strict turnover null test
- Phases 4–8: lineage-specific host-evidence grading, 99%-OTU validation, representative-sequence provenance/resolution audit, community concordance, and cross-ocean replication

## Parasite analyses

The complete workflow is under `reproducibility/syndiniales/`; the directory name is retained for continuity with the operational PR2 taxonomy. Current phylogenomics separates MALV-I (Ichthyodinida) from MALV-II/IV (Syndiniales), so manuscript text should use `MALV-I`, `MALV-II`, or `MALV parasites` when the analysis combines Group I and Group II.

- `INPUT_MANIFEST.md` records the exact unfiltered 18S and metadata inputs and SHA-256 checksums used for the archived analysis.
- `phase1/` contains the parasite-only Bray–Curtis boundary scan, permutation/bootstrap specification, PERMDISP, clade-level depth analysis, and supplementary figure script.
- `phase2/` contains unfiltered-18S input preparation, regime-specific FlashWeave inference, cross-regime retesting, formal partner × regime interaction model, compositional sensitivity analysis, taxonomy-corrected guild enrichment, and supplementary figure workflow.
- `phase3/` tests whether exact partner identities turn over more than expected after preserving parasite-specific candidate space, regime-specific availability, degree, and partner-guild composition.
- `phase4/` grades candidate-host evidence using quantitative recurrence plus **lineage-specific** independent host biology.
- `phase5/` tests whether selected taxon-aggregate associations persist at the constituent **99%-OTU** level.
- `phase6/` validates the representative-sequence artifact, extracts priority centroid sequences, and audits the taxonomic resolving power of the short marker.
- `phase7/` tests concordance between the parasite and non-parasitic unicellular-eukaryote communities after environmental adjustment and verifies that the broader 18S boundary remains after parasite features are removed.
- `phase8/` compares pre-specified East Mediterranean depth trends against independently reported photic/aphotic clade patterns at BATS.

## Critical feature provenance

The uploaded `rep-seqs-99.qza` shows that the 18S workflow used `q2-dada2 denoise_paired` followed by `q2-vsearch cluster_features_de_novo` with `perc_identity = 0.99`. Thus the current table features are **99% de-novo clusters represented by centroid sequences**, not unclustered ASVs. Any manuscript wording must use OTU/99%-OTU/centroid terminology rather than ASV terminology.

The artifact contains 40,736 representative features, exactly matching the current OTU-table feature IDs. All 78 Phase-6 priority sequence targets were recovered. The centroids are short (130–136 nt; median 132 nt), which is adequate for the operational community/clade analysis but insufficient for confident species- or strain-level host attribution.

## Integrated interpretation

**Phase 1:** the parasite community independently identifies a strong 220-m ecological transition. The same optimum is recovered in the full 310-sample analysis and the exact 249 complete cases used downstream. Six pre-specified MALV clades reproduce the independently reported BATS photic/aphotic depth directions.

**Phase 2:** association strength changes substantially across the transition. Among 207 testable parasite–protist pairs, 89 have FDR-significant partner × regime interactions (53 stronger at ≤220 m; 36 stronger at >220 m). A non-metazoan CLR sensitivity analysis retains the direction of 87/89 original significant effects and gives rho = 0.889 concordance in interaction-effect size.

A source-taxonomy audit corrected an important broad-guild classification error in the original exploratory analysis. With partner groups rebuilt from original taxonomy ranks, **Radiolaria/Acantharea are significantly enriched among associations that strengthen below 220 m** (odds ratio 3.80, Fisher P = 0.000651, BH-FDR q = 0.0169). No other guild survives multiple-testing correction. This replaces the earlier statement that no broad guild-level enrichment was detected.

**Phase 3:** altered association strength does not imply wholesale partner replacement. Cross-regime retesting recovers 95 shared associations, 93 with the same sign. Under the strict parasite-specific candidate-universe null, observed exact overlap is indistinguishable from expectation. The supported interpretation is therefore **reorganization of coupling strength within a comparatively constrained partner repertoire**.

**Phases 4–6:** host interpretation is lineage-specific. Host biology from MALV-II/*Amoebophrya* is not transferred indiscriminately to MALV-I. A particularly informative example is MALV-II clade 7: Kim et al. (2026) place *Amoebophrya sticholonchae*, directly observed infecting the radiolarian *Sticholonche zanclea*, in this clade. Our MALV-II-clade-7 association with the radiolarian *Acanthochiasma* strengthens below 220 m, providing a strong candidate-host-class hypothesis but not demonstrating infection of *Acanthochiasma*. Targeted deaggregation shows recurrence of several high-priority aggregate signals at the 99%-OTU level. The short centroid sequences prevent defensible species/strain-level phylogenetic host assignment.

**Phase 7:** parasite composition strongly tracks the broader non-parasitic unicellular-eukaryote community even after environmental, seasonal, and station adjustment (global adjusted Procrustes r ≈ 0.766, permutation P = 0.001), and the broader 18S boundary remains at 220 m after parasite features are removed. Parasites therefore track rather than generate the major vertical community transition.

The resulting narrative is: **shared vertical ecological architecture → parasite-community transition → depth-dependent restructuring of parasite–protist coupling → deep enrichment of radiolarian/Acantharean associations → largely persistent partner repertoire → a small set of lineage-informed candidate host relationships.** Statistical associations remain distinct from demonstrated infection throughout.
