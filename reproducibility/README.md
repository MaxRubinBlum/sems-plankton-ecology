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
- Figure 5: cross-domain FlashWeave associations
- Supplement: chemistry-expanded hydrography/PCA and fine-resolution 16S seasonal succession
- Syndiniales supplement Phase 1: vertical parasite-community structure and objective boundary
- Syndiniales supplement Phase 2: depth-dependent restructuring of parasite–protist associations

## Syndiniales analyses

The complete parasite workflow is under `reproducibility/syndiniales/`.

- `INPUT_MANIFEST.md` records the exact unfiltered 18S and metadata inputs and SHA-256 checksums used for the archived analysis.
- `phase1/` contains the Syndiniales-only Bray–Curtis boundary scan, permutation/bootstrap specification, PERMDISP, clade-level depth analysis, and supplementary figure script.
- `phase2/` contains the unfiltered-18S input preparation, regime-specific FlashWeave inference, cross-regime retesting, formal partner × regime interaction model, guild enrichment test, and supplementary figure script.

The Phase-2 documentation explicitly distinguishes sparse-network selection differences from statistically supported ecological turnover. Broad partner-guild directional enrichment was tested and was not significant after FDR correction; individual associations must therefore not be generalized into a community-wide host-guild shift.

The input manifests distinguish raw/source inputs from intermediate products. Files not yet committed to `data/` are named explicitly so they can be deposited without changing the analytical specification.
