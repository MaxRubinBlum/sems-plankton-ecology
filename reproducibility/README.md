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
- Supplement: chemistry-expanded hydrography/PCA and fine-resolution 16S seasonal succession

The input manifests below distinguish raw/source inputs from intermediate products. Files not yet committed to `data/` are named explicitly so they can be deposited without changing the analytical specification.