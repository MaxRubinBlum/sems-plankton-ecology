# Figure 4 — cross-domain community concordance

This directory reproduces the current Figure 4 analysis: Procrustes concordance between paired 16S and 18S rRNA gene community compositions across the full water column, after environmental/sampling adjustment, within vertical ecological regimes, and in the chemistry-complete sensitivity subset.

## Analysis

1. Intersect sample IDs present in the 16S table, 18S table and environmental metadata.
2. Convert each OTU table to sample-wise relative abundance.
3. Calculate Bray–Curtis dissimilarities separately for 16S and 18S.
4. Perform classical PCoA on each Bray–Curtis matrix and retain the first eight positive axes for the displayed analysis.
5. Compare the 16S and 18S configurations by symmetric Procrustes rotation and calculate the Procrustes correlation.
6. Assess significance by permutation (PROTEST-style permutation of sample correspondence).
7. For the adjusted analysis, regress the PCoA coordinates on log-depth, CTD temperature, salinity, dissolved oxygen, fluorescence, season and station, and perform Procrustes analysis on the residual configurations.
8. Repeat the comparison independently within the seven vertical ecological regimes. The adjusted within-regime model is fitted only where at least 20 complete environmental cases are available.
9. In the exact chemistry-complete paired subset (n = 67), compare the same hydrographic adjustment with an expanded adjustment that additionally includes nitrate + nitrite, phosphate, silicate and pH.

The near-bottom regime has 15 complete environmental cases and is therefore shown only as an unadjusted comparison.

## Main results

- Whole water column: n = 308, Procrustes r = 0.91, permutation P < 0.001.
- Adjusted complete-case dataset: n = 248, Procrustes r = 0.79, permutation P < 0.001.
- All adjusted within-regime comparisons shown in the figure have permutation P < 0.001.
- Chemistry-complete paired subset (n = 67): r = 0.835 after hydrographic/sampling adjustment and r = 0.833 after additionally including nitrate + nitrite, phosphate, silicate and pH; both permutation P = 0.0001.

The chemistry sensitivity uses the same 67 paired samples in both models, so the comparison isolates the effect of adding measured nutrient chemistry rather than changing sample composition. Measured nitrate + nitrite, phosphate, silicate and pH do not reduce the residual cross-domain concordance detectable in this subset.

A sensitivity analysis across PCoA dimensionality confirms that the result is not dependent on retaining eight axes. The adjusted correlation stabilizes at approximately r = 0.79 when 8–20 axes are retained.

## Figure conventions

Open circles represent 16S configurations and filled circles represent aligned 18S configurations. Lines connect the two representations of the same paired sample. To reduce visual clutter, the plotted figure displays a reproducible random subset of connecting lines, while all samples contribute to the Procrustes statistics. Panel a uses vertical-regime colors; the environment-adjusted panel uses neutral symbols. The final panel compares unadjusted and adjusted within-regime correlations.

## Interpretation

The analysis tests cross-domain community concordance, not direct interactions among individual taxa. Persistence of concordance after measured environmental adjustment, after adding nutrient chemistry in the chemistry-complete subset, and within vertical regimes indicates structured cross-domain covariance beyond the measured water-column gradient. These analyses do not establish causal or pairwise biological coupling, and unmeasured resource fields, particle dynamics, advection or spatial structure may still contribute to residual concordance.

## Reproduction

Run `python reproducibility/fig04/make_figure.py` from the repository root. The script reads the input paths documented in `inputs.tsv` and writes analysis summaries and figure files under `results/figure4_cross_domain_concordance/`, including `chemistry_sensitivity_summary.csv`.
