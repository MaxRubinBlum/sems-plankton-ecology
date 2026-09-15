# Reproduction validation targets

These are the current frozen manuscript checkpoints. Regenerated analyses must be compared with them before publication assets are replaced.

## Figure 1

Hydrographic PCA of standardized CTD temperature, salinity, dissolved oxygen and fluorescence for 18S-associated samples:

- PC1: 57.1%
- PC2: 25.5%

The main figure also contains 16S and 18S Bray-Curtis PCoA panels, linear sampling depth, depth-layer colors and monitoring-state symbols.

## Figure 2

The authoritative analysis is `reproducibility/fig02/matched_transition_analysis.py`. It uses the 252 samples shared by the 16S, 18S and five-component trophic tables, restricted to <=650 m, representing 52 independent profiles.

Matched pooled optima:

- 18S trophic representation: 110 m; partial R2 = 0.475
- 16S composition: 220 m; partial R2 = 0.377
- 18S composition: 220 m; partial R2 = 0.219

Matched state-specific optima:

- winter/mixed: trophic 110 m; 16S 220 m; 18S 220 m
- summer/stratified: trophic 80 m; 16S 160 m; 18S 160 m

Profile-restricted, scan-wide maximum-statistic tests give P = 0.002 for all three layers.

The pooled joint profile bootstrap uses 500 replicates. Median optima are 110 m for trophic representation and 220 m for both taxonomic layers; 95% intervals are 70-120 m, 140-220 m and 140-220 m, respectively. Trophic reorganization is shallower than the 16S optimum in 100% of replicates and shallower than the 18S optimum in 99.8%; 16S and 18S select the same optimum in 80.4%.

Candidate thresholds span 10-420 m in 10-m increments. Displayed scans require at least 12 independent profiles spanning a threshold; support values of 10, 12 and 15 recover the same optima. Bootstrap scans retain the more conservative 15-profile requirement. State-specific bootstraps use 300 replicates.

The longer 2018-2026 16S record is contextual validation only and peaks at 170 m pooled, 220 m in winter/mixed conditions and 140 m in summer/stratified conditions.

## Figure 3

Taxon selection is pooled before division by monitoring state. Each marker retains 14 taxa after a 0.20% peak-median filter and vertical-differentiation ranking. State-specific row z-scores are calculated jointly across the 14 monitoring-state by depth-class values for each taxon. Five 18S trophic components are normalized among themselves within each sample.

## Figure 4

- Paired full dataset: n = 308; Procrustes r = 0.91; P < 0.001
- Hydrography-adjusted complete cases: n = 248; r = 0.79; P < 0.001
- Chemistry-complete sensitivity subset: n = 67
- Hydrographic/sampling adjustment: r = 0.835
- Adjustment additionally including nitrate + nitrite, phosphate, silicate and pH: r = 0.833
- Both chemistry-subset permutation tests: P = 0.0001

## Figure 5 and MALV analyses

- MALV community pooled transition: 220 m
- Shared cross-regime associations: 95; same sign: 93
- Significant partner by regime interactions: 89 of 207 testable pairs
- Radiolaria/Acantharea enrichment among associations strengthening below 220 m: odds ratio 3.80; Fisher P = 0.000651; BH-FDR q = 0.0169
- Environment-adjusted MALV versus non-parasitic unicellular-eukaryote Procrustes concordance: r approximately 0.766; P = 0.001

Associations represent conditional ecological covariance and candidate host-class relationships, not demonstrated infections.

## Policy

A workflow that fails a checkpoint must not silently overwrite frozen inputs or publication assets. First investigate preprocessing, transformations, nuisance-variable coding, profile definition, randomization restrictions and software versions; document any intentional methodological change.
