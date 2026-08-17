# Reproduction validation targets

These values are frozen manuscript checkpoints. Recomputed analyses must be compared against them before replacing any publication asset.

## Figure 1

Hydrographic PCA, standardized CTD temperature + salinity + dissolved oxygen + fluorescence, 18S-associated samples:

- PC1: 57.1%
- PC2: 25.5%

Final Figure 1 additionally contains 16S and 18S Bray–Curtis PCoA panels, linear sampling depth, depth-layer colors and season symbols.

## Figure 2

Pooled strongest boundary:

- 18S trophic traits: 120 m
- 16S composition: 170 m
- 18S composition: 220 m

Season-specific strongest boundary:

- winter: traits 120 m; 16S 220 m; 18S 220 m
- summer: traits 120 m; 16S 140 m; 18S 160 m

Pooled max-statistic permutation P = 0.002 for all three data representations in the finalized analysis.

The extended pooled scan uses samples from 0–650 m and candidate thresholds from 10 to 420 m in 10-m increments. At 10 m the 16S analysis contains 134 shallow and 513 deep samples (647 total), providing a useful sample-filtering checksum.

Profile bootstrap: 250 replicates in the finalized output. Candidate thresholds require support from at least 15 independent cruise × station profiles spanning the threshold.

Frozen pooled bootstrap summaries:

- 16S median 170 m; 2.5–97.5% interval 140–170 m
- 18S median 220 m; 2.5–97.5% interval 140–220 m
- traits median 120 m; 2.5–97.5% interval 120–150 m

The independent deep-boundary control analysis contains 47 deep 18S samples from four stations, 17 near-bottom samples, partial R² = 0.163236 and station-stratified permutation P = 0.001 after accounting for station/bottom-depth context. This is a supporting robustness analysis, not the upper-boundary scan itself.

## Figure 3

Taxon selection is pooled before seasonal splitting. Main taxonomic panels retain 14 taxa per marker after a 0.20% peak-median filter and vertical-differentiation ranking. Winter and summer row z-scores are calculated jointly across the 14 season × depth values for each taxon. Five 18S trophic traits are summarized over continuous depth with local median and IQR windows.

## Policy

A script that does not reproduce these checkpoints must not silently overwrite the frozen figure inputs. Investigate preprocessing, transformation, nuisance-variable coding, profile definition, or randomization restrictions first, and document any intentional methodological change.