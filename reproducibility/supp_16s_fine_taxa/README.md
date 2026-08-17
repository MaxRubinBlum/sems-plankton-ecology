# Supplementary fine-resolution 16S seasonal succession

This supplement resolves Figure 3 prokaryotic succession below the main-text Order-level view.

## Analytical design

Feature counts are converted to within-sample relative abundance and aggregated using the terminal genus field after taxonomy validation. Generic unclassified labels, explicit higher-rank fallbacks, family names, and obvious group/clade placeholders are excluded from the strict genus-like set. Because reference taxonomies may nevertheless encode named environmental lineages in terminal fields, the manuscript should describe the displayed rows conservatively as **genus/lowest-resolved taxonomic groups**.

Selection is performed on the pooled dataset before seasonal splitting. Taxa must reach at least 0.10% median relative abundance in one vertical layer. The ranking statistic is `(max median - min median) * log(1 + 5 * max median)`. The top 30 taxa are retained and ordered by pooled abundance-weighted depth centroid.

Winter and summer medians are then calculated for this fixed taxon set. Each row is standardized jointly across all 14 season × depth values, not independently by season. Circles show maximum median relative abundance across seasons and depth layers.

## Required inputs

- cleaned 16S feature table
- cleaned 16S taxonomy table
- integrated 16S metadata with depth class and season

## Final asset

`Supplementary_16S_true_genus_winter_summer_30taxa` (SVG/PDF/600-dpi JPEG).

A taxonomy audit table should accompany the derived matrices so every displayed terminal label can be traced to the reference assignment.