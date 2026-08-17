# Figure 3 — seasonal vertical taxonomic succession and trophic organization

## Purpose

Figure 3 links the objectively detected boundary shifts in Figure 2 to the organisms and trophic strategies contributing to vertical ecological reorganization.

## Inputs

See `inputs.tsv`. Core inputs are cleaned 16S and 18S feature tables and taxonomies, integrated sample metadata containing depth class and season, and the sample-level 18S trophic-trait table.

## Taxonomic panels

### Relative abundance

For each marker dataset, feature counts are converted to sample-wise percentages by dividing each feature count by the total retained counts in that sample and multiplying by 100.

### Taxonomic aggregation

Features are aggregated primarily at Order. If Order is absent, Class and then Family are used as fallback labels. This preserves interpretable terminal groups without discarding unresolved but abundant lineages.

### Fixed taxon selection

Taxon selection is performed on the pooled dataset, before splitting samples by season. For each taxon, median relative abundance is calculated for each of seven vertical habitat classes:

1. Surface
2. Near-surface
3. DCM
4. Below DCM
5. 300–600 m
6. >600 m
7. Near-bottom

Taxa must attain at least 0.20% median relative abundance in one layer. Vertical differentiation is ranked as:

`(maximum layer median - minimum layer median) * log(1 + 5 * maximum layer median)`

The 14 highest-ranking taxa are retained for each marker dataset. This selection is independent of season.

### Taxon ordering

Rows are ordered from shallow- to deep-associated taxa using an abundance-weighted depth centroid. Representative layer depths used for ordering are 1, 50, 115, 190, 375, 1000 and 1450 m.

### Winter/summer comparison

After the fixed pooled selection and ordering, median relative abundance is recalculated separately for winter and summer within each vertical layer. The same taxa and row order are therefore displayed in both seasons.

### Heatmap scaling

For each taxon, winter and summer layer medians are concatenated and standardized jointly. Thus each row is a z-score across 14 values (7 layers × 2 seasons), allowing direct seasonal comparison. Winter and summer are NOT independently standardized.

### Abundance circles

The circle adjacent to each taxon represents its maximum median relative abundance across all depth layers and both seasons. Circle area uses square-root scaling for readability. A shared size key converts circle sizes back to percent relative abundance.

## Trophic-trait panels

Five 18S-derived trophic traits are displayed: bona fide phototrophy, constitutive mixotrophy, parasitism, radiolarian phagotrophy and diplonemid heterotrophy.

Sample-level trait values are plotted against continuous sampling depth separately for winter and summer. Local depth summaries use a moving depth window centered on a logarithmically spaced depth grid from 8 to 1700 m. The half-window is `max(35 m, 0.22 × focal depth)`. A summary is drawn only where at least six observations occur within the window. The line is the local median and the envelope is the interquartile range (25th–75th percentile). Faint points show individual observations. Winter is solid and summer dashed.

Trait panels use independent horizontal scales because the traits differ strongly in absolute prevalence. Depth is plotted on a logarithmic axis.

Horizontal dotted references at approximately 120, 170 and 220 m correspond to ecological transition depths independently identified by Figure 2 analyses; they are not fitted from the Figure 3 trait data.

## Output

Final manuscript asset:

- `Figure3_vertical_succession_winter_summer_split_v2.svg`
- corresponding vector PDF
- 600-dpi JPEG for manuscript drafting

The v2 layout includes the abundance-size key and spacing correction for the longest 18S taxonomic label.

## Reproducibility status

The complete analytical specification and inputs are now versioned. `make_figure.py` reproduces the analytical panels and should be treated as the authoritative implementation; any later aesthetic change should be made there rather than manually in the SVG.