# Figure 1 — vertical microbial organization

Figure 1 establishes the sampling structure and the broad vertical organization of prokaryotic and microbial-eukaryotic communities.

## Main panels

- Sampling design along the H transect, with linear depth, established vertical-habitat colors, and winter/summer symbols.
- Bray–Curtis PCoA of 16S community composition.
- Bray–Curtis PCoA of 18S community composition.
- 16S Shannon diversity versus depth.
- 18S Shannon diversity versus depth.

The CTD-only hydrographic PCA previously included in the main figure has moved to the supplementary hydrography workflow together with the chemistry-expanded ordination.

## Shannon diversity

Shannon diversity is calculated directly from the ASV abundance tables after converting each sample to relative abundance, using the natural-log definition H = -sum(p * ln p).

The main alpha-diversity curve combines winter and summer samples and uses robust LOWESS smoothing (`frac=0.28`, `it=2`). Thin season-specific curves are calculated independently for winter and summer (`frac=0.34`, `it=2`) and are displayed only as secondary context. The figure legend identifies the thick combined curve, thin solid winter curve, and thin dashed summer curve.

## Analytical requirements

- PCoA uses Bray–Curtis dissimilarity calculated from sample-wise relative abundances.
- PCoA variance is calculated from the positive eigenvalues, matching the validated Figure 1 checkpoint.
- Sampling depth uses a linear axis beginning at 0 m.
- Depth/habitat colors remain identical across sampling, PCoA, and Shannon panels.
- Winter and summer retain distinct point symbols.
- Shannon profiles use actual sample-level values; the smooth is descriptive and does not define ecological boundaries.
- No Figure 2 breakpoint or DCM-transition line is imposed on the Shannon panels, preserving their independence from the boundary analysis.

## Outputs to retain

- final Figure 1 SVG/PDF/600-dpi JPEG
- 16S and 18S PCoA scores and explained variance
- sample-level Shannon values
- LOWESS coordinates for combined, winter, and summer curves
- sampling-design table including station, depth, season, habitat class, and bottom depth

## Interpretation checkpoint

The alpha-diversity profiles provide an independent view of vertical organization. The 18S Shannon maximum occurs near the shallow eukaryotic turnover zone identified independently in Figure 2, whereas the 16S profile forms a broader high-diversity zone around the deeper prokaryotic transition. These patterns support, but do not define, the boundary analysis.
