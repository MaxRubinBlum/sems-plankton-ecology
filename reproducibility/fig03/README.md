# Figure 3 — vertical taxonomic succession and trophic representation

## Purpose

Figure 3 explains the biological succession underlying the transition analysis in Figure 2. Panels a–b show the full 16S and 18S records, whereas panel c displays the five-component 18S trophic representation using the same within-five normalization used for the matched Figure 2 trophic scan.

## Taxonomic panels (a–b)

Feature counts are converted to sample-wise relative abundance. Features are aggregated primarily at Order; unresolved placeholder ranks fall back to Class and then Family. Display-label cleaning changes labels only, not aggregation or abundance.

Taxon selection is performed on the pooled marker-specific record before separating monitoring states. Median relative abundance is calculated for seven vertical habitat classes: Surface, Near-surface, DCM, Below DCM, 300–600 m, >600 m and Near-bottom. Taxa must reach at least 0.20% median relative abundance in one layer. Vertical differentiation is ranked as

`(maximum layer median - minimum layer median) * log(1 + 5 * maximum layer median)`.

The 14 highest-ranking groups are retained for each marker and ordered from shallow- to deep-associated using an abundance-weighted depth centroid with representative depths 1, 50, 115, 190, 375, 1000 and 1450 m.

Winter/mixed and summer/stratified layer medians are then calculated for this fixed pooled taxon set. Each taxon is standardized jointly across all 14 state × habitat values, so winter and summer are directly comparable within a taxon. The same standardized display scale is used for panels a and b. Bubble area indicates maximum median relative abundance across all states and layers.

## Trophic panel (c)

Five 18S-derived components are shown: bona fide phototrophy, constitutive mixotrophy, parasitism, radiolarian phagotrophy and diplonemid heterotrophy.

Before plotting, the five values are normalized among themselves within each sample to sum to 100%. This is the same selected-five normalization used before Hellinger transformation in the Figure 2 trophic-transition analysis. Figure 3 therefore decomposes the same trophic representation used to estimate the 110-m transition.

Continuous winter/mixed and summer/stratified depth profiles use a logarithmically spaced depth grid from 8 to 1700 m. At each focal depth the half-window is `max(35 m, 0.22 × focal depth)`. A summary is drawn only when at least six samples occur within the window. Curves show the local median and envelopes show the interquartile range.

Dotted horizontal references mark pooled matched-sample optima independently estimated in Figure 2:

- 110 m: trophic representation
- 220 m: shared matched 16S and 18S taxonomic composition

Monitoring-state-specific transition optima remain in Figure 2 and are not added to Figure 3 to avoid visual clutter.

## Reproducibility

`make_figure.py` is the authoritative implementation. It regenerates taxon selections, seasonal medians and z-scores, normalized trophic depth summaries, and the final horizontal Figure 3 from repository inputs.

The frozen numerical outputs in `results/fig03/` correspond to the current manuscript figure. Trait proportions are marker-gene-derived representation, not direct measurements of biomass, grazing, photosynthesis or carbon flux.
