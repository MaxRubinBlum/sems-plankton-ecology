# Figure 2 — matched vertical transition analysis

Figure 2 compares the depths at which trophic representation, prokaryotic community composition and microbial-eukaryotic community composition reorganize in the southeastern Mediterranean water column.

## Current analytical design

The primary inference uses the exact sample universe shared by the 16S 99%-OTU table, the 18S 99%-OTU table and the five-component 18S trophic table. Samples deeper than 650 m are excluded from the transition scan. The matched dataset contains 252 samples from 52 independent sampling profiles.

For every candidate depth from 10 to 420 m in 10-m steps, the analysis calculates the incremental partial R² obtained by splitting samples above and below that depth after Hellinger transformation. Pooled scans account for station and monitoring state; winter/mixed and summer/stratified scans account for station. Profile identity is `cruise|station`; when cruise is missing, `year|month|station` is used so samples are not discarded or merged across years.

### Candidate-depth support

The displayed scans require at least **12 independent sampling profiles** spanning both sides of a candidate threshold. This value was chosen after explicit sensitivity analysis rather than to optimize a particular peak. Minimum-support values of 10, 12 and 15 profiles recover exactly the same transition optima for all three matched layers. Lowering the display criterion from 15 to 12 extends the seasonal scans far enough beyond the maxima to show the post-peak decline, especially for winter 18S, while retaining substantial independent-profile support.

The joint profile bootstrap retains the original more-conservative **15-profile** support rule. It resamples entire sampling profiles jointly across the three layers so pairwise differences between transition depths are estimated within the same bootstrap replicate.

### Long 16S record

A separate sensitivity scan uses the full 2018–2026 16S record. It applies the same transformations, threshold grid and nuisance structure but is **not** mixed into the matched three-layer bootstrap. The dashed 16S curve in panels a, c and d provides long-record context for the matched-period result.

### Scan-wide significance

Pooled significance is assessed using profile-restricted depth permutations and a maximum-statistic test across all admissible candidate depths. The same procedure is applied to 16S, 18S and trophic representation.

## Main results

Matched pooled optima are:

- 18S-derived trophic representation: **110 m**, partial R² = 0.475
- 16S composition: **220 m**, partial R² = 0.377
- 18S composition: **220 m**, partial R² = 0.219

Profile-restricted max-statistic tests give **P = 0.002** for all three layers.

The joint profile bootstrap gives median optima of 110 m for trophic representation and 220 m for both taxonomic layers. Trophic reorganization is shallower than the 16S optimum in 100% of bootstrap replicates and shallower than the 18S optimum in 99.8%; 16S and 18S select the same optimum in 80.4% of replicates.

State-specific matched optima are:

- winter/mixed: trophic 110 m; 16S 220 m; 18S 220 m
- summer/stratified: trophic 80 m; 16S 160 m; 18S 160 m

The longer 16S record independently peaks at 170 m pooled, 220 m in winter/mixed conditions and 140 m in summer/stratified conditions.

## DCM context

Transparent background fields in panels a, c and d show the independently derived distribution of deep-chlorophyll-maximum depths from full CTD fluorescence profiles. Up to two distinct fluorescence maxima are retained for casts with dual peaks. DCM information provides hydrographic context and does not enter transition estimation. The exact derived DCM table used to draw the figure is versioned at `results/fig02/dcm_peaks_from_ctd.csv`.

## Reproduction

Run from the repository root:

```bash
python reproducibility/fig02/matched_transition_analysis.py
python reproducibility/fig02/make_figure.py
```

The analysis script writes numerical outputs to `results/fig02/`. The plotting script reads those outputs plus the frozen CTD-derived DCM table and writes `Figure2_final_clean.svg`, `Figure2_final_clean.pdf` and `Figure2_final_clean.png`. All figure formats can be regenerated locally from the versioned analysis outputs and plotting script.

The older `recompute_boundaries.py` and `recompute_trait_boundary.py` document earlier unmatched analyses and should not be used as the authoritative implementation for the revised Figure 2.
