# Syndiniales partner restructuring — Phase 2

Phase 2 tests whether parasite–protist associations change across the **220-m boundary independently recovered in Phase 1**.

## Analysis logic

The workflow separates three questions:

1. **Network inference:** infer separate FlashWeave networks above and below 220 m.
2. **Network-selection robustness:** retest every candidate pair in both regimes so failure to be independently selected by a sparse network is not automatically interpreted as turnover.
3. **Formal regime interaction:** test partner abundance × depth regime while controlling continuous environmental and spatial covariates.

This distinction is essential because the two independently inferred sparse networks had almost no exact edge overlap, while direct cross-regime retesting recovered a substantial persistent association backbone.

## Input preparation

`prepare_phase2_inputs.py` starts from the unfiltered 18S 99%-OTU/taxonomy tables and integrated CTD metadata described in `../INPUT_MANIFEST.md`.

The archived preparation used 249 complete environmental cases. The 99%-OTUs were aggregated to the most resolved informative taxonomic label (Species → Genus → Family → Order → Class → Phylum), with broad ecological group included in the aggregation key. Parasite features were defined from the operational PR2 Syndiniales/Dino-Group annotations and *Hematodinium*.

Initial filtering retained parasite features with prevalence ≥5% and total reads ≥500, and other eukaryotes with prevalence ≥3% and total reads ≥500. The table was then split at 220 m and filtered again within each regime. Regime-specific minimum total reads were `max(100, 2 × n_samples)`.

| Regime | Samples | Features | Parasite features |
|---|---:|---:|---:|
| ≤220 m | 143 | 572 | 58 |
| >220 m | 106 | 335 | 47 |

## FlashWeave

`run_phase2.jl` uses fixed parameters: FlashWeave-Sensitive, `heterogeneous=false`, `max_k=3`, `alpha=0.01`, and `FDR=true`. Each depth regime is inferred independently with its own feature and metadata tables.

## Post-inference analysis

`analyze_phase2.py` extracts parasite ↔ unicellular-eukaryote edges, retests all candidate pairs in both regimes using environmentally residualized abundances, distinguishes persistent support from sparse-network selection differences, and fits a formal partner × regime interaction model for pairs testable on both sides of 220 m. The formal model uses log1p abundance with HC3 robust standard errors and controls for continuous depth, CTD temperature, salinity, oxygen, fluorescence, season, and station. The partner × regime term is the inferential target, with Benjamini–Hochberg correction across all 207 tests.

A non-metazoan CLR sensitivity analysis independently tests the regime-interaction result against compositional effects.

## Archived result

The upper and deep FlashWeave networks contained 210 and 128 parasite–protist edges, respectively. Only one exact edge was independently selected in both, but cross-regime retesting recovered **95 associations with residual support in both regimes** (Jaccard = 0.332 across the retested union). Thus exact sparse-network overlap dramatically overstated turnover.

Among 207 pairs testable in both regimes, **89 showed a significant partner × regime interaction after FDR correction**: 53 stronger at ≤220 m and 36 stronger at >220 m. In the non-metazoan CLR sensitivity analysis, 87/89 retained the same direction, 68/89 remained FDR-significant, and the interaction-effect changes were strongly concordant (Spearman rho = 0.889).

## Taxonomy audit and corrected guild test

An explicit audit against the original taxonomy table showed that the initial broad partner-guild heuristic misclassified or failed to classify many informative labels. In particular, several radiolarian families had been left under `Other protists`, while *Gymnoxanthella radiolariae* is a dinoflagellate despite its name. Broad guilds were therefore rebuilt from the original taxonomy ranks rather than partner-name substrings.

The correction changes the guild-level conclusion but **does not change any of the 207 formal pairwise interaction tests**. Corrected partner counts were: Radiolaria/Acantharea 67, Dinoflagellates 37, Other protists 35, MAST 24, Ciliates 12, with smaller numbers in the remaining guilds. Radiolaria/Acantharea were significantly enriched among associations that strengthened below 220 m (**odds ratio 3.80, Fisher P = 0.000651, BH-FDR q = 0.0169**). No other broad partner guild was enriched after multiple-testing correction.

A taxonomic-purity check confirmed that all corrected radiolarian labels used in this test had ≥0.8 class and phylum purity among their constituent 99%-OTUs (median class purity = 1.0).

The supported interpretation is therefore **lineage-specific restructuring of parasite–protist coupling across the boundary, including a statistically supported increase in radiolarian/Acantharean representation among deep-strengthened associations**, rather than wholesale replacement of the partner repertoire.

## Interpretation guardrails

- FlashWeave edges and abundance correlations are associations, not demonstrated infection.
- Failure to select the same sparse edge twice is not sufficient evidence of biological turnover.
- The radiolarian enrichment is a guild-level association result; it does not establish that every radiolarian partner is a host.
- Direct host evidence must be evaluated at the parasite-lineage level; host biology from MALV-II/*Amoebophrya* should not be transferred indiscriminately to MALV-I.
- Current phylogenomics separates MALV-I (Ichthyodinida) from MALV-II/IV (Syndiniales), although the operational PR2 labels used in this dataset retain Dino-Group I/II nomenclature.
