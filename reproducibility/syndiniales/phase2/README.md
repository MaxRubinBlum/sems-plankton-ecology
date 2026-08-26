# Syndiniales partner restructuring — Phase 2

Phase 2 tests whether Syndiniales–protist associations change across the **220-m boundary independently recovered in Phase 1**.

## Analysis logic

The workflow separates three questions:

1. **Network inference:** infer separate FlashWeave networks above and below 220 m.
2. **Network-selection robustness:** retest every candidate pair in both regimes so failure to be independently selected by a sparse network is not automatically interpreted as turnover.
3. **Formal regime interaction:** test partner abundance × depth regime while controlling continuous environmental and spatial covariates.

This distinction is essential because the two independently inferred sparse networks had almost no exact edge overlap, while direct cross-regime retesting recovered a substantial persistent association backbone.

## Input preparation

`prepare_phase2_inputs.py` starts from the unfiltered 18S ASV/taxonomy tables and integrated CTD metadata described in `../INPUT_MANIFEST.md`.

The archived preparation used 249 complete environmental cases. ASVs were aggregated to the most resolved informative label (Species → Genus → Family → Order → Class → Phylum), with broad ecological group included in the aggregation key. Parasite features were defined from Syndiniales/Dino-Group annotations and *Hematodinium*.

Initial filtering retained parasite features with prevalence ≥5% and total reads ≥500, and other eukaryotes with prevalence ≥3% and total reads ≥500. The table was then split at 220 m and filtered again within each regime. Regime-specific minimum total reads were `max(100, 2 × n_samples)`.

Archived inputs:

| Regime | Samples | Features | Parasite features |
|---|---:|---:|---:|
| ≤220 m | 143 | 572 | 58 |
| >220 m | 106 | 335 | 47 |

## FlashWeave

`run_phase2.jl` uses fixed parameters:

- FlashWeave-Sensitive
- `heterogeneous=false`
- `max_k=3`
- `alpha=0.01`
- `FDR=true`

Each depth regime is inferred independently with its own feature and metadata tables.

## Post-inference analysis

`analyze_phase2.py`:

- extracts parasite ↔ unicellular-eukaryote edges;
- normalizes broad partner-guild edge density by available partner features;
- retests all candidate pairs in both regimes using environmentally residualized abundances;
- distinguishes persistent support from sparse-network selection differences;
- fits a formal partner × regime interaction model for pairs testable on both sides of 220 m;
- applies Benjamini–Hochberg correction across all 207 interaction tests;
- tests broad partner-guild directional enrichment.

The formal model uses log1p abundance with HC3 robust standard errors and controls for continuous depth, CTD temperature, salinity, oxygen, fluorescence, season, and station. The partner × regime term is the inferential target.

## Archived result

The upper and deep FlashWeave networks contained 210 and 128 parasite–protist edges, respectively. Only one exact edge was independently selected in both, but cross-regime retesting recovered **95 associations with residual support in both regimes** (Jaccard = 0.332 across the retested union). Thus exact sparse-network overlap dramatically overstated turnover.

Among 207 pairs testable in both regimes, **89 showed a significant partner × regime interaction after FDR correction**: 53 stronger above and 36 stronger below 220 m. Broad guild-level directional enrichment was not significant after correction. The supported conclusion is therefore **lineage-specific restructuring of parasite–protist coupling across the boundary**, not wholesale replacement of one partner guild by another.

## Run order

```bash
python prepare_phase2_inputs.py \
  --otu /path/to/unfiltered_18S_otu.csv \
  --taxonomy /path/to/unfiltered_18S_taxonomy.csv \
  --metadata /path/to/18S_samples_CTD_chemistry.csv \
  --outdir work/phase2_input

julia --project=. run_phase2.jl work/phase2_input results/syndiniales_phase2/flashweave

python analyze_phase2.py \
  --input-dir work/phase2_input \
  --upper-edges results/syndiniales_phase2/flashweave/phase2_upper_le_220m.edgelist \
  --deep-edges results/syndiniales_phase2/flashweave/phase2_deep_gt_220m.edgelist \
  --outdir results/syndiniales_phase2

python make_summary_figure.py \
  --results results/syndiniales_phase2 \
  --out-prefix results/syndiniales_phase2/Supplementary_Figure_ParasitePartner_Phase2
```

## Interpretation guardrails

- FlashWeave edges are conditional associations, not demonstrated infection.
- Failure to select the same sparse edge twice is not sufficient evidence of biological turnover.
- Broad guild-level directional enrichment was non-significant after FDR; individual radiolarian or dinoflagellate examples must not be generalized into a community-wide host-guild shift.
- Candidate hosts should remain described as associated protists unless independently supported by microscopy, culture, FISH, single-cell data, or sequence-level host references.
