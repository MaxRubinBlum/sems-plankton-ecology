# Syndiniales vertical ecology — Phase 1

This directory reproduces the Syndiniales-only vertical analysis used to test whether parasite community composition has an objective depth transition.

## Inputs

The analysis uses the **unfiltered 18S table** (Metazoa retained), its PR2 taxonomy table, and the eukaryote sample metadata. These source files are not committed because the repository currently stores only the Metazoa-filtered 18S table. Exact source checksums used for the archived analysis are listed in `../INPUT_MANIFEST.md`.

Expected inputs are an unfiltered 18S ASV count table (`OTUID` rows, samples in columns), an unfiltered taxonomy table containing `Class`, `Order`, `Family`, `Genus`, and `Species`, and metadata containing `sample-id`, `select`, `station`, `season`, and `depth`. Only samples with `select == "a"` and non-missing depth are analyzed.

## Workflow

`run_phase1.py` performs the analysis:

1. Define Syndiniales from taxonomy (`Class` contains `Syndiniales`).
2. Calculate within-Syndiniales composition, ASV richness, Shannon diversity, and fraction of all 18S reads.
3. Calculate Bray–Curtis distances.
4. Scan candidate boundaries from 50–400 m at 10-m increments and estimate incremental partial R² for the upper/deep split after controlling station and season.
5. Test the maximum scan statistic using 199 depth permutations within station × season profiles.
6. Bootstrap station × season profiles 150 times to quantify breakpoint uncertainty.
7. Repeat the boundary scan independently by season while controlling station.
8. Produce Bray–Curtis PCoA and a profile-restricted permutation PERMDISP across the selected boundary.
9. Aggregate Syndiniales to the best-resolved clade label and calculate clade-specific depth associations.
10. `make_summary_figure.py` regenerates the supplementary figure from saved result tables.

Random seed: `20260826`.

## Archived result

The archived run contained 310 samples and 1,621 Syndiniales ASVs. The objective scan peaked at **220 m** (partial R² = 0.211697; max-statistic permutation P = 0.005). The profile bootstrap median was 220 m with a 95% interval of 220–240 m. Winter and summer independently peaked at 220 m. PERMDISP was significant (P = 0.006), so centroid and dispersion effects must be distinguished in interpretation.

## Run

```bash
python run_phase1.py \
  --otu /path/to/unfiltered_18S_otu.csv \
  --taxonomy /path/to/unfiltered_18S_taxonomy.csv \
  --metadata /path/to/metadata_euk.csv \
  --outdir results/syndiniales_phase1

python make_summary_figure.py \
  --results results/syndiniales_phase1 \
  --out-prefix results/syndiniales_phase1/Supplementary_Figure_Syndiniales_Phase1
```

## Interpretation guardrails

- The boundary is selected objectively from community dissimilarity, not chosen from individual taxa.
- The permutation P-value is based on the maximum statistic across the scan.
- The significant PERMDISP result is retained and reported.
- 18S relative abundance is not equivalent to parasite cell abundance because rRNA gene copy number varies among eukaryotes.
