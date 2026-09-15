# Processed data

Derived files in this directory are generated from the analysis-ready source tables under `data/raw/`.

## Metadata

Integrated metadata are stored under `data/processed/metadata/`. The Figure 4 workflow uses `data/processed/metadata/18S_samples_CTD_chemistry.csv`, containing sample identifiers, depth, vertical class, monitoring state, station and measured environmental variables.

## Frozen Figure 4 dataset

The exact paired Figure 4 inputs are stored under `data/processed/fig04/`:

- `metadata_fig4_308.csv`
- `otu_table_prok_fig4_308.csv.gz`
- `otu_table_euk_fig4_308.csv.gz`

The validated analysis contains 308 paired samples and 248 complete cases for the hydrographic adjustment. Preserve sample identifiers exactly so metadata and abundance-table columns remain aligned.
