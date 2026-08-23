# SEMS plankton ecology data package

This package contains the complete cleaned data layer currently used by the analysis workflows.

## Full cleaned source data
- `data/raw/16S/otu_table_prok_clean.csv.gz`
- `data/raw/16S/taxonomy_prok_clean.csv`
- `data/raw/18S/otu_table_euk_clean.csv.gz`
- `data/raw/18S/taxonomy_euk_clean.csv`

The abundance tables are gzip-compressed CSV files. pandas and R read them directly.

## Metadata
- `data/processed/metadata/metadata_16S_clean.csv`
- `data/processed/metadata/metadata_18S_clean.csv`
- `data/processed/metadata/16S_samples_CTD_chemistry.csv`
- `data/processed/metadata/18S_samples_CTD_chemistry.csv`

The integrated CTD/chemistry tables are the environmental metadata used by the current water-column analyses.

## Frozen Figure 4 inputs
- `data/processed/fig04/otu_table_prok_fig4_308.csv.gz`
- `data/processed/fig04/otu_table_euk_fig4_308.csv.gz`
- `data/processed/fig04/metadata_fig4_308.csv`

These preserve the exact paired 308-sample input state used for the cross-domain concordance analysis.

`VALIDATION.json` records identifier/count checks made when this package was assembled.
`SHA256SUMS.csv` records file checksums.
