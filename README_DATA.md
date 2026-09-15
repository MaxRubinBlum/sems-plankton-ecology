# SEMS plankton ecology data

This repository contains the cleaned abundance, taxonomy and environmental tables used by the documented analyses. These are analysis-ready source tables, not sequencing read archives.

## Amplicon tables

### 16S

- `data/raw/16S/otu_table_prok_clean.csv.gz`
- `data/raw/16S/taxonomy_prok_clean.csv`
- `data/raw/16S/taxonomy_prok_GTDB226_crosswalk.csv`
- `data/raw/16S/taxonomy_prok_GTDB226_rule_summary.csv`

### 18S

- `data/raw/18S/otu_table_euk_clean.csv.gz`
- `data/raw/18S/taxonomy_euk_clean.csv`

The abundance tables are gzip-compressed CSV files and can be read directly by pandas or R.

## Metadata

- `data/processed/metadata/metadata_16S_clean.csv`
- `data/processed/metadata/metadata_18S_clean.csv`
- `data/processed/metadata/16S_samples_CTD_chemistry.csv`
- `data/processed/metadata/18S_samples_CTD_chemistry.csv`

The integrated CTD and chemistry tables contain the environmental metadata used in the water-column analyses.

## Frozen Figure 4 inputs

- `data/processed/fig04/otu_table_prok_fig4_308.csv.gz`
- `data/processed/fig04/otu_table_euk_fig4_308.csv.gz`
- `data/processed/fig04/metadata_fig4_308.csv`

These preserve the exact 308-sample paired input state used for cross-domain concordance.

## Validation and provenance

- `reproducibility/VALIDATION.md` lists the current numerical checkpoints.
- Figure-specific `inputs.tsv` files identify analytical inputs.
- `reproducibility/syndiniales/INPUT_MANIFEST.md` records checksums and provenance for the MALV workflow.

Sample identifiers must remain unchanged because they link abundance tables to environmental metadata.
