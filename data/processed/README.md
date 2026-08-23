# Processed data manifest

Processed files are generated from the raw data under `data/raw/` and are not substitutes for the raw source tables.

Figure 4 requires:

- `data/processed/18S_samples_CTD_chemistry.csv` — integrated sample metadata containing at least `sample-id`, `depth`, `ds2`, `season`, `station`, `ctd_temperature`, `ctd_salinity`, `ctd_oxygen`, and `ctd_fluorescence`.

The validated Figure 4 analysis expects 308 samples shared among the 16S table, 18S table and this metadata table, with 248 complete cases for the five environmental variables used in the adjusted analysis.

The current working-session source for this processed table was `18S_samples_CTD_chemistry.csv` from the environmental-integration workflow. Preserve sample IDs exactly so they match the 16S and 18S abundance-table column names.
