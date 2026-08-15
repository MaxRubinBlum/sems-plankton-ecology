# Raw-data manifest

Raw inputs for this project are kept unchanged and should be added under the following paths:

```text
data/raw/
├── 16S/
│   ├── otu_table.csv
│   ├── taxonomy_SILVA.csv
│   └── taxonomy_GTDB226.tsv
├── 18S/
│   ├── otu_table.csv
│   └── taxonomy.csv
└── metadata/
    ├── metadata_16S.csv
    ├── metadata_18S.csv
    └── waterstations.xlsx
```

Current source filenames in the working analysis are:

- `otu_table_prok.csv`
- `taxonomy_prok.csv`
- `gtdb.tsv`
- `metadata_prok.csv`
- `otu_table_euk.csv`
- `taxonomy_euk.csv`
- `metadata_euk.csv`
- `waterstations.xlsx`

The two taxonomy systems for 16S must remain separate and be compared by ASV/feature ID. Processed filtering products belong in `data/processed/`, not here.
