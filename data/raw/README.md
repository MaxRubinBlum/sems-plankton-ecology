# Analysis-ready source data

This directory contains the cleaned source tables used by the repository workflows. Files are retained unchanged; filtering, integration and figure-specific transformations are written to `data/processed/`, `results/` or documented figure directories.

## Contents

```text
data/raw/
├── 16S/
│   ├── otu_table_prok_clean.csv.gz
│   ├── taxonomy_prok_clean.csv
│   ├── taxonomy_prok_GTDB226_crosswalk.csv
│   └── taxonomy_prok_GTDB226_rule_summary.csv
├── 18S/
│   ├── otu_table_euk_clean.csv.gz
│   └── taxonomy_euk_clean.csv
└── README.md
```

The original SILVA-based 16S taxonomy and GTDB release 226 crosswalk are retained separately and linked by feature identifier. The 16S and 18S features are 99% de-novo OTU clusters represented by centroid sequences; preserve their identifiers exactly.

Raw sequencing reads are not stored in this GitHub repository. Public sequence-archive accession details will be added to the main README when finalized.
