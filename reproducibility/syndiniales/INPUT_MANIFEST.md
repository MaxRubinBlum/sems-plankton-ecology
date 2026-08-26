# Input manifest for Syndiniales Phases 1–2

Exact source files used in the archived 26 August 2026 analysis:

| Role | Source filename used locally | SHA-256 |
|---|---|---|
| Unfiltered 18S ASV table | `otu_table(20260826-074715).csv` | `d049bf7d59443269cbd064f57da0c8f71bbeb8dc864c00671a11841dc55db849` |
| Unfiltered 18S taxonomy | `taxonomy(20260826-074700).csv` | `a390dacdb0070989fd3bb3c2be3ca8371fceac0543505f87415b51f07995be3e` |
| Eukaryote metadata used in Phase 1 | `metadata_euk.csv` | `0b2a25569bb70975a0b319503effdf7ffd9ea18633bf2d4c930cccd31258bb40` |
| Integrated 18S CTD metadata used to prepare FlashWeave input | `18S_samples_CTD_chemistry.csv` | `6ca5abd68c0f56a9f4714af20b4eaa298685d70bc2be293fa85619ca7c3cd4b7` |

The repository's existing `data/raw/18S/otu_table_euk_clean.csv.gz` and `taxonomy_euk_clean.csv` are **not interchangeable with these unfiltered files** for the complete parasite workflow because Metazoa were removed from the clean dataset. Phase 1 itself uses Syndiniales only for composition, but its reported fraction of total 18S reads also depends on the unfiltered denominator. Phase 2 was prepared from the unfiltered dataset before candidate partners were classified.
