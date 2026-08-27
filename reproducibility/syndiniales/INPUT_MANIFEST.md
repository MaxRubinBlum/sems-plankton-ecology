# Input manifest for Syndiniales Phases 1–8

Exact source files used in the archived August 2026 analyses:

| Role | Source filename used locally | SHA-256 |
|---|---|---|
| Unfiltered 18S **99%-OTU** abundance table | `otu_table(20260826-074715).csv` | `d049bf7d59443269cbd064f57da0c8f71bbeb8dc864c00671a11841dc55db849` |
| Unfiltered 18S taxonomy for the 99%-OTU centroids | `taxonomy(20260826-074700).csv` | `a390dacdb0070989fd3bb3c2be3ca8371fceac0543505f87415b51f07995be3e` |
| Eukaryote metadata used in Phase 1 | `metadata_euk.csv` | `0b2a25569bb70975a0b319503effdf7ffd9ea18633bf2d4c930cccd31258bb40` |
| Integrated 18S CTD metadata used to prepare FlashWeave input | `18S_samples_CTD_chemistry.csv` | `6ca5abd68c0f56a9f4714af20b4eaa298685d70bc2be293fa85619ca7c3cd4b7` |
| Matching representative-sequence artifact used in Phase 6 | `rep-seqs-99.qza` | `751d987161559e47a7c39774567e4c966817741c331d7f58ec2803cb5c656670` |

## Feature provenance

QIIME 2 provenance in `rep-seqs-99.qza` shows that the original reads were denoised with `q2-dada2 denoise_paired` and hashed feature IDs, after which the features were clustered by `q2-vsearch cluster_features_de_novo` with `perc_identity = 0.99`.

Consequently, the current abundance table contains **99% de-novo OTU clusters represented by centroid sequences**. These features must not be described as unclustered ASVs. The sequence artifact contains 40,736 representative features and its feature-ID set exactly matches the 40,736 IDs in the abundance table.

The repository's existing `data/raw/18S/otu_table_euk_clean.csv.gz` and `taxonomy_euk_clean.csv` are **not interchangeable with these unfiltered files** for the complete parasite workflow because Metazoa were removed from the clean dataset. Phase 1 itself uses Syndiniales only for composition, but its reported fraction of total 18S reads also depends on the unfiltered denominator. Phase 2 was prepared from the unfiltered dataset before candidate partners were classified.
