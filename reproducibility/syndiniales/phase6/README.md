# Representative-sequence validation — Phase 6

Phase 6 validates the matching 18S representative-sequence artifact and evaluates how far the short marker can support candidate-host interpretation.

## Representative-sequence input

Uploaded artifact: `rep-seqs-99.qza`

SHA-256: `751d987161559e47a7c39774567e4c966817741c331d7f58ec2803cb5c656670`

QIIME 2 type: `FeatureData[Sequence]`

The artifact contains 40,736 representative sequences. Its feature-ID set exactly matches the 40,736 features in the current 18S OTU table.

## Provenance

QIIME 2 provenance shows two successive steps:

1. `q2-dada2 denoise_paired`, with hashed feature IDs.
2. `q2-vsearch cluster_features_de_novo`, with `perc_identity = 0.99`.

The current abundance table therefore consists of **99% de-novo OTU clusters represented by centroid sequences**, not unclustered ASVs.

## Target extraction

All 78 requested high-priority OTUIDs were recovered; none were missing and all 78 centroid sequences were unique. Sequence lengths were 130–136 nt (median 132 nt).

`extract_repseq_targets.py` validates the artifact against the OTU table, checks provenance, and extracts a taxonomy-linked target FASTA and metadata table.

## Sequence-resolution audit

Because these centroids are only ~132 nt, Phase 6 does not treat them as sufficient for strain- or species-level phylogenetic placement. A local nearest-sequence audit against all 40,736 centroids showed:

- median nearest-neighbour normalized sequence similarity = 0.985;
- 85.9% of priority centroids had a top neighbour in the same annotated class;
- 79.5% had a top neighbour in the same annotated order.

Several parasite centroids nevertheless had equally or more similar neighbours carrying different low-level labels. This is consistent with the limited resolving power of the short 18S marker and argues against upgrading an association to a species-specific host claim from sequence similarity alone.

## External biological validation

The most informative validation therefore comes from **lineage-specific host biology**, not forced placement of a ~132-nt fragment. Kim et al. (2026) place the type species *Amoebophrya sticholonchae*, directly observed infecting the radiolarian *Sticholonche zanclea*, within MALV-II clade 7. This provides independent same-clade host-class precedent for the Group-II-clade-7/radiolarian association in our dataset. It does not demonstrate that our associated radiolarian taxon is infected.

Current phylogenomics also separates MALV-I from MALV-II/IV: Holt et al. (2023) retain Syndiniales for MALV-II/IV and place MALV-I with Ichthyodinida. Consequently, *Amoebophrya* host records from Group II must not be transferred directly to Group I associations.

## Interpretation

Phase 6 supports three statements:

1. every priority sequence-cluster target is traceable to the exact abundance feature analyzed;
2. the short centroids support operational higher-level clade assignments but not confident species/strain-level host identification;
3. candidate-host interpretation should combine association statistics with lineage-specific independent biology.

References: Kim et al. 2026, *Journal of Phycology*, `10.1111/jpy.70160`; Holt et al. 2023, *Nature Communications* 14:7049, `10.1038/s41467-023-42807-0`.
