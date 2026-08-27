# Sequence-level validation targets — Phase 6

Phase 6 now has the matching 18S representative-sequence artifact and completes the sequence-extraction/provenance step for high-priority candidate host associations.

## Representative-sequence input

Uploaded artifact: `rep-seqs-99.qza`

SHA-256: `751d987161559e47a7c39774567e4c966817741c331d7f58ec2803cb5c656670`

QIIME 2 type: `FeatureData[Sequence]`

The artifact contains 40,736 representative sequences. Its feature-ID set exactly matches the 40,736 features in the current 18S OTU table.

## Provenance

QIIME 2 provenance shows two successive steps:

1. `q2-dada2 denoise_paired`, with hashed feature IDs.
2. `q2-vsearch cluster_features_de_novo`, with `perc_identity = 0.99`.

This is analytically important: the current table consists of **99% de-novo OTU clusters represented by centroid sequences**. Although the centroid IDs originate from hashed DADA2 feature IDs, the abundance features analyzed after clustering must not be described as unclustered ASVs.

This provenance correction also applies to Phase 5. Its numerical result is unchanged, but the targeted deaggregation is correctly described as **99%-OTU-level validation**.

## Target extraction

The Phase-6 priority list contains 78 requested OTUIDs derived from the highest-priority candidate associations.

- requested IDs: 78
- recovered from `rep-seqs-99.qza`: 78
- missing: 0
- unique recovered sequences: 78
- sequence lengths: 130–136 nt; median 132 nt

`extract_repseq_targets.py` validates the artifact against the OTU table, checks provenance, extracts the requested centroid sequences, and writes a taxonomy-linked FASTA and metadata table.

The archived target FASTA is intentionally small and contains only the 78 priority centroid sequences rather than the complete 40,736-sequence artifact.

## What Phase 6 establishes

The sequence file now confirms that every sequence-cluster target used for follow-up is recoverable and traceable to the exact abundance feature analyzed in Phases 1–5. It also prevents overstatement of the Phase-5 resolution: support is at the 99%-OTU-centroid level, not the unclustered-ASV level.

## External reference matching

The next sequence step is external reference matching/phylogenetic placement of the priority Syndiniales centroids against curated PR2/GenBank references, especially characterized Group I/II, *Amoebophrya*, and *Euduboscquella* sequences. This step should report alignment coverage and identity explicitly because the available V4-region centroids are short (~132 nt), and species-level assignments from short exact/near-exact matches require caution.

No infection claim should be upgraded solely from sequence identity: host interpretation still requires the independent biological evidence framework defined in Phase 4.
