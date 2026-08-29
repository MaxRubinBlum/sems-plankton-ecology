# Parasite candidate-host evidence grading — Phase 4

Phase 4 separates **lineage-specific candidate-host associations** from host-plausible and mechanistically unresolved ecological associations. It does not infer infection from co-occurrence alone.

## Inputs

- Phase-2 cross-regime retesting table
- Phase-2 formal partner × regime interaction results
- Phase-2 non-metazoan CLR sensitivity results
- original 18S taxonomy ranks for partner-guild assignment
- published host biology coded as an explicit evidence class, not as a statistical prior

## Quantitative inclusion

A pair is considered positive in a regime when it has residual support from Phase 2, residual Spearman rho >= 0.20, and at least 20 co-present samples in that regime. The quantitative evidence score adds support for replication across regimes, FlashWeave selection, a significant regime interaction, and same-direction CLR sensitivity.

## Taxonomy correction

The original partner-name heuristic was audited against the source taxonomy table. The corrected workflow assigns broad partner guilds from phylum/class/order/family/genus ranks. This is important for labels such as radiolarian families that lack `RAD` in the displayed name, and for names such as *Gymnoxanthella radiolariae*, which is a dinoflagellate rather than a radiolarian.

## Lineage-specific host evidence

Host precedent is no longer transferred across all operational `Dino-Group` parasites.

- MALV-II/Dino-Group II and MALV-IV correspond to Syndiniales in current phylogenomics.
- MALV-I/Dino-Group I is phylogenetically distinct and has been proposed as Ichthyodinida; broad host-class precedent is allowed, but *Amoebophrya* host records from MALV-II are not treated as direct evidence for MALV-I clades.
- Group II clade 7 receives lineage-specific radiolarian host-class support because *Amoebophrya sticholonchae*, a parasite directly observed infecting the radiolarian *Sticholonche zanclea*, falls within MALV-II clade 7 (Kim et al. 2026). This supports the host class, not the specific identity of an associated radiolarian taxon in our data.
- Several Group II clades represented in the *A. ceratii* complex receive lineage-specific dinoflagellate host-class precedent.
- *Hematodinium*/*Syndinium* protist associations remain ecological rather than candidate final-host links because these lineages have established metazoan hosts.

## Corrected archived result

The corrected evidence table contains 228 positive associations:

- **6 Tier-1 lineage-specific candidate-host associations**
- 61 Tier-2 host-plausible associations
- 50 Tier-3 robust ecological associations
- 111 exploratory associations

Tier labels are evidence summaries, **not probabilities of infection**.

The strongest qualitative addition is the Group II clade 7–radiolarian relationship: our clade-7 association with *Acanthochiasma* strengthens below 220 m, while independent biological work documents a different radiolarian host for a clade-7 *Amoebophrya*. This combination makes the association unusually compelling but still does not demonstrate infection of *Acanthochiasma*.

## References

- Kim et al. 2026. *Journal of Phycology*. `10.1111/jpy.70160`.
- Holt et al. 2023. *Nature Communications* 14:7049. `10.1038/s41467-023-42807-0`.
- Anderson et al. 2024. *ISME Communications* 4, ycae014. `10.1093/ismeco/ycae014`.
- Anderson & Harvey 2020. *mSphere* 5:e00209-20.

## Interpretation guardrails

Candidate-host terminology is restricted to associations with quantitative support plus relevant **lineage-specific** biological precedent. Short environmental rRNA sequences and statistical associations cannot by themselves establish a host identity or infection state.
