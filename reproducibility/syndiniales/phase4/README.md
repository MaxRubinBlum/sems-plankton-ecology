# Syndiniales candidate-host evidence grading — Phase 4

Phase 4 separates **host-compatible candidate associations** from robust but mechanistically unresolved ecological associations. It does not infer infection from co-occurrence alone.

## Inputs

- Phase-2 cross-regime retesting table (`cross_regime_edge_retesting.csv`)
- Phase-2 formal partner × regime interaction results
- Phase-2 non-metazoan CLR sensitivity results
- Published host biology used only as an explicit evidence class, not as a statistical prior

## Quantitative inclusion

A pair is considered positive in a regime when it has residual support from Phase 2, residual Spearman rho >= 0.20, and at least 20 co-present samples in that regime.

The transparent evidence score adds support for replication across regimes, FlashWeave selection, a significant regime interaction, and same-direction CLR sensitivity. Literature evidence is coded separately.

## Literature classes

- `documented-host-genus`: partner genus has direct Syndiniales infection precedent (e.g. *Prorocentrum*, *Gymnodinium*, *Gyrodinium*).
- `documented-host-group`: the broader taxonomic group contains verified Syndiniales hosts.
- `cross-study`: the same Syndiniales clade and candidate host class recurs independently in the BATS depth-network study (Anderson et al. 2024).
- `emerging-host-group`: repeated sequence/network evidence exists, but direct infection is unconfirmed (notably radiolarians).
- `ecological`: known parasite biology argues against interpreting the protist partner as the final host (e.g. *Hematodinium*/*Syndinium* protist links).

## Archived result

The refined table contains 228 positive candidate associations: 8 Tier-1 strong candidate-host associations, 31 Tier-2 host-plausible associations, 84 Tier-3 robust ecological associations, and 105 exploratory associations.

Tier labels are evidence summaries, **not probabilities of infection**.

## References used for host-precedent coding

- Anderson SR et al. 2024. ISME Communications 4, ycae014. doi:10.1093/ismeco/ycae014
- Anderson SR & Harvey EL. 2020. mSphere 5, e00209-20.
- Bråte J et al. 2012. Protist 163:767–777. doi:10.1016/j.protis.2012.04.004
- Coats DW et al. 2012. J Eukaryot Microbiol 59:1–11. doi:10.1111/j.1550-7408.2011.00588.x
- Maranda L. 2001. J Phycol 37:245–248. doi:10.1046/j.1529-8817.2001.037002245.x

## Interpretation guardrail

Candidate-host terminology is restricted to associations with quantitative support plus relevant biological precedent. Radiolarian associations remain candidate-host hypotheses because direct Syndiniales infection of radiolarians has not been demonstrated.
