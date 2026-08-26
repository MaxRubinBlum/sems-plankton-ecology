# Targeted ASV-level validation — Phase 5

Phase 5 tests whether high-priority taxon-aggregate parasite–protist associations are supported by specific original ASV pairs rather than being created only by taxonomic aggregation.

## Procedure

1. Map each candidate aggregate feature back to its constituent 18S OTUIDs using the same aggregation logic used for FlashWeave input.
2. Within the 249 complete environmental cases, retain constituent ASVs with prevalence >=2% and >=50 total reads.
3. For broad aggregates, retain ASVs accounting for approximately 90% of reads, capped at 20 ASVs per aggregate.
4. Residualize log1p ASV abundances against the exact metadata matrix supplied to FlashWeave.
5. Test parasite-ASV × partner-ASV residual Spearman correlations where co-presence >=20 samples.
6. Apply one BH-FDR correction across all targeted ASV-pair tests; support requires q<0.05 and rho>=0.20.

## Archived result

102 constituent ASV pairs were testable. 46 retained positive support after global FDR. Seven of nine aggregate candidates with sufficient constituent-ASV co-presence contained at least one supported specific ASV pair.

This is a targeted sensitivity analysis, not an exhaustive search of all possible ASV pairs.
