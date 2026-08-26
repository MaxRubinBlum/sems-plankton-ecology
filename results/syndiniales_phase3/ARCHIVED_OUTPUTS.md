# Syndiniales Phase 3 archived outputs

Phase 3 evaluates exact partner-set overlap across the 220-m boundary using the Phase-2 cross-regime residual-supported edge table and a strict parasite-specific candidate-universe null.

## Input used for archived run

- `cross_regime_edge_retesting.csv`
- SHA-256: `d46a8058ad611bd387a8647a8ba3a06b13a5ea3b4162d8cf4396748d7d31b1fc`

## Archived output checksums

The archived analysis generated the following tables:

- `syndiniales_lineage_partner_turnover.csv` — SHA-256 `cbf172a8bb454d6c92e46c083da1c7dda70cfef361e9160b585c30ff4eaaff4f`
- `shared_partner_edges_and_sign_behavior.csv` — SHA-256 `8adcb7a5aadaae6188ec99cbc3ed305d52b84051f47cd4863a7caff86f55ce0e`
- `partner_turnover_candidate_universe_null_tests.csv` — SHA-256 `a61752ec8c3d89b302c4a04533e151cb56b7276bd521b128c452fc6256b6ca44`
- `global_partner_turnover_candidate_universe_null_test.csv` — SHA-256 `ef777190a513368a838504a846c09dfe75ad88e8bd38102f2a06e326ec43788e`
- `summary.csv` — SHA-256 `9fce86b692a3f4d875a3fc30a2aab1bc20f7d765b3667cce53329fa16d75539e`

## Validated archived result

- 241 residual-supported associations at `<=220 m`
- 140 residual-supported associations at `>220 m`
- 95 supported in both regimes
- 93/95 shared associations retained the same residual sign
- 45 parasite lineages were comparable across regimes
- observed shared partner edges = 95
- strict-null expectation = 97.12655, 95% interval 92–102
- lower-tail P for excess turnover = 0.273736
- upper-tail P for excess conservation = 0.835808
- 0/45 lineages showed significant excess turnover after BH-FDR
- 0/45 lineages showed significant excess conservation after BH-FDR

The intended interpretation is **association-strength restructuring within a constrained partner repertoire**, not wholesale partner replacement.
