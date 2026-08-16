# Protist functional-trait framework

## Aim

Functional annotation of Eastern Mediterranean 18S taxa follows an evidence-first workflow. Functional guilds are **not assigned directly from PR2 taxonomy**. Instead, taxa receive published biological traits and guilds are derived later from multivariate trait space.

## Annotation hierarchy

1. **Ramond et al. trait backbone** — match named taxa to the published marine-protist trait dataset used as the basis of Avrahami et al. (2025).
2. **External evidence extension** — for ecologically important taxa absent from Ramond, add only traits supported by peer-reviewed literature or established curated resources.
3. **Unresolved traits remain NA** — depth occurrence, co-occurrence with bacteria, or taxonomic intuition are not used to infer trophic traits.

Each extension records the taxon/rank to which the evidence applies, trait state, reference/DOI, confidence, and an explicit `do_not_infer` boundary.

## Deep-ocean extension

The deep-sea extension prioritizes taxa by mean abundance below 300 m and depth enrichment. Major targets include Syndiniales/MALV, diplonemids and Eupelagonemidae, Radiolaria/Acantharia, Labyrinthulomycetes, MAST, Picozoa and other abundant dark-ocean protists.

Evidence is propagated only to the taxonomic rank supported by a source. For example, a trait demonstrated for Eupelagonemidae can be applied at family level, whereas photosymbiosis demonstrated for a particular acantharian subclade is not generalized to all Acantharia.

## Traits versus ecological covariates

Depth affinity and 16S–18S association strength are **not** inputs to the trait clustering. They are response/context variables used later to test whether functional assemblages change vertically or exhibit cross-domain coupling. This avoids circular inference.

## Planned analysis

1. Collapse ASVs to the most reliable available taxonomic rank.
2. Match to the Ramond backbone.
3. Merge literature-curated deep-sea trait extensions.
4. Report trait-annotation coverage by sample and depth layer.
5. Calculate Gower distances using curated traits.
6. Ordinate functional space and evaluate cluster solutions/stability.
7. Name clusters only after inspecting their trait composition.
8. Test functional composition across depth/hydrography and integrate with 16S–18S coupling analyses.

## Provenance

The files under `data/traits/` are version-controlled analytical inputs. Every manual biological annotation must be traceable to a reference. The exploratory rule-based guild assignments generated during initial analysis are deliberately excluded from this workflow.

### Key methodological references

- Ramond P. et al. (2019). A trait-based approach for marine protists. *Environmental Microbiology*.
- Avrahami Y. et al. (2025). Seasonal Transition in the Dominance of Photoautotrophic and Heterotrophic Protists in the Photic Layer of a Subtropical Marine Ecosystem. *Environmental Microbiology Reports* 17:e70126.
