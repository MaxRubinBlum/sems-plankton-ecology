# SEMS plankton ecology

Reproducible analysis of long-term 16S and 18S rRNA amplicon datasets from Eastern Mediterranean Sea water-column transects.

The project is organized around vertical ecosystem structure across prokaryotes, microbial eukaryotes and metazoan 18S signatures, with additional analyses of cross-domain ecological associations and depth-specific taxonomic succession.

## Repository layout

- `data/raw/` — original input tables and metadata; never edited in place
- `data/processed/` — reproducibly filtered/derived tables
- `data/traits/` — evidence-backed protist trait annotations and curation inputs
- `scripts/` — analysis and data-processing code
- `results/` — derived statistics, tables and figures
- `docs/` — methods/provenance notes

## Core depth grouping

The primary ecological depth grouping is the `ds3` field already present in the metadata:

1. `A_surface`
2. `B_nearsurface`
3. `C_DCM`
4. `D_below_DCM`
5. `F_300-600`
6. `G_below_600` (includes near-bottom samples)

Where analyses require an explicit near-bottom category, scripts may use the seven-layer `ds2` field, including `H_near_bottom`; the selected field must be reported with each result.

## Taxonomy

Two prokaryotic taxonomic assignments are retained for comparison and traceability:

- the original SILVA-based taxonomy
- GTDB release 226 taxonomy generated with the QIIME 2 sklearn classifier

The GTDB workflow is documented in `docs/taxonomy.md`.

18S taxonomy and the original metazoan-containing 18S table are retained so that microbial-eukaryote and metazoan signals can be analyzed separately.

## Protist functional traits

Functional ecology is built from published traits rather than direct taxonomy-to-guild rules. Ramond et al. provides the backbone used by Avrahami et al. (2025); abundant dark-ocean taxa missing from that framework are extended only where peer-reviewed evidence supports a trait at the relevant taxonomic rank. Missing traits remain `NA`.

Depth distribution and 16S–18S co-occurrence are not used to assign function, preventing circular inference. Functional groups are named only after multivariate trait-space analysis. See `docs/protist_traits.md`.

## Data policy

Raw tables are immutable. Cleaning and taxonomic filtering are performed only by scripts and written to `data/processed/`.

Repository is currently private while analyses and manuscripts are in development.
