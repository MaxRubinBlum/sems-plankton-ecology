# SEMS plankton ecology

Reproducible analysis of long-term 16S and 18S rRNA amplicon datasets from Eastern Mediterranean Sea water-column transects.

The project is organized around vertical ecosystem structure across prokaryotes, microbial eukaryotes and metazoan 18S signatures, with additional analyses of cross-domain ecological associations and depth-specific taxonomic succession.

## Repository layout

- `data/raw/` — original input tables and metadata; never edited in place
- `data/processed/` — reproducibly filtered/derived tables
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

## Taxonomy

Two prokaryotic taxonomic assignments are retained for comparison and traceability:

- the original SILVA-based taxonomy
- GTDB release 226 taxonomy generated with the QIIME 2 sklearn classifier

The GTDB workflow is documented in `docs/taxonomy.md`.

18S taxonomy and the original metazoan-containing 18S table are retained so that microbial-eukaryote and metazoan signals can be analyzed separately.

## Data policy

Raw tables are immutable. Cleaning and taxonomic filtering are performed only by scripts and written to `data/processed/`.

Repository is currently private while analyses and manuscripts are in development.
