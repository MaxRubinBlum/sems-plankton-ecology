# SEMS plankton ecology

Data, code and analytical provenance for the manuscript **Coordinated microbial community turnover across an ultra-oligotrophic water column**.

The study combines long-term 16S and 18S rRNA gene amplicon surveys with hydrographic and chemical observations from the southeastern Mediterranean Sea. The repository supports the analyses of vertical community transitions, trophic representation, cross-domain concordance and MALV parasite associations.

## Start here

1. Read [README_DATA.md](README_DATA.md) for the data inventory.
2. Install the Python dependencies:

   ```bash
   python -m pip install -r reproducibility/requirements.txt
   ```

3. Read [reproducibility/README.md](reproducibility/README.md) for the analysis map and figure-specific commands.
4. Compare regenerated results with the frozen checkpoints in [reproducibility/VALIDATION.md](reproducibility/VALIDATION.md).

Figure 5 and parts of the MALV workflow additionally require Julia and FlashWeave; their environment files are stored beside the relevant scripts.

## Repository structure

- `data/raw/` — cleaned source abundance and taxonomy tables used by the analyses
- `data/processed/` — integrated metadata and frozen derived inputs
- `data/figure_inputs/` — compact figure-level inputs
- `reproducibility/` — authoritative figure workflows and validation checkpoints
- `results/` — versioned numerical outputs and selected figure products
- `methods/` — detailed analytical methods
- `docs/` — provenance notes, taxonomy decisions and archived analytical checkpoints
- `scripts/` — supporting data-processing and comparison utilities

## Core depth grouping

The `ds3` metadata field defines the primary ecological depth classes:

1. `A_surface`
2. `B_nearsurface`
3. `C_DCM`
4. `D_below_DCM`
5. `F_300-600`
6. `G_below_600` (including near-bottom samples)

## Taxonomy and feature terminology

The repository retains the original SILVA-based prokaryotic taxonomy and a GTDB release 226 crosswalk for comparison and traceability. The 18S taxonomy follows PR2.

The amplicon features used here are **99% de-novo OTU clusters represented by centroid sequences**, not unclustered ASVs. See [docs/TAXONOMY_CROSSWALK.md](docs/TAXONOMY_CROSSWALK.md) and [reproducibility/syndiniales/INPUT_MANIFEST.md](reproducibility/syndiniales/INPUT_MANIFEST.md).

## Reuse and citation

Code is released under the [MIT License](LICENSE). Cite the associated manuscript when using the analyses or data; full publication details and accession identifiers will be added after acceptance. Marker-gene read proportions are compositional proxies and should not be interpreted as cell abundance, biomass or process rates.
