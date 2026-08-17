# Figure 1–2 publication checkpoint

## Figure 1

Current publication asset: `Figure1_publication_with_loadings`.

Architecture:
- sampling depth and station bathymetric context on a linear depth axis;
- hydrographic PCA based on CTD temperature, salinity, dissolved oxygen and fluorescence;
- hydrographic PCA includes loading vectors for Temperature, Salinity, O₂ and Fluorescence;
- 16S and 18S Bray–Curtis PCoAs;
- depth layer encoded by color and season by symbol (winter circle; summer triangle);
- no figure header or panel headers/letters in the final artwork, for manuscript layout editing.

Hydrographic PCA: PC1 = 57.1%, PC2 = 25.5%.

The manuscript legend is maintained in `docs/figure_legends_fig1_fig2.md`.

## Figure 2

Current publication asset: `Figure2_publication`.

Architecture:
- pooled biological boundary scan;
- pooled profile-bootstrap breakpoint uncertainty;
- winter/summer boundary scans;
- winter/summer profile-bootstrap breakpoint uncertainty;
- winter and summer DCM depth distributions are shown as transparent CTD-derived density fields behind the boundary-scan panels only;
- DCM peaks are extracted directly from CTD fluorescence profiles and may include a second peak when supported by peak prominence;
- no figure header or panel headers/letters in the final artwork.

Pooled strongest breakpoints:
- 18S trophic/ecological traits: 120 m;
- 16S composition: 170 m;
- 18S composition: 220 m.

Season-specific strongest breakpoints:
- winter: traits 120 m, 16S 220 m, 18S 220 m;
- summer: traits 120 m, 16S 140 m, 18S 160 m.

All three pooled max-statistic permutation tests: P = 0.002.

Figure inputs committed under `data/figure_inputs/` include the breakpoint summary and CTD-derived winter/summer DCM peaks.

## Supplementary environmental PCA

Supplementary Figure X is the chemistry-matched hydrographic/biogeochemical PCA.

Complete-case sample sizes: 45 18S samples and 48 16S samples. The 18S visualization uses standardized temperature, salinity, dissolved oxygen, fluorescence, NO₃+NO₂, PO₄, silicate, pH and total alkalinity. PC1 = 63.0% and PC2 = 22.6% for the 18S chemistry-matched subset.

The supplementary legend is maintained in `docs/supplementary_figure_legends.md`, and PCA loadings are committed under `data/figure_inputs/chemistry_hydrography_PCA_loadings.csv`.
