# Figure 1 — hydrographic template

Figure 1 establishes the environmental structure in which the microbial samples were collected.

## Inputs

Primary inputs are the integrated sample/CTD tables for 16S and 18S samples and the station metadata. CTD variables used in the hydrographic ordination are temperature, salinity, dissolved oxygen and fluorescence. The chemistry-expanded ordination is retained as a supplementary analysis rather than replacing the CTD-only main panel.

## Analytical requirements

- Sampling depth in the vertical environmental panel is plotted on a linear depth scale, beginning at 0 m.
- Summer and winter samples use different symbols while retaining the established depth/habitat color scheme.
- The hydrographic PCA is calculated after centering and scaling the CTD variables temperature, salinity, dissolved oxygen and fluorescence so variables measured in different units contribute comparably.
- PCA is fitted only to complete observations for the variables entering that ordination; sample filtering and final n must be exported with the scores.
- The plotted PCA scores, variable loadings, explained variance and sample identifiers must be written to derived tables.
- Chemistry is not silently mixed into the main hydrographic PCA. The expanded CTD+chemistry analysis is a separate supplementary workflow with its own complete-case set and documented variables.

## Outputs to retain

- final Figure 1 SVG/PDF/600-dpi JPEG
- PCA scores
- PCA loadings
- explained variance
- sample-level environmental table used for plotting
- station/depth/season mapping used for aesthetics

## Methods traceability

The manuscript Methods should state that hydrographic PCA used standardized CTD temperature, salinity, dissolved oxygen and fluorescence associated with the microbial samples. The supplementary chemistry-expanded analysis must list its additional variables explicitly and report any reduction in sample number due to chemistry availability.