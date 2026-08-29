# Hydrographic-coordinate analysis of seasonal ecological boundaries

## Objective

To test whether seasonal shifts in ecological transition depths reflected changes in the physical position of the water-column structure, biological breakpoint distributions were compared in raw depth space and after mapping them into several hydrographic coordinates.

## Hydrographic data

Full-resolution HaiSec CTD cast data were used rather than only environmental values interpolated to individual metabarcoding samples. CTD variables included pressure, temperature, salinity, dissolved oxygen and fluorescence. Casts were linked to the cruises and seasons represented in the 16S and 18S datasets. Analyses of the upper-water-column transitions were restricted to approximately 10–420 dbar, matching the range used in the ecological breakpoint scans.

Each CTD cast was interpolated to a common 5-dbar grid. Potential temperature referenced to the surface was calculated using the UNESCO 1983 adiabatic-temperature-gradient formulation, and potential-density anomaly (sigma-theta) was calculated from EOS-80 density evaluated at the potential temperature. This sigma-theta estimate was used as a density coordinate for comparison among seasons; it is not a TEOS-10 Absolute Salinity/Conservative Temperature calculation.

## Multivariate hydrographic coordinate

Temperature, salinity, dissolved oxygen and fluorescence were standardized to zero mean and unit variance across the relevant CTD profiles. Principal-component analysis was then applied to these standardized variables. Hydrographic PC1 was retained as a multivariate coordinate describing the dominant combined physical-biogeochemical gradient.

## Biological breakpoint distributions

Biological transition depths were taken from the existing profile-aware threshold analyses for:

- 16S community composition;
- 18S community composition; and
- 18S trophic structure.

For each dataset and season, the threshold scan identified the candidate depth producing the largest partial R-squared after accounting for sampling profile. Profile-bootstrap distributions of the optimal threshold were used to quantify breakpoint uncertainty and seasonal displacement. Existing threshold-scan and bootstrap procedures were retained unchanged; the present analysis only re-expressed those inferred breakpoint depths in alternative hydrographic coordinates.

Syndiniales were evaluated separately in the parasite-focused analysis and were not included in the numerical comparison among hydrographic coordinates reported here.

## Mapping breakpoints into hydrographic space

For every bootstrap-derived breakpoint depth, the corresponding hydrographic state was obtained from the full-resolution CTD profiles represented by that biological dataset and season. Median values across applicable casts were calculated for:

- pressure/depth;
- salinity;
- sigma-theta;
- dissolved oxygen; and
- Hydrographic PC1.

Thus, the same biological bootstrap distribution was expressed in each candidate environmental coordinate without re-optimizing the biological breakpoint separately for each variable.

## Comparison of seasonal alignment

For each biological dataset and hydrographic coordinate, winter and summer bootstrap distributions were compared using three descriptive measures:

1. **Standardized seasonal separation**: the absolute difference between winter and summer bootstrap means divided by the pooled standard deviation. Smaller values indicate closer seasonal alignment.
2. **Normalized Wasserstein distance**: the Wasserstein distance between winter and summer bootstrap distributions divided by their pooled standard deviation.
3. **95% interval overlap**: the overlap of the winter and summer 2.5–97.5 percentile intervals, expressed relative to the total span of both intervals.

The standardized seasonal separation was used for the summary comparison in the supplementary figure; the other metrics were retained as sensitivity descriptors.

## Interpretation

This analysis tests whether seasonal movement of an ecological transition in meters becomes smaller when the transition is expressed in a hydrographic coordinate. A reduction in seasonal separation indicates that the biological transition tracks that component of the hydrographic structure more closely than a fixed depth. It does not by itself demonstrate a causal effect of that variable or identify a discrete named water mass.

## Output files

Results are stored in `results/hydrography/`:

- `biological_boundaries_in_hydrographic_coordinates.csv`
- `hydrographic_coordinate_seasonal_alignment.csv`
- `hydrographic_coordinate_bootstrap_intervals.csv`
