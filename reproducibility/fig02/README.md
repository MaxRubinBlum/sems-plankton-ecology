# Figure 2 — objective ecological boundary detection

Figure 2 tests whether community turnover identifies a reproducible deep ecological boundary and whether its position changes seasonally.

## Core principle

Boundary depth is estimated from the community data rather than imposed from predefined depth bins. Global and season-specific analyses are retained. The breakpoint analysis and boundary scan are complementary views of the same ecological transition and must use documented, reproducible candidate-depth rules.

## Inputs

- cleaned 16S and 18S community tables
- sample metadata containing depth, station, season, bottom depth and near-bottom status
- 18S-derived trophic-trait matrix for the trait boundary analysis
- CTD profiles for the supplementary DCM-depth analysis; DCM is derived from profiles rather than sample-category labels and profiles with dual fluorescence maxima must retain both peaks in the source-derived DCM table

## Boundary scan

For each candidate depth, samples/profiles are separated into shallower and deeper sets and ecological dissimilarity/turnover across the candidate boundary is quantified. Candidate depths must satisfy the support criterion used in the finalized analysis (at least 15 independent profiles contributing support on the relevant sides of the candidate boundary). The exact score definition, candidate grid and distance transformation are to be encoded in the analysis script rather than inferred from the figure.

Global scans are followed by winter- and summer-specific scans using the same algorithm. The seasonal analyses are not independently tuned to maximize visual separation.

## Statistical inference

The finalized analysis uses profile-level resampling rather than treating individual depths from the same station profile as independent replicates. Uncertainty in boundary position is estimated by profile bootstrap. Significance of the maximum boundary statistic is assessed with a max-statistic permutation procedure so inference accounts for scanning across multiple candidate depths. Random seeds and number of bootstrap/permutation iterations must be fixed and reported by the script.

## Covariate-controlled test

The deep-boundary effect is additionally tested while controlling for station, station bottom depth and near-bottom status. This analysis addresses whether the detected transition is simply produced by geography/bathymetry or by sampling close to the seabed. Model formula, response/distance representation, permutation restrictions and sample counts must be exported with the results.

## Seasonal DCM context

Seasonal DCM position is derived from CTD fluorescence profiles and belongs in the Supplement rather than the main boundary figure. Dual DCM peaks are retained where detected. The DCM distributions may be shown as transparent seasonal bands in boundary plots for context, but DCM depth is not used to define the ecological boundary.

## Required outputs

- candidate-depth boundary scores for 16S, 18S and trophic traits
- global and winter/summer scans
- bootstrap distributions and confidence intervals for peak depth
- max-statistic permutation null distributions and p-values
- breakpoint estimates and uncertainty
- covariate-controlled boundary-test table
- CTD-derived DCM peak table and seasonal summary
- exact sample/profile counts supporting each candidate depth

## Methods traceability

The manuscript must clearly distinguish (i) objective community boundary detection, (ii) breakpoint estimation, (iii) covariate-controlled confirmation and (iv) independent hydrographic/DCM context. This prevents the ecological boundary from being presented as a visually chosen depth threshold.