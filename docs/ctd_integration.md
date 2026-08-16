# CTD integration checkpoint

HaiSec CTD casts are integrated with the 16S/18S molecular metadata to attach in-situ environmental conditions and a cast-specific fluorescence maximum (used here as a DCM proxy).

## Current matching

- 18S: 249/310 samples have successful CTD interpolation; 231 have DCM-relative depth.
- 16S: 512/789 samples have successful CTD interpolation.
- The principal unmatched block is 2026, outside the current HaiSec 35–54 CTD archive.

## Variables

For each matched sample the workflow interpolates temperature, salinity, dissolved oxygen and fluorescence at sampling depth. The DCM proxy is the maximum of a 5-point rolling-median fluorescence profile between 2 and 200 dbar (only casts reaching at least 30 dbar). `delta_depth_from_dcm_m = sample depth - DCM depth`.

## Preliminary ecological result

In summer, absolute depth remains slightly more strongly associated with the major 18S trophic axes than DCM-relative depth, indicating that vertical succession is not simply tracking the moving DCM. Nevertheless, the fluorescence maximum marks a strong trophic transition: phototrophs concentrate around it, Radiolaria/phagotrophy increase immediately below it, and diplonemid heterotrophy becomes prominent deeper in the water column.

After controlling for absolute depth, local fluorescence remains positively associated with bona fide phototrophy and mixotrophy. These are exploratory associations and should be followed by multivariable community/environment models before mechanistic interpretation.

## Next analyses

1. Whole-community 16S and 18S environmental partitioning using depth + CTD variables + season/geography.
2. Objective ecological-boundary analysis comparing 16S, 18S and functional turnover against hydrographic transitions.
3. Publication figure integrating bathymetry, CTD sections and molecular sampling positions.
