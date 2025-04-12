# Population spatial frequency Toolbox

The pSF Toolbox streamlines the population spatial frequency tuning (pSFT) approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

We provide code for estimating pSFT, simulating pSFT, and generating the stimulus presentation used for pSFT mapping via [Psychtoolbox](http://psychtoolbox.org).

We also provide sample data and estimates for two subjects — a hundred voxels per ROI (V1–V3), as well as a script for how you might use our pipeline.

The example data contains a structure, `example_data`, with fields `measured_BOLD` [voxel×TR×ROI] and `I` (1xTR); `example_estimates` contains an array [voxel x estimate x ROI x subject]:
- col 1: $R^2$
- col 2: pSFT peak
- col 3: pSFT bandwidth (cpd)
- col 4: BOLD amplitude
- col 5: BOLD baseline
- col 6: pRF eccentricity
- col 7: pRF polar angle
- col 8: pRF size

## Fit_pSFT.m
This is the main function that estimates pSFT from a measured BOLD time series, input spatial frequeny time series, and hemodynamic impluse response function. As long as the inputs are in the required format (e.g., the shape of the measured BOLD time series matrix has voxels along the rows and time points along the columns), you will receive a structure, `pSFT`, that contains:
- estimated pSFT parameters 
- estimated pSFT curves
- estimated neural time series
- estimated BOLD time series
- estimated $R^2$ values
- estimated SSE values
- fmincon exit flags



### fitVoxels.m
Here you'll find where the magic happens


## simulate_pSFT.m