# Population spatial frequency Toolbox

The pSF Toolbox streamlines the population spatial frequency tuning (pSFT) approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

We provide a suite of scripts that estimate and simulate pSFT, and generate the stimulus presentation used for pSFT mapping via [Psychtoolbox](http://psychtoolbox.org).

We also provide sample data and estimates for two subjects — a hundred voxels per ROI (V1–V3), as well as a script for how you might use our pipeline.

The sample data contains a structure with fields `measured_BOLD` [voxel×TR×ROI] and `I` (1xTR). 
Sample estimates contains an array [voxel x estimate x ROI x subject]:
- col 1: $R^2$
- col 2: pSFT peak
- col 3: pSFT bandwidth (cpd)
- col 4: BOLD amplitude
- col 5: BOLD baseline
- col 6: pRF eccentricity
- col 7: pRF polar angle
- col 8: pRF size

## estimatePSF.m
This is the main high-level function for estimating pSFT parameters. It takes the measured BOLD time series, the stimulus spatial frequency time series, and a hemodynamic impulse response function (HIRF) as input to return a structure, `pSFT`, containing:
- estimated pSFT parameters (peak SF, bandwidth, BOLD amplitude, baseline)
- estimated pSFT curves
- estimated neural time series
- estimated BOLD time series
- $R^2$ values
- SSE values
- `fmincon` exit flags

Toggles and settings to make note of that should be defined before entering this function (see `example_pipeline.m`):
- Parallelization
- Coarse and fine grid search
- Plot pSF and voxel time series (measured & estimate)
- Spatial frequency range used to generate tuning curves
- Initial pSF parameters
- pSF parameter bounds


### fitVoxels.m
This function performs the voxel-wise pSF parameter estimation using `fmincon`.

## example_pipeline.m
An example workflow: loading the sample data (post voxel selection), defining the HRF using `defineHRF.m`, calling `estimatePSF.m` to perform the analysis, and visualizing the results.

## Utility Functions
*   `chunkTimeSeries.m`: Splits large time series datasets into smaller, more manageable chunks for parallel processing or memory efficiency.
*   `cpd2oct.m`: Calculates the pSFT bandwidth in octave units from the half-maximum spatial frequencies in cycles per degree (cpd).
*   `gridSearch.m`: Implements a grid search algorithm to find suitable starting parameters for the more complex non-linear optimization performed in `fitVoxels.m`.
*   `calcFit.m`: Computes goodness-of-fit statistics, R-squared ($R^2$) and Sum of Squared Errors (SSE), to evaluate how well the estimated BOLD data — derived from the pSFT model — explains the measured BOLD data.
*   `defineHRF.m`: Creates a canonical hemodynamic response function (HRF) model, which is used to convolve the predicted neural response to estimate the BOLD signal.
*   `logGauss.m`: Defines the log Gaussian function used to model the population spatial frequency tuning curve.

## simulate_pSFT.m