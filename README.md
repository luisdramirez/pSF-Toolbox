# Population spatial frequency Toolbox

The pSF Toolbox streamlines the population spatial frequency (pSF) tuning  approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

We provide a suite of scripts that estimate and simulate pSF (`/estimate-pSF`), and for stimulus presentation via [Psychtoolbox](http://psychtoolbox.org) when measuring pSF (`/measure-pSF`). 

We include an example workflow (`/estimate-pSF/example_pipeline.m`) for estimating pSF from a sample dataset (SF input and measured BOLD time series for two subjects — 100 voxels per ROI (V1–V3)).

The sample data is organized as a structure array with fields `I` and `measured_BOLD`. **Note that the shape of the time series data must have time along the first dimension (e.g., time point x voxel)!**

## estimatePSF.m
This is the main high-level function for estimating pSF parameters. It takes the measured BOLD time series, the stimulus spatial frequency time series, and a hemodynamic impulse response function (HIRF) as input to return a structure `pSF` containing:
- estimated pSF parameters (peak SF, bandwidth, BOLD amplitude, baseline)
- estimated pSF curves
- estimated neural time series
- estimated BOLD time series
- $R^2$ values
- SSE values
- `fmincon` exit flags

Toggles and settings to make note of that should be defined before entering this function (see `example_pipeline.m`):
- Parallelization (true/false)
- Coarse grid search (true/false)
- Fine grid search (true/false)
- Spatial frequency range used to generate tuning curves 
- Initial pSF parameters
- pSF parameter bounds

### fitVoxels.m
This function performs the voxel-wise pSF parameter estimation using `fmincon`.

## example_pipeline.m
An example workflow.

## Utility Functions
*   `logGauss.m`: Defines the log Gaussian function used to model the population spatial frequency tuning curve.
*   `calcFit.m`: Computes goodness-of-fit statistics, R-squared ($R^2$) and Sum of Squared Errors (SSE), to evaluate how well the estimated BOLD data — derived from the pSF model — explains the measured BOLD data.
*   `defineHRF.m`: Creates a canonical hemodynamic response function (HRF) model, which is used to convolve the predicted neural response to estimate the BOLD signal.
*   `gridSearch.m`: Implements a grid search algorithm to find suitable starting parameters for the more complex non-linear optimization performed in `fitVoxels.m`.
*   `chunkTimeSeries.m`: Splits large time series datasets into smaller, more manageable chunks for parallel processing.
*   `cpd2oct.m`: Calculates the pSF bandwidth in octave units from the half-maximum spatial frequencies in cycles per degree (cpd).


## simulate_pSF.m