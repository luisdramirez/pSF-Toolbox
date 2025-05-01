# Population spatial frequency Toolbox

The pSF Toolbox streamlines the population spatial frequency (pSF) tuning approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

We provide a suite of scripts that estimate and simulate pSF (`/estimate-pSF`), and for stimulus presentation via [Psychtoolbox](http://psychtoolbox.org) when measuring pSF (`/measure-pSF`). 

We include an example workflow (`/estimate-pSF/example_pipeline.m`) for estimating pSF from a sample dataset (SF input and measured BOLD time series for two subjects — 100 voxels per ROI (V1–V3)).

`sample_data` is organized as a structure array with fields `I` and `measured_BOLD`. 

**!! Note that the shape of the time series data must have time along the first dimension (e.g., time point x voxel) !!**

### `/estimate-pSF`
This directory contains scripts for estimating pSF parameters from fMRI data.
-   `estimatePSF.m`: Main function for pSF parameter estimation. Takes measured BOLD, stimulus SF time series, and HIRF as input. Returns estimated parameters, curves, time series, and goodness-of-fit metrics. See more info in dedicated section below. 
-   `example_pipeline.m`: Demonstrates a complete workflow for estimating pSF parameters using sample data. Includes setting up estimation parameters (parallelization, grid search, parameter bounds, HRF definition) and visualizing results.
-   `functions/`: Contains core functions used by the estimation scripts:
    -   `fitVoxels`: Performs voxel-wise parameter estimation using `fmincon`.
    -   `logGauss`: Defines the log Gaussian function for the pSF tuning curve.
    -   `calcFit`: Computes SSE.
    -   `defineHRF`: Creates a canonical HRF model.
    -   `gridSearch`: Implements grid search for initial parameter estimates.
    -   `chunkTimeSeries`: Splits time series for parallel processing.
    -   `cpd2oct`: Converts bandwidth from cpd to octaves.
-   `simulate_pSF.m`: Useful for generating synthetic tuning curves.

### `/measure-pSF`
This directory provides scripts for presenting stimuli using Psychtoolbox to measure pSF tuning.
-   `run_scan.m`: Main script for running the pSF experiment. Handles stimulus presentation, timing, and response collection (if applicable). Requires configuration based on experimental setup.
-   `functions/`: Contains supporting functions for stimulus generation, display, and experimental control.
-   `stimuli/`: Likely contains stimulus definitions or generation scripts (e.g., grating parameters).
-   `data/`: Intended for storing acquired experimental data.

### estimatePSF
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
