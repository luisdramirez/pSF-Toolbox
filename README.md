# The Population Spatial Frequency Toolbox

The pSF-Toolbox streamlines the population spatial frequency tuning (pSFT) approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

**Requirements**
- Psychtoolbox-3 must be installed for stimulus presentation.
- The Optimization Toolbox for MATLAB must be installed for parameter optimization.
- The shape of the BOLD percent signal change time series data must have time along the first dimension and voxels along the second (i.e., time points x voxels).
- **Strongly recommended**: The Parallel Computing Toolbox for MATLAB must be installed for parallelization.

We provide a suite of scripts for (1) stimulus presentation via Psychtoolbox-3 to measure pSFT with fMRI (see `/measure-pSF`) and (2) voxel-wise parameter optimization (see `/estimate-pSF`). 
 
### `/measure-pSF`
This directory contains scripts for executing the experiment via Psychtoolbox.
We provide an example scan session script for data acquisition (see `/measure-pSF/run_session.m`) that can be modified with respect to the experimental setup. For example, the input device name, toggles (e.g., save run info), subject ID, directories, and screen parameters should be verified by the user. 

Critical functions include `prepareScan` and `presentStimuli`. 
Users will find key stimulus and timing parameters inside `prepareScan`. For example, to adjust the size of the stimulus, the user must change `p.aperture_radius_deg` (or `p.aperture_radius_px`); to match the fMRI scan length, `t.TR` must match the duration of the repetition time. 

`presentStimuli` will output a structure that compiles all the experiment structures (scan parameters `p`, timing parameters `t`, window paramters `w`, frame sequences `frames`, and behavioral data `behav_data`).  

**Directory contents**
-   `/stimuli`: Stimulus textures will be stored here by default.
-   `/data`: Experimental run info will be stored here by default.
-   `/functions`: Contains supporting functions for stimulus generation, display, and experimental control.
    -   `checkPTB`: Verifies Psychtoolbox installation.
    -   `prepareScan`: Initializes parameters, stimuli, timing, and Psychtoolbox window for a scan.
    -   `createTextures`: Generates bandpass-filtered noise textures for stimuli.
    -   `createApertures`: Creates stimulus apertures.
    -   `genFrames`: Defines the sequence of events and timing for each experimental frame.
    -   `presentStimuli`: Draws stimuli frame by frame. Compiles run information (e.g., parameters, behavioral data) into struct `run_info`.


### `/estimate-pSF`
This directory contains scripts for estimating pSFT parameters from fMRI data.
We include an example workflow for estimating pSFT from a sample dataset that contains SF input and measured BOLD time series from two subjects — 100 voxels in V1, V2, and V3 (see `/estimate-pSF/example_pipeline.m`). `sample_data` is a structure array with fields `I` and `measured_BOLD`.

Below are toggles and parameters that must be defined before entering `estimatePSF.m` (see `example_pipeline.m`).

Toggles:
- Parallelization (true/false)
- Coarse grid search (true/false)
- Fine grid search (true/false)

Parameters:
- Spatial frequencies used to generate tuning curves (`p.sfs`)
- Initial pSFT parameters (`p.init_params`)
- pSFT parameter bounds (`p.pSFT_bounds`)


**Directory contents**
-   `example_pipeline`: Demonstrates a complete workflow for estimating pSFT parameters using sample data. Includes setting up estimation parameters (parallelization, grid search, parameter bounds, HRF definition) and visualizing results.
-   `estimatePSF`: This is the main high-level function for estimating pSFT parameters. It takes the measured BOLD time series, the stimulus spatial frequency time series, and a hemodynamic impulse response function (HIRF) as input to return a structure `pSFT` containing:
    - estimated pSFT parameters (peak SF, bandwidth, BOLD amplitude, baseline)
    - estimated pSFT curves
    - estimated neural time series
    - estimated BOLD time series
    - $R^2$ values
    - SSE values
    - `fmincon` exit flags
-   `/functions`: Contains supporting functions:
    -   `fitVoxels`: Performs voxel-wise parameter estimation using `fmincon`.
    -   `logGauss`: Defines the log Gaussian function for the pSFT tuning curve.
    -   `calcFit`: Computes SSE.
    -   `defineHRF`: Creates a canonical HRF model.
    -   `gridSearch`: Implements grid search for initial parameter estimates.
    -   `chunkTimeSeries`: Splits time series for parallel processing.
    -   `cpd2oct`: Converts bandwidth from cpd to octaves.
    -   `checkRequiredToolboxes`: Verifies that the required MATLAB Toolboxes are installed.
-   `simulate_pSF.m`: Useful for generating synthetic tuning curves.



