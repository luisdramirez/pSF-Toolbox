# The Population Spatial Frequency Toolbox

The pSF-Toolbox streamlines the population spatial frequency tuning (pSFT) approach, originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019), with accessibility in mind. 

We provide a suite of scripts for **(1)** stimulus presentation via Psychtoolbox-3 to measure pSFT with fMRI (see `/measure-pSFT`) and **(2)** voxel-wise parameter optimization (see `/estimate-pSFT`). 

**Requirements**
- Psychtoolbox-3 must be installed for stimulus presentation.
- The Optimization Toolbox for MATLAB must be installed for parameter optimization.
- The shape of the BOLD percent signal change time series data must have time along the first dimension and voxels along the second (i.e., time x voxels).
- **Strongly recommended**: The Parallel Computing Toolbox for MATLAB must be installed for parallelization.

The repository can be cloned as is.

**Citation**: If you use this toolbox in your research, please cite using the information provided in `CITATION.cff`. 
 
## measure-pSFT
This directory contains scripts for executing the experiment via Psychtoolbox.

We provide an example scan session script for data acquisition (see `/measure-pSFT/run_session.m`) that can be modified with respect to the experimental setup. For example, the input device name, toggles (e.g., save run info), subject ID, directories, and screen parameters should be verified by the user.  

Critical functions include `prepareScan` and `presentStimuli`. 

Users will find key stimulus and timing parameters inside `prepareScan`. For example, to adjust the size of the stimulus, the user must change `p.aperture_radius_deg` (or `p.aperture_radius_px`); to match the fMRI scan length, `t.TR` must match the duration of the repetition time. 

`presentStimuli` will output a structure, `run_info` that compiles all the experiment's structures (scan parameters `p`, timing parameters `t`, window paramters `w`, frame sequences `frames`, and behavioral data `behav_data`). While already in `run_info`, the matrix containing the SF input time series for every block is stored as a separate `.mat` file for convenience, as the time series across multiple blocks and runs should be concatenated as an input vector (i.e., time x 1) into the pSFT optimization pipeline (see `generateSFTimeSeries`).


**Directory contents**
-   `run_session.m`: Example scan session script for data acquisition
-   `/data`: Experimental run info will be stored here by default.
-   `/stimuli`: Stimulus textures will be stored here by default.
    - `verify_stimuli.m`: Analyzes spectral energy of stimuli. Can generate and save experimental stimuli as well. 
-   `/functions`: Contains supporting functions for stimulus generation, display, and experimental control.
    -   `checkPTB`: Verifies Psychtoolbox installation.
    -   `prepareScan`: Initializes parameters, stimuli, timing, and Psychtoolbox window.
    -   `createTextures`: Creates bandpass-filtered noise textures for stimuli.
    -   `createApertures`: Creates stimulus apertures.
    -   `genFrames`: Generates the sequence of events and timing for each experiment frame.
    -   `presentStimuli`: Draws stimuli frame by frame. Compiles run information (e.g., parameters, behavioral data) into struct `run_info`.
    -   `generateSFTimeSeries`: Generates a vectorized SF input time series.
    -   `stimulusParams`: Defines stimulus parameters (aperture radius, contrast, noise filter count, noise sample count).


## estimate-pSFT
This directory contains scripts for estimating pSFT parameters from fMRI data.

We include an example workflow for estimating pSFT from a sample dataset that contains concatenated SF input and measured BOLD time series across 9 scan runs from two subjects — 100 voxels in V1, V2, and V3 (see `/estimate-pSFT/example_pipeline.m`). `sample_data` is a structure array with fields `I` and `measured_BOLD`. Note that this example pipeline assumes that each subject has the same number of regions of interest (ROIs).

`estimatePSFT` is the main high-level function for estimating pSFT parameters. 

It takes the stimulus spatial frequency time series, measured BOLD time series, and a hemodynamic response function (HRF) as input to return a structure `pSFT` containing:
- estimated pSFT parameters (peak SF, bandwidth, BOLD amplitude, baseline)
- estimated pSFT curves
- estimated neural time series
- estimated BOLD time series
- $R^2$ values
- SSE values
- `fmincon` exit flags

Below are toggles and parameters that must be defined before entering `estimatePSFT` (see `example_pipeline`).

Toggles:
- Parallelization (true/false)
- Coarse grid search (true/false)
- Fine grid search (true/false)

Parameters:
- Spatial frequencies used to generate tuning curves (`p.sfs`)
- Initial pSFT parameters (`p.init_params`)
- pSFT parameter bounds (`p.pSFT_bounds`)

**Directory contents**
-   `example_pipeline.m`: Demonstrates a complete workflow for estimating pSFT parameters using sample data. Includes setting up estimation settings (parallelization, grid search, parameter bounds, HRF definition), visualizing results, and saving results. Results are saved under `/estimates/` with the file format `all_pSFT_n_yyyy-MM-dd_HH-mm-ss.mat` that contains the subject x ROI struct array `all_pSFT`.
-   `/functions`: Contains supporting functions:
    -   `estimatePSFT`: Main high-level function for estimating pSFT parameters.
    -   `fitVoxels`: Performs voxel-wise parameter estimation using `fmincon`, called within estimatePSFT.
    -   `logGauss`: Defines the log Gaussian function used for the pSFT model.
    -   `evaluatePSFT`: Evaluates goodness of fit (SSE) between measured and estimated BOLD.
    -   `generateBOLD`: Generates BOLD response from neural response by convolving with HRF and applying amplitude and baseline scaling.
    -   `defineHRF`: Creates a canonical HRF model based on Boynton & Heeger 1996 Journal of Neuroscience.
    -   `gridSearch`: Implements grid search for initial parameter estimates.
    -   `chunkTimeSeries`: Splits voxel time series into chunks for parallel processing.
    -   `cpd2oct`: Converts pSFT bandwidth from cycles per degree of visual angle to octaves.
    -   `checkRequiredToolboxes`: Verifies that the required MATLAB Toolboxes are installed.
    -   `R2`: Calculates coefficient of determination (R²) between measured and estimated time series.
    -   `SSE`: Calculates sum of squared errors between measured and estimated time series.
    -   `plotSettings`: Returns a structure of standardized plot settings (fonts, colors, sizes) for visualization.
-   `/validate-pSFT`: Contains validation scripts for testing the pSFT estimation pipeline.
    -   `simulate_pSFT.m`: Simulates pSFT curves and demonstrates how changes in peak SF (μ) and bandwidth (σ) affect the tuning function.
    -   `validate_pSFT.m`: Validates the pSFT estimation pipeline by simulating stimulus and BOLD time series from pre-defined pSFT parameters, then estimating parameters and comparing estimates to ground truth. Results are saved under `/validate-pSFT/estimates/` and `/validate-pSFT/figures/`.
    -   `/functions`: Contains modular functions for validation:
        -   `generateSimulatedBOLD.m`: Generates clean BOLD responses from ground truth pSFT parameters and adds Gaussian noise at multiple SNR levels.
        -   `computeValidationMetrics.m`: Computes validation metrics (RMSE, correlation, R²) for each SNR level, both per-subject/ROI and pooled across all voxels, with bootstrap confidence intervals.

### validate-pSFT

This subdirectory contains scripts for validating and simulating pSFT. 

This sub-module demonstrates the effectiveness of the model fitting pipeline. Inspired by the validation framework of Lerma-Usabiaga et al. (2020), `validate_pSFT.m` generates synthetic BOLD data from known parameters (log-Gaussian model convolved with a canonical HRF) across three SNR levels: 5.29, -0.51, and -4.29 dB. The script runs the estimation pipeline on this data and benchmarks performance by comparing recovered parameters to ground truth. Metrics include RMSE and Pearson correlations with bootstrap confidence intervals. The module automatically saves timestamped results and visualization reports.

`simulate_pSFT.m` demonstrates how the log Gaussian pSFT model responds to changes in peak spatial frequency and bandwidth parameters. 

## License

The Population Spatial Frequency Toolbox is released under the GNU General Public License v3.0 or later. See the [LICENSE.txt](LICENSE.txt) file for details.



