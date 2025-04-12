# Population spatial frequency Toolbox

The pSF Toolbox was designed to streamline the population spatial frequency tuning (pSFT) approach originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019) with accessibility in mind. 

We provide code for estimating pSFT, simulating pSFT, and generating the stimulus presentation used for pSFT mapping.

We also provide a "toy" dataset from one of our studies that used the pSFT approach.

## Fit_pSFT.m
This is the main function that estimates pSFT from a measured BOLD time series, input spatial frequeny time series, and hemodynamic impluse response function. As long as the inputs are in the required format (e.g., the shape of the measured BOLD time series matrix has voxels along the rows and time points along the columns), you will receive a structure, `pSFT`, that contains:
- estimated pSFT parameters 
- estimated SFT curves
- estimated neural time series
- estimated BOLD time series
- estimated $R^2$ values
- estimated SSE values
- fmincon exit flags

### fitVoxels.m
Here you'll find where the magic happens


