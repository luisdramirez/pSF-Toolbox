# Population spatial frequency Toolbox

The pSF Toolbox was designed to streamline the population spatial frequency tuning approach originally developed by [Aghajari, Vinke, & Ling 2020 (Journal of Neurophysiology)](https://doi.org/10.1152/jn.00291.2019) with accessibility in mind. 

## Fit_pSFT.m
This is the main function that estimates population spatial frequency tuning (pSFT) from a measured BOLD time series, input spatial frequny time series, and hemodynamic impluse response function. As long as the inputs are in the required format (e.g., the shape of the measured BOLD time series matrix has voxels along the rows and time points along the columns), you will receive a structure, `pSFT`, that contains:
- estimated pSFT parameters 
- estimated SFT curves
- estimated neural time series
- estimated BOLD time series
- estimated R^2 values
- estimated SSE values
- fmincon exit flags

### fitVoxels.m
Here you'll find where the magic happens


