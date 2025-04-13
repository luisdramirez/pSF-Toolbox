% example_pipeline

%% Prepare workspace
clear all; close all; clc;

addpath('functions');
addpath('data');
addpath('estimates')

%% Toggles

toggles.parallelization = false;
toggles.coarse_grid_search = true;
toggles.fine_grid_search = true;
toggles.make_voxel_plots = false; % Plot individual voxel fits
toggles.disp_on = true; % Display progress

%% Parallelization setup for parfor loop

p.num_cores = 8;
p.num_chunks = p.num_cores-1;

if toggles.parallelization
    maxNumCompThreads(p.num_cores);
    if p.num_cores > 1
        parpool('local', p.num_chunks)
    end
    toggles.make_voxel_plots = false; % Note: figures can't be drawn during parfor loops, so this will force plot generation off.
end

%% Spatial frequency parameters
% These are hardcoded SF parameters used for generating pSFT curves.
% Ideally, the range should match the range of SFs presented.

p.sf_min = 0.1;
p.sf_max = 12;

p.sf_count = 100;
p.sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), p.sf_count);

%% Initalize pSFT model parameters
% [mu, sigma, beta, beta_0]

p.init_params = [1 1 1 0]; 

p.pSFT_bounds(1,:) = [6, 4, 25, 10];
p.pSFT_bounds(2,:) = [0.009, 0.1, -25, -10]; 

%% Load data

load('data/sample_data.mat');

%% Estimate PSF

num_subjs = length(sample_data.measured_BOLD);
num_ROIs = size(sample_data.measured_BOLD{1}, 3);
total_elapsed_time = 0;

HIRF = defineHRF();

for subj = 1:num_subjs
    for roi = 1:num_ROIs

        if toggles.disp_on
            disp(['Estimating pSF for subject ' num2str(subj) ', V' num2str(roi) ' ...']);
        end

        [pSFT, p, ~] = estimatePSF(sample_data.measured_BOLD{subj}(:,:,roi), sample_data.I{subj}, HIRF, p, toggles);

        total_elapsed_time = total_elapsed_time + p.elapsed_time;

    end
end

if toggles.disp_on
    disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']);
end

if toggles.parallelization
    delete(gcp);
end