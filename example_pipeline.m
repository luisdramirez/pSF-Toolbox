% example_pipeline


%% Prepare workspace
clear all; close all; clc;

addpath('functions');
addpath('data');
addpath('estimates')

%% Load data

load('data/sample_data.mat');


%% Estimate PSF

num_subjs = 

HIRF = defineHIRF();
pSFT = estimatePSF(example_data.measured_BOLD, example_data.I, HIRF);

