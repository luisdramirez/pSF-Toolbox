% estimatePSF 
%   Estimates pSF parameters from voxel time series via fmincon()
%
%   Inputs:
%       measured_BOLD - voxel time series [voxels x time]
%       I - input time series [1 x time]
%
%   Outputs:
%       pSFT - struct containing: 
%       -   estimated pSFT parameters [voxels x 4]
%       -   estimated SFT curves [voxels x sf_count]
%       -   estimated neural time series [voxels x time]
%       -   estimated BOLD time series [voxels x time]
%       -   estimated R^2 values [voxels x 1]
%       -   estimated SSE values [voxels x 1]
%       -   fmincon exit flags [voxels x 1]

function [pSFT, p, toggles] = estimatePSF(measured_BOLD, I, HIRF)

    tic;
    
    addpath('/functions');
    addpath('/script_modules');
    
    %% Toggles
    
    toggles.parallelization = 1;
    
    toggles.coarse_grid_search = 1;
    toggles.fine_grid_search = 1;
    
    toggles.make_voxel_pSFT_plots = 0;  % Plot individual voxel fits
    
    %% Parallelization setup for parfor loop
    
    num_cores = 8;
    num_chunks = num_cores-1;
    
    if toggles.parallelization
        maxNumCompThreads(num_cores);
        if num_cores > 1
            parpool('local', num_chunks)
        end
        toggles.make_voxel_pSFT_plots = 0; % note: figures can't be drawn during parfor loops, so this will force plot generation off.
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
  
    %% Initialize 'chunks' structure array 
        
    chunk_size = cell(1, num_vox_chunks);    

    chunks = struct('vox_indices', chunk_size, ...
        'param_est', chunk_size, ...    
        'est_SFT', chunk_size, ...
        'est_R', chunk_size, ...
        'est_BOLD', chunk_size, ...
        'r2', chunk_size, ...
        'sse', chunk_size, ...
        'start_values', chunk_size, ...
        'start_sse', chunk_size, ...
        'exitflag', chunk_size, ...
        'measured_BOLD', chunk_size);

    %% Chunk time series data

    [measured_BOLD_chunks, vox_chunk_indices] = chunkTimeSeries(num_chunks, measured_BOLD);

    for chunk = 1:num_chunks
        chunks(chunk).measured_BOLD = measured_BOLD_chunks{chunk};
        chunks(chunk).vox_indices = vox_chunk_indices(chunk,:);
    end
    
    %% Loop through chunks of voxels
    
    if toggles.parallelization
        parfor chunk = 1:num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk).measured_BOLD, I, HIRF, p, toggles);
        end
        delete(gcp);
    else
        for chunk = 1:num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk).measured_BOLD, I, HIRF, p, toggles);
        end
    end
    
    %% Resplice chunks

    field_names = fieldnames(chunks);

    for i_field = 1:length(field_names)
        pSFT.(field_names{i_field}) = cat(1, chunks.(field_names{i_field}));
    end
    
    %% Print elapsed time

    disp(['Elapsed time: ~' num2str(round(toc/60,1)) ' minute(s).']);

end