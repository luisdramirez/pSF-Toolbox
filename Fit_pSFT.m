% Fit_pSFT 
%   Estimates pSFT parameters from voxel time series via fmincon()
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
%       -   exit flags [voxels x 1]

function [pSFT, script_t] = Fit_pSFT(measured_BOLD, I, HIRF)

    tic;

    curr_time = datetime;
    curr_time.Format = 'MM.dd.yyyy_HHmm';
    curr_time = char(curr_time);
    
    % Use time as rng seed
    script_t.rng_seed = sum(100*clock);
    rng(script_t.rng_seed);
    
    % Add paths
    addpath('/functions');
    addpath('/script_modules');
    
    %% Toggles
    
    % Parallelization
    toggles.parallelization = 1;
    
    % Grid searches
    toggles.coarse_grid_search = 1;
    toggles.fine_grid_search = 1;
    
    % Plot individual voxel fits
    toggles.make_voxel_pSFT_plots = 0;
    
    %% Parallelization setup for parfor loop
    
    num_cores = 16;
    num_chunks = num_cores-1;
    
    if toggles.parallelization
        maxNumCompThreads(num_cores);
        if num_cores > 1
            parpool('local', num_chunks)
        end
        toggles.make_voxel_pSFT_plots = 0; % note: figures can't be drawn during parfor loops, so this will force plot generation off in case it's on.
    end
    
    %% Spatial frequency parameters
    % These are hardcoded SF parameters used for generating pSFT curves.
    % Ideally, the range should match the range of SFs presented.
    
    p.sf_min = 0.1;
    p.sf_max = 12;
    
    p.sf_count = 100;
    p.sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), p.sf_count); % Values used to generate smooth tuning curves
    
    %% Initalize pSFT model parameters and bounds
    % [mu, sigma, beta, beta_0]
    
    p.init_params = [1 1 1 0]; 
    
    p.pSFT_bounds(1,:) = [6, 4, 25, 10];
    p.pSFT_bounds(2,:) = [0.009, 0.1, -25, -10]; 
  
    %% Initialize chunks structure array 
        
    chunk_size = cell(1, num_vox_chunks);    

    chunks = struct('vox_inds', chunk_size, ...
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
        
    total_num_vox = size(measured_BOLD,1);
    
    [measured_BOLD_chunks, vox_chunk_indices] = chunkTimeSeries(total_num_vox, num_chunks, measured_BOLD);

    for chunk = 1:num_chunks
        chunks(chunk).measured_BOLD = measured_BOLD_chunks{chunk};
        chunks(chunk).vox_inds = vox_chunk_indices(chunk,:);
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

    elapsed_time = toc;
    script_t.elapsed_time = round(elapsed_time/60);
    disp(' '); disp(['Elapsed time: ~' num2str(script_t.elapsed_time) ' minute(s).']); disp(' ');


end