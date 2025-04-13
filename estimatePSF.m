% estimatePSF 
%   Estimates pSF parameters from voxel time series via fmincon()
%
%   Inputs:
%       measured_BOLD - voxel time series [voxels x time]
%       I - input time series [1 x time]
%       HIRF - hemodynamic impulse response function [1 x time]
%       p - structure of parameters (see example_pipeline.m)
%       toggles - structure of toggles (see example_pipeline.m)
%
%   Outputs:
%       pSF - struct containing: 
%       -   estimated pSF parameters [voxels x 4]
%       -   estimated pSF curves [voxels x sf_count]
%       -   estimated neural time series [voxels x time]
%       -   estimated BOLD time series [voxels x time]
%       -   estimated R^2 values [voxels x 1]
%       -   estimated SSE values [voxels x 1]
%       -   fmincon exit flags [voxels x 1]

function pSF = estimatePSF(measured_BOLD, I, HIRF, p, toggles)

    if toggles.disp_on, disp('Initializing estimatePSF ...'); end

    addpath('functions');
      
    %% Initialize 'chunks' structure array 
        
    chunks_size = cell(1, p.num_chunks);    

    chunks = struct('vox_indices', chunks_size, ...
        'param_est', chunks_size, ...    
        'est_SFT', chunks_size, ...
        'est_R', chunks_size, ...
        'est_BOLD', chunks_size, ...
        'r2', chunks_size, ...
        'sse', chunks_size, ...
        'start_values', chunks_size, ...
        'start_sse', chunks_size, ...
        'exitflag', chunks_size, ...
        'measured_BOLD', chunks_size);

    %% Chunk time series data

    if toggles.disp_on, disp('Chunking time series data...'); end

    [measured_BOLD_chunks, vox_chunk_indices] = chunkTimeSeries(p.num_chunks, measured_BOLD);

    for chunk = 1:p.num_chunks
        chunks(chunk).measured_BOLD = measured_BOLD_chunks{chunk};
        chunks(chunk).vox_indices = vox_chunk_indices(chunk,:);
    end
    
    %% Loop through chunks of voxels
    
    if toggles.disp_on, disp('Fitting voxels...'); end
    
    if toggles.parallelization
        parfor chunk = 1:p.num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk), I, HIRF, p, toggles);
        end
    else
        for chunk = 1:p.num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk), I, HIRF, p, toggles);
        end
    end
    
    %% Resplice chunks

    if toggles.disp_on, disp('Resplicing chunks...'); end

    field_names = fieldnames(chunks);

    for i_field = 1:length(field_names)
        pSF.(field_names{i_field}) = cat(1, chunks.(field_names{i_field}));
    end

end