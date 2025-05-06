% estimatePSF 
%   Estimates pSFT parameters from voxel time series via fmincon
%
%   Inputs:
%       I - input time series [time x 1]
%       measured_BOLD - voxel time series [time x voxels]
%       HIRF - hemodynamic impulse response function [time x 1]
%       p - structure of parameters (see example_pipeline.m)
%       toggles - structure of toggles (see example_pipeline.m)
%
%   Outputs:
%       pSFT - struct containing: 
%       -   estimated pSFT parameters [4 x voxels]
%       -   estimated pSFT curves [sf_count x voxels]
%       -   estimated neural time series [time x voxels]
%       -   estimated BOLD time series [time x voxels]
%       -   estimated R^2 values [1 x voxels]
%       -   estimated SSE values [1 x voxels]
%       -   fmincon exit flags [1 x voxels]
%       -   measured BOLD time series [time x voxels]

function pSFT = estimatePSF(I, measured_BOLD, HIRF, p, toggles)

    if toggles.disp_on, disp('Initializing estimatePSF ...'); end

    addpath('functions');
      
    %% Initialize 'chunks' structure array 
        
    chunks_size = cell(1, p.num_chunks);    

    chunks = struct('vox_indices', chunks_size, ...
        'measured_BOLD', chunks_size, ...
        'param_est', chunks_size, ...    
        'est_SFT', chunks_size, ...
        'est_R', chunks_size, ...
        'est_BOLD', chunks_size, ...
        'r2', chunks_size, ...
        'sse', chunks_size, ...
        'start_values', chunks_size, ...
        'start_sse', chunks_size, ...
        'exitflag', chunks_size);

    %% Chunk time series data

    if toggles.disp_on, disp('Chunking time series data...'); end

    [measured_BOLD_chunks, vox_chunk_indices] = chunkTimeSeries(measured_BOLD, p.num_chunks);

    for chunk = 1:p.num_chunks
        chunks(chunk).vox_indices = vox_chunk_indices(chunk,:);
        chunks(chunk).measured_BOLD = measured_BOLD_chunks{chunk};
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
        pSFT.(field_names{i_field}) = cat(2, chunks.(field_names{i_field}));
    end

end