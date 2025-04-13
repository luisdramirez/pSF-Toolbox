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
%       pSFT - struct containing: 
%       -   estimated pSFT parameters [voxels x 4]
%       -   estimated SFT curves [voxels x sf_count]
%       -   estimated neural time series [voxels x time]
%       -   estimated BOLD time series [voxels x time]
%       -   estimated R^2 values [voxels x 1]
%       -   estimated SSE values [voxels x 1]
%       -   fmincon exit flags [voxels x 1]

function [pSFT, p, toggles] = estimatePSF(measured_BOLD, I, HIRF, p, toggles)

    tic;

    if toggles.disp_on
        disp('Initializing estimatePSF ...');
    end
      
    %% Initialize 'chunks' structure array 
        
    chunk_size = cell(1, p.num_chunks);    

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

    if toggles.disp_on
        disp('Chunking time series data...');
    end

    [measured_BOLD_chunks, vox_chunk_indices] = chunkTimeSeries(p.num_chunks, measured_BOLD);

    for chunk = 1:p.num_chunks
        chunks(chunk).measured_BOLD = measured_BOLD_chunks{chunk};
        chunks(chunk).vox_indices = vox_chunk_indices(chunk,:);
    end
    
    %% Loop through chunks of voxels
    
    if toggles.disp_on
        disp('Fitting voxels...');
    end
    
    if toggles.parallelization
        parfor chunk = 1:p.num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk), I, HIRF, p, toggles);
        end
        delete(gcp);
    else
        for chunk = 1:p.num_chunks
            chunks(chunk) = fitVoxels(chunks(chunk), I, HIRF, p, toggles);
        end
    end
    
    %% Resplice chunks

    if toggles.disp_on
        disp('Resplicing chunks...');
    end

    field_names = fieldnames(chunks);

    for i_field = 1:length(field_names)
        pSFT.(field_names{i_field}) = cat(1, chunks.(field_names{i_field}));
    end
    
    %% Print elapsed time

    p.elapsed_time = round(toc/60,1);
    if toggles.disp_on
        disp(['Elapsed time: ~' num2str(p.elapsed_time) ' minute(s).']);
        disp(' ');
    end

end