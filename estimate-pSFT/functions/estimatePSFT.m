% estimatePSFT 
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

function pSFT = estimatePSFT(I, measured_BOLD, HIRF, p, toggles)

    if toggles.disp_on, disp('Initializing estimatePSFT ...'); end

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

%% 

           % The Population Spatial Frequency Toolbox
% Copyright (C) 2025 Luis D. Ramirez, Feiyi Wang, Emily Wiecek, Louis N. Vinke, and Sam Ling
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% Contact luisdramirez95@gmail.com for any questions or comments.