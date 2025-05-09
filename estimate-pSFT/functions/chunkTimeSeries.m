% chunkTimeSeries - Create chunks of time series data 
%   Assumes the time series data shape is a 2D matrix with time along the first dimension (ie time x voxels).
%
%   Syntax
%       [time_series_chunks, chunk_indices] = chunkTimeSeries(time_series_data, num_chunks)
%
%   Input Arguments
%       time_series_data – time series data [time x voxels]
%       num_chunks – number of chunks to create
%
%   Output Arguments
%       time_series_chunks – cell array of time series data chunks
%       chunk_indices – matrix of chunk indices

function [time_series_chunks, chunk_indices] = chunkTimeSeries(time_series_data, num_chunks)

    %% Create chunk indices

    chunk_indices = nan(num_chunks,2); % [starting row index, ending row index]
    
    num_voxels = size(time_series_data,2);
    chunk_size = floor(num_voxels/num_chunks);
    chunk_starts = 1:chunk_size:num_voxels;

    for i_chunk = 1:num_chunks
        if i_chunk == num_chunks
            chunk_indices(i_chunk,:) = [chunk_starts(i_chunk), num_voxels];
        else
            chunk_indices(i_chunk,:) = [chunk_starts(i_chunk), chunk_starts(i_chunk+1)-1];
        end
    end
    
    %% Create chunks of time series data
    
    time_series_chunks = cell(1,num_chunks);
    
    for i_chunk = 1:num_chunks
        time_series_chunks{i_chunk} = time_series_data(:,chunk_indices(i_chunk,1):chunk_indices(i_chunk,2));            
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
% Contact luisdramirez95@gmail.com for any questions or comments.