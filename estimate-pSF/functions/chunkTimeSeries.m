% chunkTimeSeries - Create chunks of time series data 
%   Assumes the time series data shape is a 2D matrix with time along the columns.
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