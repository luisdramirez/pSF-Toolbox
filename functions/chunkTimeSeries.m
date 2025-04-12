% chunkTimeSeries - Create chunks of time series data 
%   Assumes the time series data shape is a 2D matrix with time along the columns.
%
%   Syntax
%       [time_series_chunks, chunk_indices] = chunkTimeSeries(num_rows, num_chunks, time_series_data)
%
%   Input Arguments
%       num_rows – number of rows in the time series data
%       num_chunks – number of chunks to create
%       time_series_data – time series data [voxels x time]
%
%   Output Arguments
%       time_series_chunks – cell array of time series data chunks
%       chunk_indices – matrix of chunk indices

function [time_series_chunks, chunk_indices] = chunkTimeSeries(num_chunks, time_series_data)

    %% Create chunk indices

    chunk_indices = nan(num_chunks,2); % [starting row index, ending row index]
    
    num_rows = size(time_series_data,1);
    chunk_size = floor(num_rows/num_chunks);
    chunk_starts = 1:chunk_size:num_rows;

    for i_chunk = 1:num_chunks
        if i_chunk == num_chunks
            chunk_indices(i_chunk,:) = [chunk_starts(i_chunk), num_rows];
        else
            chunk_indices(i_chunk,:) = [chunk_starts(i_chunk), chunk_starts(i_chunk+1)-1];
        end
    end
    
    %% Create chunks of time series data
    time_series_chunks = cell(1,num_chunks);
    
    for i_chunk = 1:num_chunks
        time_series_chunks{i_chunk} = time_series_data(chunk_indices(i_chunk,1):chunk_indices(i_chunk,2),:);            
    end

end