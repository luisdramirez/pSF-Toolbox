function frames = genFrames(p, t)

    frames.onsets = 0:t.frame_dur:t.scan_dur-t.frame_dur;
    frames.count = length(frames.onsets);

    %% Initialize frames

    frames.block = zeros(1, frames.count);
    frames.update_block = zeros(1, frames.count);
    frames.update_noise_sf = zeros(1, frames.count);
    frames.update_noise_sample = zeros(1, frames.count);
    frames.target_on = zeros(1, frames.count);
    frames.response_window = zeros(1, frames.count);
    frames.target_onset_indices = zeros(1, frames.count);
    frames.target_response_map = zeros(1, frames.count);

    %% Initialize parameters

    blank_period_frames_count = round(t.blank_period_dur * t.frame_rate);
    block_frames_count = round(t.block_dur * t.frame_rate);
    noise_SF_update_frames = round(t.frame_rate / t.noise_SF_update_freq); % Frames between SF updates (1 Hz -> every 20 frames)
    noise_sample_update_frames = round(t.frame_rate / t.noise_sample_update_freq); % Frames between sample updates (10 Hz -> every 2 frames)
    target_frames_count = round(t.fixation_target_dur * t.frame_rate);
    min_gap_frames_count = round(t.min_time_btwn_fixation_targets * t.frame_rate);
    delay_frames_count = round(t.fixation_target_delay * t.frame_rate);
    resp_window_frames_count = round(t.fixation_task_resp_window * t.frame_rate);
    min_onset_distance = target_frames_count + min_gap_frames_count;

    %% Generate frames

    target_counter = 0;

    for n_block = 1:p.num_blocks

        block_start_frame = n_block * blank_period_frames_count + (n_block - 1) * block_frames_count + 1;
        block_end_frame = block_start_frame + block_frames_count - 1; % n_block * blank_period_frames_count + n_block * block_frames_count;
        frame_range = block_start_frame:block_end_frame;
        
        frames.block(frame_range) = 1;
        if n_block < p.num_blocks
            % update block number at the end of the rest period, just before the next block starts
            frames.update_block(frame_range(end) + blank_period_frames_count) = 1; 
        end

        frames.update_noise_sf(frame_range(1:noise_SF_update_frames:end)) = 1;
        frames.update_noise_sample(frame_range(1:noise_sample_update_frames:end)) = 1;

        %% Fixation Task
        
        % Determine valid range for target ONSETS within this block
        earliest_onset_frame = block_start_frame + delay_frames_count;
        latest_onset_frame = block_end_frame - target_frames_count + 1;
        possible_onset_frames = earliest_onset_frame:latest_onset_frame;
        
        selected_target_onsets = [];
        
        while ~isempty(possible_onset_frames)
            % Select a random potential onset frame
            rand_idx = randi(length(possible_onset_frames));
            selected_onset = possible_onset_frames(rand_idx);
            selected_target_onsets = [selected_target_onsets, selected_onset]; %#ok<AGROW>
            
            % Determine the 'keep out' zone for *other* potential onsets
            % Remove frames too close (before or after) the selected onset
            keep_out_start = selected_onset - min_onset_distance + 1;
            keep_out_end = selected_onset + min_onset_distance - 1;
            
            % Remove the frames within the keep out zone from the possible onsets
            possible_onset_frames(possible_onset_frames >= keep_out_start & possible_onset_frames <= keep_out_end) = [];
        end
        
        % Assign target_on and response_window based on selected onsets
        selected_target_onsets = sort(selected_target_onsets);
        for i = 1:length(selected_target_onsets)
            onset_frame = selected_target_onsets(i);

            if datasample([0 1], 1)
                target_counter = target_counter + 1;
                
                frames.target_onset_indices(onset_frame) = target_counter;

                % Target period
                target_end_frame = onset_frame + target_frames_count - 1;
                % Ensure target does not exceed total scan duration (shouldn't happen with current logic, but good practice)
                target_end_frame = min(target_end_frame, frames.count); 
                frames.target_on(onset_frame:target_end_frame) = 1;
                
                % Response window period
                resp_window_end_frame = onset_frame + resp_window_frames_count - 1;
                % Ensure response window does not exceed total scan duration
                resp_window_end_frame = min(resp_window_end_frame, frames.count);
                frames.response_window(onset_frame:resp_window_end_frame) = 1;
                
                % Map these response window frames to the current target's unique index
                % If windows overlap, the later target's index will overwrite the earlier one's
                frames.target_response_map(onset_frame:resp_window_end_frame) = target_counter; 
            end
        end

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