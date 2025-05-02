% presentStimuli
%   Presents SF bandpass filtered noise stimuli to the subject for a fixed number of blocks.
%   The subject responds to a change in luminance at fixation with a button press.
%   The subject is provided a blank period between blocks.
%
% Syntax
%   run_info = presentStimuli(p, w, t, stimuli)
%
% Input Arguments
%   p – structure containing scan parameters
%   w – structure containing window parameters
%   t – structure containing timing parameters
%   stimuli – structure containing stimuli
%   frames – structure containing frame parameters

% Output Arguments
%   run_info – structure containing p, w, t, frames, behav_data

function run_info = presentStimuli(p, w, t, stimuli, frames)
        
    %% Generate presentation sequences

    p.I = nan(p.noise_filter_count, p.num_blocks);
    for n_block = 1:p.num_blocks
        p.I(:,n_block) = datasample(1:p.noise_filter_count, p.noise_filter_count, 'Replace', false); 
    end

    %% Wait for trigger 

    DrawFormattedText(w.window, '~', w.ppd, w.ppd, w.white, w.gray);
    Screen('Flip', w.window);

    disp('Waiting for trigger...')
    KbTriggerWait(p.trigger_key, p.device_number);
    t.trigger_time_stamp = GetSecs;
    disp('Trigger detected!')

    PsychHID('KbQueueCreate', p.device_number);
    PsychHID('KbQueueStart', p.device_number);

    %% Initialize behavioral data structure

    num_total_targets = max(frames.target_onset_indices);
    behav_data.detection = false(1, num_total_targets);       
    behav_data.response_time = nan(1, num_total_targets);  % RT relative to target onset

    %% Stimulus presentation

    n_block = 1;
    sf_indx = 0;
    noise_sample_indx = 0;

    for n_frame = 1:frames.count

        if n_frame == 1
            frames.flip_times = frames.onsets + GetSecs;
            t.run_start = frames.flip_times(1);
        end

        %% Update noise SF

        if frames.update_noise_sf(n_frame) && sf_indx < p.noise_filter_count
            sf_indx = sf_indx + 1;
            noise_sample_indx = 0;
        end

        if frames.update_noise_sample(n_frame) && noise_sample_indx < p.noise_sample_count
            noise_sample_indx = noise_sample_indx + 1;
        end

        %% Update fixation dot

        if frames.target_on(n_frame)
            fixation_color = p.fixation_task_luminance;
        else
            fixation_color = w.black;
        end

        %% Draw stimuli
        
        if frames.block(n_frame)

            Screen('DrawTexture', w.window, stimuli.pSF_textures_made(p.I(sf_indx, n_block), noise_sample_indx), [], stimuli.pSF_stimuli_patch);
            Screen('DrawTexture', w.window, stimuli.stimulus_aperture_made, [], stimuli.stimulus_aperture_patch);
            Screen('DrawTexture', w.window, stimuli.fixation_aperture_made, [], stimuli.fixation_aperture_patch);
            Screen('FillOval', w.window, fixation_color, stimuli.fixation_dot_patch);
       
        else
            
            Screen('FillOval', w.window, w.white, stimuli.fixation_dot_patch);
        
        end

        [frames.VBLTimestamp(n_frame), ...
            frames.StimulusOnsetTime(n_frame), ...
            frames.FlipTimestamp(n_frame), ...
            frames.Missed(n_frame)] = ...
            Screen('Flip', w.window, frames.flip_times(n_frame));

        %% Check response

        [key_pressed, first_press] = PsychHID('KbQueueCheck', p.device_number);
        which_keys = find(first_press);
        
        if key_pressed && ~isempty(intersect(which_keys, p.keypress_numbers))

            % Check if this frame is within *any* response window
            if frames.response_window(n_frame) == 1
                % Get the index of the target this window belongs to (handles overlaps)
                current_target_index = frames.target_response_map(n_frame);

                % Check if the pressed key is a valid response key
                if current_target_index > 0
                    
                    % Find the onset frame for this target to calculate RT
                    % (Could pre-calculate this, but finding it here is okay too)
                    target_onset_frame_for_rt = find(frames.target_onset_indices == current_target_index, 1);
                    
                    if ~isempty(target_onset_frame_for_rt)
                        behav_data.detection(current_target_index) = true;
                        behav_data.response_time(current_target_index) = GetSecs - frames.StimulusOnsetTime(target_onset_frame_for_rt);
                    end
                end
            end
        end

        %% Update block info 

        if frames.update_block(n_frame)
            n_block = n_block + 1;
            sf_indx = 0;
            PsychHID('KbQueueFlush', p.device_number);
        end

    end

    t.run_end = GetSecs;    
    behav_data.overall_hit_rate = sum(behav_data.detection) / num_total_targets;
    disp(['Hit rate: ' num2str(round(100*behav_data.overall_hit_rate)) '%']);

    %% Compile structures 
        
    run_info.p = p;
    run_info.t = t;
    run_info.w = w;
    run_info.frames = frames;
    run_info.behav_data = behav_data;

end