% prepareScan
%   Prepares the scan for the experiment.
%
% Syntax
%   [p, w, t, stimuli, frames] = prepareScan(p, w, t, dirs, toggles)
%
% Input Arguments
%   p – structure containing scan parameters
%   w – structure containing window parameters
%   t – structure containing timing parameters
%   dirs – structure containing directories
%   toggles – structure containing toggles
%
% Output Arguments
%   p – structure containing scan parameters
%   w – structure containing window parameters
%   t – structure containing timing parameters
%   stimuli – structure containing textures and stimuli for drawing
%   frames – structure containing frame deadlines and event control

function [p, w, t, stimuli, frames] = prepareScan(p, w, t, dirs, toggles)
    
    %% Scan parameters

    p.num_blocks = 6;
    
    % +++ Stimuli parameters +++
    p.aperture_radius_deg = 9; % defines the apparent stimulus size
    p.aperture_radius_px = round(p.aperture_radius_deg * w.ppd);
    p.stimulus_radius_px = round(p.aperture_radius_px * 1.25); % default = 1.1

    p.stimulus_contrast = 0.9; % default = 0.9
    p.noise_filter_count = 40; % default = 40
    p.noise_sample_count = 10; % default = 10

    % ------------------------------

    % +++ Fixation parameters +++
    p.fixation_task_luminance = 70; % default = 70 

    p.fixation_aperture_deg = 1;
    p.fixation_aperture_px = round(p.fixation_aperture_deg * w.ppd);
    if ~mod(p.fixation_aperture_px, 2), p.fixation_aperture_px = p.fixation_aperture_px + 1; end 

    % The fixation aperture radius defines the inner annulus radius.
    p.fixation_aperture_radius_deg = 0.32; % default = 0.32°
    p.fixation_aperture_radius_px = round(p.fixation_aperture_radius_deg * w.ppd); 
    if ~mod(p.fixation_aperture_radius_px, 2), p.fixation_aperture_radius_px = p.fixation_aperture_radius_px + 1; end 

    p.fixation_dot_deg = 0.15; % default = 0.15°
    p.fixation_dot_px = round(p.fixation_dot_deg * w.ppd); 
    if ~mod(p.fixation_dot_px, 2), p.fixation_dot_px = p.fixation_dot_px - 1; end 
    
    % ------------------------------

    % +++ Timing parameters +++
    t.TR = 1; % fMRI TR duration in seconds
    t.blank_period_dur = t.TR * 10; % default = 10 s
    t.block_dur = t.TR * p.noise_filter_count;
    t.scan_dur = (t.block_dur * p.num_blocks) + (t.blank_period_dur * (p.num_blocks + 1));

    t.noise_sample_update_freq = 10; % Hz, default = 10
    t.noise_SF_update_freq = 1; % Hz, default = 1

    t.fixation_target_dur = 0.250; % default = 0.250
    t.min_time_btwn_fixation_targets = 0.600; % default = 0.600
    t.fixation_target_delay = 1; % default = 1
    t.fixation_task_resp_window = 1; % default = 1

    t.frame_rate = 20; % Hz, default = 20
    t.frame_dur = 1/t.frame_rate;
   
    % ------------------------------

    %% Device input

    p.trigger_key = KbName('5%');
    p.keypress_numbers = [KbName('1!') KbName('2@') KbName('3#') KbName('4$')];

    p.device_number = 0;
    [Kb_indices, product_names, ~] = GetKeyboardIndices;

    for i = 1:length(product_names)
        if strcmp(product_names{i}, p.device_string)
            p.device_number = Kb_indices(i);
            break;
        end
    end

    if p.device_number == 0
        error(['Input device ' p.device_string ' not detected.']);
    else
        disp(['Input device ' p.device_string ' detected.']);
    end

    %% Load/create textures

    % pSF textures
    texture_filepath = fullfile(dirs.stimuli_dir, ...
        sprintf('bandpass_filtered_noise_%s_%s_%s.mat', ...
                strrep(sprintf('%0.2f', w.ppd), '.', ''), ...
                strrep(sprintf('%0.6f', w.px_size), '.', ''), ...
                strrep(sprintf('%0.2f', p.stimulus_radius_px), '.', '')));

    if exist(texture_filepath, 'file') == 2
        disp('Loading textures...');
        load(texture_filepath, 'textures', 'filters');
    else
        disp('Creating textures...');
        [textures, filters] = createTextures(p.stimulus_radius_px, p.stimulus_contrast, p.noise_filter_count, p.noise_sample_count, w.ppd, w.px_size, toggles.save_textures, texture_filepath);
    end

    stimuli.pSF_textures = textures;
    stimuli.filters = filters;

    disp('Textures in memory.');
    clear textures stimulus filters;

    % Apertures
    [stimuli.stimulus_aperture, stimuli.fixation_aperture] = createApertures(p, w);

    %% Generate frames

    frames = genFrames(p, t);

    %% Open window 

    if toggles.macOS, Screen('Preference', 'SkipSyncTests', 1); end

    w.window = PsychImaging('OpenWindow', w.use_screen, w.gray, [0 0 w.screen_width_px w.screen_height_px]);
    HideCursor;

    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

    Screen('BlendFunction', w.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

    w.default_CLUT = Screen('ReadNormalizedGammaTable', w.window);
    if toggles.gamma_correction
        Screen('LoadCLUT', w.window, w.corrected_CLUT);
    end

    %% Make stimuli for drawing

    stimuli.stimulus_aperture_made = Screen('MakeTexture', w.window, stimuli.stimulus_aperture);
    stimuli.fixation_aperture_made = Screen('MakeTexture', w.window, stimuli.fixation_aperture);

    stimuli.pSF_textures_made = nan(p.noise_filter_count, p.noise_sample_count);
    for n_filter = 1:p.noise_filter_count
        for n_sample = 1:p.noise_sample_count
            stimuli.pSF_textures_made(n_filter, n_sample) = Screen('MakeTexture', w.window, stimuli.pSF_textures(:, :, n_filter, n_sample));
        end
    end

    %%  Patches
    
    stimuli.pSF_stimuli_patch = CenterRectOnPoint([0 0 w.screen_height_px w.screen_height_px], w.centerX, w.centerY);
    stimuli.stimulus_aperture_patch = CenterRectOnPoint([0 0 w.screen_width_px w.screen_height_px], w.centerX, w.centerY);
    stimuli.fixation_aperture_patch = CenterRectOnPoint([0 0 p.fixation_aperture_px p.fixation_aperture_px], w.centerX, w.centerY);
    stimuli.fixation_dot_patch = CenterRectOnPoint([0 0 p.fixation_dot_px p.fixation_dot_px], w.centerX, w.centerY);

end