
%% Prepare workspace

clear all; 
close all; 
clc;

% Grab date
t.the_date = string(datetime('now', 'Format', 'yyyyMMdd'));
t.the_time = string(datetime('now', 'Format', 'HHmm'));

% Generate unique seed for rng
rng('shuffle');
t.my_rng_seed = rng;

%% Input device name

p.device_string = 'Apple Internal Keyboard / Trackpad';
% p.device_string = 'USB-HID Keyboard'; 
% p.device_string ='Current Designs, Inc. 932';

%% Toggles

toggles.macOS = true; % If true, will skip psychtoolbox sync test 
toggles.gamma_correction = false;
toggles.save_textures = false;
toggles.save_run_info = false;

%% Set subject ID and number of runs

p.subj_ID = '000';
num_runs = 9; % default = 9

%% Set directories 

dirs.script_dir = pwd;
dirs.functions_dir = fullfile(dirs.script_dir, 'functions');
dirs.stimuli_dir = fullfile(dirs.script_dir, 'stimuli');
dirs.data_dir = fullfile(dirs.script_dir, 'data');
dirs.subj_dir = fullfile(dirs.data_dir, ['S' p.subj_ID]);
dirs.corrected_CLUT_dir = [];

if toggles.gamma_correction && ~exist(dirs.corrected_CLUT_dir, 'dir')
    error(['Corrected CLUT directory not found: ' dirs.corrected_CLUT_dir]);
end

if exist(dirs.functions_dir, 'dir')
    addpath(dirs.functions_dir);
else
    error(['Functions directory not found: ' dirs.functions_dir]);
end

if exist(dirs.stimuli_dir, 'dir')
    addpath(dirs.stimuli_dir);
else
    mkdir(dirs.stimuli_dir);
end

if exist(dirs.data_dir, 'dir')
    addpath(dirs.data_dir);
else
    mkdir(dirs.data_dir);
end

%% Check PTB

checkPTB();

%% Verify subject directory and run number

subj_files = dir([dirs.subj_dir, '/*.mat']);
subj_filenames = {subj_files.name};
if isempty(subj_filenames)
    if ~exist(dirs.subj_dir,'dir'), mkdir(dirs.subj_dir); end
    p.run_num = 1;
else
    most_recent_file = strsplit(subj_filenames{end}, '_');
    p.run_num = str2double(most_recent_file(end)) + 1;
end

%% Screen parameters

screens = Screen('Screens'); % Grab the available screens
w.use_screen = max(screens); % Use the most external monitor

w.view_distance = 57; % cm
w.screen_width = 30; % cm
w.screen_width_px = 756; 
w.screen_height_px = 491; 
w.centerX = round(w.screen_width_px/2); 
w.centerY = round(w.screen_height_px/2);

w.visual_angle = 2 * atan2d(w.screen_width/2,  w.view_distance); % deg of visual angle
w.ppd = round(w.screen_width_px/w.visual_angle); % px/deg
w.px_size = w.screen_width/w.screen_width_px; % cm/px

w.gray = 127;
w.white = ones(1,3) * 255; 
w.black = zeros(1,3); 

%% Prepare scan

t.TR = 1; % fMRI TR duration in seconds

[p, w, t, stimuli, frames] = prepareScan(p, w, t, dirs, toggles);

%% Present stimuli

for n_run = 1:num_runs

    disp(['Entering run ' num2str(n_run) '...']);

    run_info = presentStimuli(p, w, t, stimuli, frames);

    disp(['Run ' num2str(n_run) ' complete.']);
    
    if toggles.save_run_info
        save(fullfile(dirs.subj_dir, ['S' p.subj_ID '_run' num2str(n_run) '.mat']), 'run_info');
        disp(['Run ' num2str(n_run) ' saved.']);
    end
end

%% Restore setup and close screen

Screen('LoadNormalizedGammaTable', w.window, w.default_CLUT);
Screen('CloseAll'); ShowCursor;
PsychHID('KbQueueStop', p.device_number);
PsychHID('KbQueueRelease', p.device_number);