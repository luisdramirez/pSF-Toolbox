
%% Prepare workspace

clear all;
close all;
clc;
commandwindow;

% Grab date
t.the_date = string(datetime('now', 'Format', 'yyyyMMdd'));
t.the_time = string(datetime('now', 'Format', 'HHmm'));

% Generate unique seed for rng
rng('shuffle');
t.my_rng_seed = rng;

%% Input device name

p.device_string = 'Apple Internal Keyboard / Trackpad';

%% Toggles

toggles.macOS = true; % If true, will skip psychtoolbox sync test
toggles.gamma_correction = false;
toggles.save_textures = true;
toggles.save_run_info = true;

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

if ~exist(dirs.subj_dir,'dir')
    mkdir(dirs.subj_dir);
    p.run_num = 1;
else
    % Find existing run files (exclude _I.mat files)
    subj_files = dir(fullfile(dirs.subj_dir, ['S' p.subj_ID '_run*.mat']));
    subj_filenames = {subj_files.name};
    main_files = subj_filenames(~contains(subj_filenames, '_I.mat'));

    if isempty(main_files)
        p.run_num = 1;
    else
        % Extract run numbers from filenames
        run_nums = zeros(1, length(main_files));
        for i = 1:length(main_files)
            tokens = regexp(main_files{i}, '_run(\d+)\.mat$', 'tokens');
            if ~isempty(tokens)
                run_nums(i) = str2double(tokens{1}{1});
            end
        end
        p.run_num = max(run_nums) + 1;
    end
end

if p.run_num > num_runs
    disp(['All ' num2str(num_runs) ' runs already completed for subject S' p.subj_ID '.']);
    return;
end

disp(['Starting at run ' num2str(p.run_num) ' of ' num2str(num_runs) '.']);

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

w.gray = 127;
w.white = ones(1,3) * 255;
w.black = zeros(1,3);

%% Prepare scanning session

[p, w, t, stimuli, frames] = prepareScan(p, w, t, dirs, toggles);

%% Present stimuli

for n_run = p.run_num:num_runs

    disp(' ');
    disp(['Entering run ' num2str(n_run) ' of ' num2str(num_runs) '...']);

    run_info = presentStimuli(p, w, t, stimuli, frames);

    % Check if run was aborted
    if run_info.aborted
        disp('Run aborted by user. Exiting session.');
        break;
    end

    disp(['Run ' num2str(n_run) ' complete.']);

    %% Save run info and SF time series

    if toggles.save_run_info
        save(fullfile(dirs.subj_dir, ['S' p.subj_ID '_run' num2str(n_run) '.mat']), 'run_info');

        I = generateSFTimeSeries(run_info.p, run_info.t, run_info.filters);
        save(fullfile(dirs.subj_dir, ['S' p.subj_ID '_run' num2str(n_run) '_I.mat']), 'I');

        disp(['Run ' num2str(n_run) ' saved.']);
    end

end

%% Restore setup and close screen

Screen('LoadNormalizedGammaTable', w.window, w.default_CLUT);
Screen('CloseAll'); ShowCursor;
PsychHID('KbQueueStop', p.device_number);
PsychHID('KbQueueRelease', p.device_number);

disp('Session complete.');

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