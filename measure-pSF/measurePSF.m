% measurePSF

%% Prepare workspace

clear all; close all; clc;
commandwindow;

% Grab date
t.the_date = string(datetime('now', 'Format', 'yyyyMMdd'));
t.the_time = string(datetime('now', 'Format', 'HHmm'));

% Generate unique seed for random number generator (rng)
rng('shuffle');
t.my_rng_seed = rng;

% Set directories
project_dir = pwd;
stimuli_dir = 'stimuli';

addpath([project_dir '/functions']);
addpath([project_dir '/' stimuli_dir]);

%% Screen parameters

w.use_screen = 0;
w.view_distance = 57;
w.screen_width = 30;
w.screen_width_px = 800;
w.screen_height_px = 600;

w.visual_angle = (2*atan2d(w.screen_height/2,  w.view_distance));
w.ppd = round(w.screen_height_px/w.visual_angle);
w.px_size = w.screen_height/w.screen_height_px;

w.gray = 127;
w.bg_color = w.gray;

%% Load stimuli

texture_filename = 'bandpass_filtered_noise.mat';

cd([project_dir '/' stimuli_dir])
if exist(texture_filename, 'file') == 2
    load(texture_filename);
else
    textures = makeStimuli(image_size_px, ppd, px_size);
end
cd(project_dir)

%% Create aperture textures

aperture_radius_deg = w.visual_angle/2;
aperture_radius_px = round(aperture_radius_deg * w.ppd);

% Create cartesian coordinates for texture: -x to x, -y to y
[x,y] = meshgrid(-w.screen_width_px/2:w.screen_width_px/2-1,-w.screen_height_px/2:w.screen_height_px/2-1);

% Calculate eccentricity (e.g., radius) of each pixel coordinate in the cartesian grid, from the center of the grid
r = sqrt((x).^2 + (y).^2); % solving for r in: x^2 + y^2 = r^2

% Create a matrix of 1s
aperture_mask = ones(w.screen_height_px, w.screen_width_px);

% Set values inside desired eccentricity to be fully transparent
aperture_mask(r <= aperture_radius_px) = 0;

% Create aperture texture
aperture_tex = nan(w.screen_height_px, w.screen_width_px, 2);
aperture_tex(:,:,1) = ones(size(aperture_mask)) * w.gray;
aperture_tex(:,:,2) = imgaussfilt(aperture_mask,  0.1*w.ppd)*255;



%% Open window

[window, w.screen_size_px] = PsychImaging('OpenWindow', w.use_screen, w.bg_color, [0 0 w.screen_width_px w.screen_height_px]);
t.monitor_refresh_dur = Screen('GetFlipInterval', window);
HideCursor;

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Make stimuli

% Fixation space

% Make aperture
aperture = Screen('MakeTexture', window, aperture_tex);

% pSF stimuli
pSF_stimuli = nan(filter.count, stimulus.noise_sample_count);
for n_filter = 1:filter.count
    for n_sample = 1:stimulus.noise_sample_count
        pSF_stimuli(n_filter, n_sample) = Screen('MakeTexture', window, textures(:, :, n_filter, n_sample));
    end
end

%% Generate presentation sequences

p.num_blocks = 6;

p.sf_pres_order = nan(p.num_blocks, filter.count);
for n_block = 1:p.num_blocks
    p.sf_pres_order(n_block, :) = datasample(1:filter.count, filter.count, 'Replace', false); 
end

%% Stimulus presentation

for n_block = 1:p.num_blocks




end

%% Measure behavioral performance

behav_data.detection_accuracy = [];

%% Save data

if toggles.save_run

    disp('Saving run...')
    
    cd(data_dir)
    scan_run.p = p;
    scan_run.t = t;
    scan_run.w = w;
    scan_run.frames = frames;
    scan_run.pres_timing = pres_timing;
    scan_run.exe_time = exe_time;
    scan_run.behav_data = behav_data;

    save_filename = ['measurePSF_S' p.subj_ID '_Run' num2str(p.run_num) '.mat'];
    save(save_filename, 'scan_run');    
    disp(['Run ' num2str(p.run_num) ' saved!']);

end