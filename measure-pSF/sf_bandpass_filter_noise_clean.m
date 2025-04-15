% Bandpass filter spatial frequency content of broadband, Gaussian white noise

%{

pSFT FBA Experiment
------------------------------

min SF level = 0.1 cpd
max SF level = 12 cpd
Bandpass filter width = 0.1 cpd
Logarithmically spaced samples = 40
Noise samples per SF = 10
Michelson Contrast = 100%
Size = full screen height

============================
Written by Luis D. Ramirez
luisdr@bu.edu
Boston University
Created on September 13, 2022

%}

%% Prepare workspace

% Clear workspace
clear all; close all; clc;
KbName('UnifyKeyNames');
commandwindow;

% Set directories
script_dir = pwd;
save_dir = 'probe_stimuli';

%% !! Toggles !! 

% Location  
which_setup = 2; % 0 = MacBook, 1 = Display++, 2 = Scanner

% Probe  
probe_side = 'R'; % Probe hemifield location: 'L' = left, 'R' = right
apply_aperture = 0; % Applies the wedge aperture to the image (!! turn this off if creating stimuli for the experiment !!)

% Outputs  
save_textures = 1; % Save .mat files of the bandpass filtered noise textures
plot_textures = 0; % Produces a figure for each unique texture (be mindful of this when generating large quantities of images)
cycle_textures = 0; % Note that this will infinitely cycle through all the noise textures. Press Ctrl + C to break this loop.
open_window = 0;

% Sync Test  
if which_setup == 0
    Screen('Preference', 'SkipSyncTests', 1); % 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

%% Display parameters
% Window selection: viewing distance, and width and height

% Display size & viewing distance 
if which_setup == 0 % If using my Macbook
    display_setup = 'Macbook';
    screens = Screen('Screens'); % Grab the available screens
    w.use_screen = min(screens); % If there are two or more displays, 'max' should grab the most external display.
    w.refresh_rate = 60; % display refresh rate in Hertz, Hz
    w.view_distance = 57; % in centimeters, cm
    w.screen_width = 30; % in centimeters, cm
    w.screen_width_px = 1024; % in pixels, px 
    w.screen_height_px = 640; % in pixels, px
    if open_window
        w.screen_width_px = 2048; % in pixels, px 
        w.screen_height_px = 1280; % in pixels, px 
    end
elseif which_setup == 1 % If using Display++ for staircasing
    display_setup = 'Display++';
    w.use_screen = 0;
    w.refresh_rate = 120; % display refresh rate in Hertz, Hz
    w.view_distance = 114;  % cm
    w.screen_width = 52.3; % cm
    w.screen_width_px = 1440; % px 
    w.screen_height_px = 1080; % px 
elseif which_setup == 2 % If using the scanner
    display_setup = 'Scanner';
    w.use_screen = 0;
    w.refresh_rate = 60; % display refresh rate in Hertz, Hz
    w.view_distance = 99;   % in cm, ideal distance: 1 cm equals 1 visual degree (at 57 cm)
    w.screen_width = 42.7;  % cm horizontal display size - 41.5 in scanner, 44.5 without eye-tracking; OLD SCREEN: 51.2
    w.screen_height = 32;   % cm 29.5 before 13/12/18, 29 old screen
    w.screen_width_px = 1024; % px
    w.screen_height_px = 768; % px
end

% Set screen resolution 
if which_setup == 0
%     old_res = Screen('Resolution', w.use_screen, w.screen_width_px, w.screen_height_px); % note: macOS 12 is incompatible with some Screen() resolution functions
else
%     old_res = Screen('ConfigureDisplay','Scanout', w.use_screen, w.use_screen, w.screen_width_px, w.screen_height_px);
end

% Visual angle & pixel info 
% Calculate visual angle of entire screen, pixels per degree of visual angle, and size of pixel in visual degrees
if which_setup ~= 2
    w.visual_angle = (2*atan2d(w.screen_width/2,  w.view_distance)); % Visual angle of the whole screen
    w.ppd = round(w.screen_width_px/w.visual_angle); % Pixels per degree of visual angle
    w.px_size = w.screen_width/w.screen_width_px; % size of pixel in cm
else % Use screen height if @ scanner to measure visual angle, pixel per degree, and pixel size
    w.visual_angle = (2*atan2d(w.screen_height/2,  w.view_distance)); 
    w.ppd = round(w.screen_height_px/w.visual_angle);
    w.px_size = w.screen_height/w.screen_height_px; 
end

% Gray value 
gray = 127;

% Nyquist frequency
% The nyquist frequency is 1 cycle in 2 pixels (0.5 cycles/pixel) – the highest frequency possible
f_Nyquist = 1/(2*w.px_size);  

%% Open window

if open_window
    bg_color = gray;
    [window, w.screen_rect_sz_px] = Screen('OpenWindow', w.use_screen, bg_color, [0 0 w.screen_width_px, w.screen_height_px]); %HideCursor;
    Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);   % Enable alpha blending
    w.centerX_pixels = w.screen_rect_sz_px(3)/2; w.centerY_pixels = w.screen_rect_sz_px(4)/2; % center X, Y coordinates in pixels

    % Text settings
    Screen('TextStyle', window, 1); % 0=normal, 1=bold, 2=italic
    Screen('TextSize', window, 18);

    loading_text_string = 'Loading stimuli...';
    loading_text_boundary = Screen('TextBounds', window, loading_text_string);
    loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX_pixels, w.centerY_pixels);
    Screen('DrawText', window, loading_text_string, loading_text_patch(1),  loading_text_patch(2), 255);
    Screen('Flip', window);
end

%% Noise stimulus parameters

% Probe contrast
noise_contrast = 1;

% Noise sample count
% default = 10
noise_sample_count = 10;

% Probe size
noise_height_px = w.screen_height_px;
noise_height_deg = noise_height_px/w.ppd;

if which_setup == 2
    noise_width_px = noise_height_px;
else
    noise_width_px = w.screen_width_px;
end

noise_width_deg = noise_width_px/w.ppd;

% Make probe size even 
if mod(noise_width_px,2), noise_width_px=noise_width_px+1; end
if mod(noise_height_px,2), noise_height_px=noise_height_px+1; end

%% == Bandpass filter parameters ==

% Filter count
% default = 40
sf_bp_filter_count = 40; 

% Probe minimum SF
sf_bp_filter_min = 0.1; 

% Probe max SF
sf_bp_filter_max = 12; 

% Bandpass filter width
sf_bp_filter_width = 0.1; 

% Gaussian SD 
sf_bp_filter_gauss_sd = 0.1*w.ppd; % used to smooth the bandpass filters

% Filter centers, lower and upper bounds
% The central spatial frequencies of the bandpass filters are logarithmically sampled from a given min. and max.
sf_bp_filter_centers = 10.^linspace(log10(sf_bp_filter_min), log10(sf_bp_filter_max), sf_bp_filter_count);
% Lower and upper bounds of each filter is its center ± the desired filter width
sf_bp_filter_lower_bound = sf_bp_filter_centers - sf_bp_filter_width;
sf_bp_filter_upper_bound = sf_bp_filter_centers + sf_bp_filter_width;

%% Make the probe aperture
% This aperture is applied to the texture itself and is not a separate mask.
% Alpha level for aperture: 
% 0 = completely transparent
% 255 = completely opaque

aperture_angle_cap = 100;

[x_cartesian_space, y_cartesian_space] = meshgrid(-(noise_width_px/2):(noise_width_px/2)-1, -(noise_height_px/2):(noise_height_px/2)-1); % create coordinates for texture: -x to x, -y to y;
probe_radii_coord = sqrt(x_cartesian_space.^2 + y_cartesian_space.^2); % Calculating the eccentricity (radius) of each point in the meshgrid relative to the center of the 2D image
aperture_angle = nan(noise_height_px, noise_width_px);
for x_cord = 1:noise_height_px
    for y_cord = 1:noise_width_px
        aperture_angle(x_cord, y_cord) = rad2deg(cart2pol(x_cartesian_space(x_cord, y_cord), y_cartesian_space(x_cord, y_cord)));
    end
end

aperture = zeros(noise_height_px, noise_width_px);
aperture(abs(aperture_angle) >= aperture_angle_cap) = 1;
% aperture(probe_radii_coord <= aperture_diam_px/2 & abs(aperture_angle) <= aperture_angle_cap) = 1;

% Apply gaussian filter to smooth aperture
gauss_filt_sd = 0.1*w.ppd;
aperture = imgaussfilt(aperture,gauss_filt_sd);
aperture = aperture(1:noise_height_px, 1:noise_width_px);
if strcmp(probe_side,'R')
    aperture = fliplr(aperture);
end

% apply an alpha layer to gray box that is fully transparent on the probe side
% aperture_stim(:,:,1) = ones(size(aperture))*gray;
% aperture_stim(:,:,2) = imgaussfilt(aperture*255, gauss_filt_sd);
% aperture_stim_made = Screen('MakeTexture', window, aperture_stim);

%% !! Create bandpass filtered noise !!

noise_textures = nan(noise_height_px, noise_height_px, sf_bp_filter_count, noise_sample_count);

for n_bp_filter = 1:sf_bp_filter_count
    for n_noise_sample = 1:noise_sample_count
        
        init_noise_tex = 2*rand(noise_height_px,noise_height_px)-1; % Create noise
        fft_noise = fftshift(fft2(init_noise_tex));   % FFT noise and shift 0Hz component to the center

        sf_bp_filter = Bandpass2(size(init_noise_tex), sf_bp_filter_lower_bound(n_bp_filter)/f_Nyquist, sf_bp_filter_upper_bound(n_bp_filter)/f_Nyquist); % Create the bandpass filter
        sf_bp_filter = imgaussfilt(sf_bp_filter, sf_bp_filter_gauss_sd); % Smooth bandpass filter
        sf_bp_filter_norm = (sf_bp_filter - min(sf_bp_filter(:))) ./ (max(sf_bp_filter(:)) - min(sf_bp_filter(:))); % Normalize filter

        filtered_noise = sf_bp_filter_norm .* fft_noise; % Apply the bandpass filter to FFT noise
        filtered_noise = ifft2(ifftshift(filtered_noise)); % Inverse shift and inverse fft to get the real image
        filtered_noise = (filtered_noise - min(filtered_noise(:)))./max(filtered_noise(:)); % Normalize the texture
        filtered_noise = noise_contrast * filtered_noise;

        if open_window
            filtered_noise = abs(filtered_noise*gray + gray); % Convert real component of texture to grayscale
            if apply_aperture, filtered_noise(:,:,2) = aperture.*255; end % apply aperture layer
            noise_textures_made(n_bp_filter, n_noise_sample) = Screen('MakeTexture', window, filtered_noise);
        else
            if apply_aperture, filtered_noise = filtered_noise.*aperture; end % apply aperture layer
            noise_textures(:,:,n_bp_filter, n_noise_sample) = abs(filtered_noise*gray+gray); % convert real component of texture to grayscale
        end
    end
end

if plot_textures
    for n_bp_filter = 1:sf_bp_filter_count
        figure('Name',['Center SF: ' num2str(round(sf_bp_filter_centers(n_bp_filter),2)) 'cpd'],'Color',[1 1 1])
        current_tex = noise_textures(:,:,n_bp_filter, 1);
        imshow(current_tex,[0 255]);
    end
end

if cycle_textures
    figure
    while 1 % press Ctrl + C to break this loop
        for n_bp_filter = 1:sf_bp_filter_count
            for n_noise_sample = 1:noise_sample_count
                current_tex = noise_textures(:,:,n_bp_filter, n_noise_sample);
                imshow(current_tex,[0 255]); title(['Center SF: ' num2str(sf_bp_filter_centers(n_bp_filter)) 'cpd']);
            end
        end
    end
end

%% Initialize device input
% Check which device number the keyboard is assigned to

%{
device_number = 0;
[keyboard_indices, product_names, ~] = GetKeyboardIndices;

if which_setup == 0 % If using my Macbook
    device_string = 'Apple Internal Keyboard / Trackpad'; % Macbook
    %     device_string =  'USB-HID Keyboard'; % external keyboard
    keypress_numbers = [KbName('LeftArrow') KbName('RightArrow')];

elseif which_setup == 1 % If using Display++ for staircasing
    device_string = 'Logitech USB Keyboard';
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    keypress_numbers = [KbName('LeftArrow') KbName('RightArrow')];

elseif which_setup == 2 % If scanning
%     device_string = {'Current Designs, Inc. 932'};
    device_string = 'AT Translated Set 2 keyboard';    % Prepare Psychtoolbox for Imaging %
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % Assign trigger key
    trigger_key = KbName('5%');  % KbName('=+');
    keypress_numbers = [KbName('1!') KbName('2@')];
end

for i = 1:length(product_names)
    if strcmp(product_names{i}, device_string)
        device_number = keyboard_indices(i);
        break;
    end
end

if device_number == 0
    error('No device by that name was detected');
end
%}

%% Draw to window

if open_window
    
    fixation_dot_degrees = 0.15;
    fixation_dot_pixels = fixation_dot_degrees*w.ppd;
    
    fixation_patch = CenterRectOnPoint([0 0 fixation_dot_pixels fixation_dot_pixels], w.centerX_pixels, w.centerY_pixels); % defining size and location of probe stimuli
    
    % Probe stimulus eccentricity
    p.probe_eccen_px = noise_width_px/2;

    if strcmp(probe_side,'L') % Probe will be on the left
        probe_patch = CenterRectOnPoint([0 0 noise_width_px noise_height_px], w.centerX_pixels-p.probe_eccen_px, w.centerY_pixels); % defining size and location of probe stimuli
        params_text_side = 'right';
    elseif strcmp(probe_side,'R') % Probe will be on the right
        probe_patch = CenterRectOnPoint([0 0 noise_width_px noise_height_px], w.centerX_pixels+p.probe_eccen_px, w.centerY_pixels); % defining size and location of probe stimuli
        params_text_side = 'left';
    end

    % Initialize keyboard input
    esc_pressed = 0;
    PsychHID('KbQueueCreate', device_number);
    PsychHID('KbQueueStart', device_number);
    
    % Infinite loop
    while ~esc_pressed
        for n_bp_filter = 1:sf_bp_filter_count
            for n_noise_sample = 1:noise_sample_count
                
                % Draw and Flip 
                Screen('DrawTexture', window, noise_textures_made(n_bp_filter, n_noise_sample),[],probe_patch);
                Screen('FillOval', window, 0, fixation_patch)
                DrawFormattedText(window, ...
                    ['sf = ' num2str(sf_bp_filter_centers(n_bp_filter)) ' cpd (' num2str(n_bp_filter) ')' '\n' '\n' 'sample # = ' num2str(n_noise_sample)], ...
                    params_text_side, 'center');
                
                Screen('Flip', window);
             

                % Monitor response 
                [key_pressed, first_press] = PsychHID('KbQueueCheck', device_number);
                which_press = find(first_press);

                if key_pressed
                    if which_press(1) == KbName('Escape')
                        esc_pressed = 1;
                    end
                end
                
                % Simulate experiment timing (10 Hz)
                pause(0.1);

            end
        end
    end
    
    % Clear screen
    sca;
end

%% !! Save textures !!

if save_textures

    % Store noise textures info into structure
    img.contrast = noise_contrast; % probe contrast
    img.min_SF = sf_bp_filter_min; % probe min. SF
    img.max_SF = sf_bp_filter_max; % probe max. SF
    img.filter_count = sf_bp_filter_count; % # of SF bandpass filters
    img.filters = sf_bp_filter_centers; % SF bandpass filter centers
    img.filter_width = sf_bp_filter_width; % SF bandpass filter width
    img.gauss_sd_width = sf_bp_filter_gauss_sd; % Gaussian filter width (used to smooth bandpass filter)
    img.noise_sample_count = noise_sample_count; % Noise sample count
    img.size_px = [noise_width_px noise_height_px]; % Noise size (pixels)
    img.size_deg = [noise_width_deg noise_height_deg]; % Noise size (visual degrees)
    img.ppd = w.ppd; % Pixels per degree of visual angle
    img.view_dist = w.view_distance; % Viewing distance (centimeters)

    % Set filename
    cd(save_dir)
    if apply_aperture
        filtered_textures_filename = ['filtered_noise_tex_' display_setup '_' num2str(noise_width_px) 'pxWidth_' num2str(100*noise_contrast) '%Contrast_' probe_side '.mat'];
    else
        filtered_textures_filename = ['filtered_noise_tex_' display_setup '_' num2str(noise_width_px) 'pxWidth_' num2str(100*noise_contrast) '%Contrast.mat'];
    end
    
    % Save matfile 
    save(filtered_textures_filename,'noise_textures','img','-v7.3'); disp([filtered_textures_filename ' saved in ' save_dir '.']);
    cd(script_dir)

end
