% verifyStimuli
%   Verifies the textures by computing the power spectrum.
%   Plots power spectra and textures, and saves the figures.
%   Can also generate and save the textures if they do not exist.
%
% Syntax
%   verifyStimuli()
%
% Input Arguments
%   None
% Output Arguments
%   None

function verifyStimuli()

%% Toggles

toggles.save_textures = true;
toggles.save_energy = true;
toggles.save_texture_figures = true;
toggles.save_energy_figures = true;

%% Set directories

script_dir = fileparts(mfilename('fullpath')); % Updated to be robust
stimuli_dir = script_dir;
parent_dir = fileparts(script_dir);
functions_dir = fullfile(parent_dir, 'functions');
texture_figures_dir = fullfile(stimuli_dir, 'figures', 'textures');
energy_figures_dir = fullfile(stimuli_dir, 'figures', 'energy');

% Add functions to path
if exist(functions_dir, 'dir')
    addpath(functions_dir);
else
    error(['Functions directory not found: ' functions_dir]);
end

% Create texture figures directory
if ~exist(texture_figures_dir, 'dir')
    mkdir(texture_figures_dir);
end

% Create energy figures directory
if ~exist(energy_figures_dir, 'dir')
    mkdir(energy_figures_dir);
end

%% Screen parameters
% Note: If the parameters here match your experimental setup, then this script can generate the textures for your setup.

w.view_distance = 57; % cm
w.screen_width = 30; % cm
w.screen_width_px = 756;
w.screen_height_px = 491;
w.centerX = round(w.screen_width_px/2);
w.centerY = round(w.screen_height_px/2);

w.visual_angle = 2 * atan2d(w.screen_width/2,  w.view_distance); % degrees of visual angle (deg)
w.ppd = round(w.screen_width_px/w.visual_angle); % pixels per degree (px/deg)
w.px_size = w.screen_width/w.screen_width_px; % cm per pixel (cm/px)

%% Stimuli parameters

p = struct();
p = stimulusParams(p, w);

%% Load textures if they exist, otherwise generate them

% Construct filename to verify (matches prepareScan.m logic)
% Example output filename: bandpass_filtered_noise_42_123_090_40_10.mat
texture_filename = sprintf('bandpass_filtered_noise_%s_%s_%s_%s_%s.mat', ...
    strrep(sprintf('%0.0f', w.ppd), '.', ''), ...
    strrep(sprintf('%0.0f', p.stimulus_diameter_px), '.', ''), ...
    strrep(sprintf('%0.2f', p.stimulus_contrast), '.', ''), ...
    strrep(sprintf('%0.0f', p.noise_filter_count), '.', ''), ...
    strrep(sprintf('%0.0f', p.noise_sample_count), '.', ''));

texture_filepath = fullfile(stimuli_dir, texture_filename);

fprintf('Target texture file to load: %s\n', texture_filename);

if ~exist(texture_filepath, 'file')
    response = input([texture_filepath ' not found. Would you like to create it? (y/n) '], 's');
    if response == 'y'
        [textures, filters] = createTextures(p.stimulus_diameter_px, p.stimulus_contrast, p.noise_filter_count, p.noise_sample_count, w.ppd, toggles.save_textures, texture_filepath);
    else
        error('Texture file not found and user did not want to create it.');
    end
else
    load(texture_filepath);
    disp('Textures loaded.')
end

%% Check if energy profiles file exists

% Construct filename similar to texture filename
% Example output filename: energy_42_123_090_40_10.mat
energy_filename = sprintf('energy_%s_%s_%s_%s_%s.mat', ...
    strrep(sprintf('%0.0f', w.ppd), '.', ''), ...
    strrep(sprintf('%0.0f', p.stimulus_diameter_px), '.', ''), ...
    strrep(sprintf('%0.2f', p.stimulus_contrast), '.', ''), ...
    strrep(sprintf('%0.0f', p.noise_filter_count), '.', ''), ...
    strrep(sprintf('%0.0f', p.noise_sample_count), '.', ''));

energy_filepath = fullfile(stimuli_dir, energy_filename);

fprintf('Target energy file: %s\n', energy_filename);

energy_file_exists = exist(energy_filepath, 'file');
if energy_file_exists
    fprintf('Energy file found.\n');

    % Continue asking until a valid response is given
    valid_response = false;
    while ~valid_response
        response = input('Overwrite (ow), Load (l), or Exit (x)? ', 's');

        if strcmpi(response, 'x')
            fprintf('Exiting...\n');
            return;
        elseif strcmpi(response, 'l')
            % Load existing energy data
            load(energy_filepath, 'energy_data');
            all_energy_profiles = energy_data.energy_profiles;
            all_bin_centers = energy_data.bin_centers;
            all_target_centers = energy_data.target_centers;
            all_peak_freqs = energy_data.peak_freqs;
            n_textures = energy_data.n_textures;

            fprintf('Energy profiles loaded.\n');

            % Ask about regenerating figures
            regen_textures = input('Regenerate texture figures? (y/n) ', 's');
            regen_individual = input('Regenerate individual profile figures? (y/n) ', 's');
            regen_all = input('Regenerate all profile figure? (y/n) ', 's');

            % Set toggles based on user input
            if strcmpi(regen_textures, 'y')
                toggles.save_texture_figures = true;
            else
                toggles.save_texture_figures = false;
            end

            if strcmpi(regen_individual, 'y')
                toggles.save_energy_figures = true;
                regenerate_individual = true;
            else
                toggles.save_energy_figures = false;
                regenerate_individual = false;
            end

            if strcmpi(regen_all, 'y')
                regenerate_all_profiles = true;
            else
                regenerate_all_profiles = false;
            end

            % Skip verification loop
            skip_verification = true;
            valid_response = true;
        elseif strcmpi(response, 'ow')
            % Overwrite - proceed as normal
            fprintf('Proceeding with verification (will overwrite existing file).\n');
            skip_verification = false;
            regenerate_individual = false;
            regenerate_all_profiles = false;
            valid_response = true;
        else
            fprintf('Invalid response. Please enter "ow" to overwrite, "l" to load, or "x" to exit.\n');
        end
    end
else
    fprintf('Energy file not found. Proceeding with verification.\n');
    skip_verification = false;
    regenerate_individual = false;
    regenerate_all_profiles = false;
end

%% Prepare figure handle

if toggles.save_texture_figures || toggles.save_energy_figures
    fg = figure('Visible','on','Color',[1 1 1]);
    set(0,'CurrentFigure',fg);
    plot_settings = plotSettings();
end

%% Run texture energy verification

if ~skip_verification
    fprintf('Starting texture verification...\n');
    verification_start = tic;

    n_textures = size(textures,3);

    % Initialize storage for all energy profiles
    all_energy_profiles = cell(n_textures, 1);
    all_bin_centers = cell(n_textures, 1);
    all_target_centers = zeros(n_textures, 1);
    all_peak_freqs = zeros(n_textures, 1);

    for sf_idx = 1:n_textures

        fprintf('Processing texture %d of %d... \n', sf_idx, n_textures);
        loop_start = tic;

        %% Pull texture

        texture = textures(:,:,sf_idx,1);

        % Save texture figure
        if toggles.save_texture_figures
            imagesc(texture), colormap(gray), axis square, box off, axis off;

            figure_name = ['Texture ' num2str(sf_idx)];
            figure_path = fullfile(texture_figures_dir, [figure_name '.pdf']);
            saveas(gcf, figure_path);

            idx = strfind(figure_path, '/measure-pSFT/');
            if ~isempty(idx)
                disp_path = figure_path(idx(1):end);
            else
                disp_path = figure_path;
            end
            disp(['Saved ' figure_name ' in ' disp_path]);
            clf
        end

        %% Compute power spectrum
        % FFT calculates energy at fixed intervals of ppd/stimulus_diameter_px.
        % If stimulus_diameter_px doesn't contain an integer number of cycles, the target frequency will fall "between" two bins.
        % The nearest bin is chosen, resulting in offsets.
        % To avoid this, the frequency domain is padded for interpolation.

        % Subtract mean to remove DC component
        texture = texture - mean(texture(:));

        % Calculate the native frequency resolution (cycles per degree per bin)
        native_res = w.ppd / size(texture, 1);

        % If target frequency is not a multiple of resolution, or resolution is too coarse
        if mod(filters.centers(sf_idx), native_res) > 1e-5 || size(texture, 1) < 512
            % Pad to a high-power of 2 to interpolate the spectrum for peak precision
            padding = 2^nextpow2(size(texture, 1) * 2);
        else
            padding = size(texture, 1); % Stay in native resolution
        end

        fft_texture = fftshift(fft2(texture, padding, padding)); % convert to frequency domain and shift zero frequency to center
        power_spectrum = abs(fft_texture).^2; % compute power spectrum by squaring the magnitude

        %% Create coordinate space

        [rows, cols] = size(fft_texture);

        % Frequency axes range from -0.5 to just below 0.5 (Nyquist).
        % The +0.5 and -0.5 frequencies are aliases, so only -0.5 is included.

        u = (-cols/2:cols/2-1) / cols;
        v = (-rows/2:rows/2-1) / rows;
        [U, V] = meshgrid(u, v);

        %% Compute radial frequency

        rho_cpp = sqrt(U.^2 + V.^2); % cycles per pixel
        rho_cpd = rho_cpp * w.ppd; % (cycles/pixel) * (pixels/degree) = (cycles/degree)

        %% Get target center frequency

        target_center = filters.centers(sf_idx);

        %% Radially average power spectrum

        max_freq = max(rho_cpd(:));

        % Use overlapping bins: fine step size for smooth curve, wider bin for averaging
        step_sz = 0.001;    % Fine step between bin centers (controls curve resolution)
        bin_width = 0.05;   % Width of each bin (controls smoothing amount)

        % Discretize rho_cpd into fine, non-overlapping bins (size = step_sz)
        fine_edges = 0:step_sz:max_freq+step_sz;
        [~, ~, bin_idx] = histcounts(rho_cpd(:), fine_edges);

        % Accumulate power and counts in bins
        % Filter out invalid bins (0)
        valid_mask = bin_idx > 0;
        bin_idx = bin_idx(valid_mask);
        power_vals = power_spectrum(valid_mask);

        num_fine_bins = length(fine_edges) - 1;
        fine_power = accumarray(bin_idx, power_vals, [num_fine_bins 1]);
        fine_counts = accumarray(bin_idx, 1, [num_fine_bins 1]);

        % Apply sliding window smoothing (boxcar average)
        window_sz = round(bin_width / step_sz);

        smoothed_power = movsum(fine_power, window_sz);
        smoothed_counts = movsum(fine_counts, window_sz);

        % Compute mean energy profile
        energy_profile = smoothed_power ./ smoothed_counts;

        % Generate bin centers
        bin_centers = (fine_edges(1:end-1) + fine_edges(2:end)) / 2;

        % Normalize energy profile
        energy_profile = energy_profile / max(energy_profile);

        %% Find peak frequency from energy profile and compute offset from target

        [~, peak_idx] = max(energy_profile);
        peak_freq = bin_centers(peak_idx);

        % Compute offsets
        percent_offset = 100 * (peak_freq - target_center) / target_center;
        octave_offset = log2(peak_freq / target_center);

        if abs(percent_offset) > 5
            warning('Peak SF fidelity failure: Texture %d has a peak offset of %.2f%%. Consider adjusting diameter_px.', sf_idx, percent_offset);
        end

        % Store energy profile data for combined plot
        all_energy_profiles{sf_idx} = energy_profile;
        all_bin_centers{sf_idx} = bin_centers;
        all_target_centers(sf_idx) = target_center;
        all_peak_freqs(sf_idx) = peak_freq;

        %% Plot power spectrum

        if toggles.save_energy_figures

            plot(bin_centers, energy_profile, 'k-', 'LineWidth', 1); % Measured energy profile
            hold on;

            xline(target_center, 'g-', 'LineWidth', 1.5); % Target center
            xline(peak_freq, 'r-', 'LineWidth', 1.5); % Actual peak

            % Format figure
            title(sprintf('Target: %.3f cpd | Peak: %.3f cpd | Offset: %+.3f%% (%+.3f oct)', target_center, peak_freq, percent_offset, octave_offset));
            xlabel('Spatial frequency (cpd)');
            ylabel('Normalized power');
            xticks([0.1 0.5 1 5 10 15]);
            xlim([0.1 15]);
            ylim([0 1]);
            set(gca, 'TickDir', 'out', 'XScale', 'log');
            box off;
            legend({'Energy Profile', 'Target', 'Peak'}, 'Location', 'best');

            % Save figure
            figure_name = ['Filter ' num2str(sf_idx) ' Power Spectrum'];
            figure_path = fullfile(energy_figures_dir, [figure_name '.pdf']);
            saveas(gcf, figure_path);

            idx = strfind(figure_path, '/measure-pSFT/');
            if ~isempty(idx)
                disp_path = figure_path(idx(1):end);
            else
                disp_path = figure_path;
            end
            disp(['Saved ' figure_name ' in ' disp_path]);
            clf

        end

        loop_duration = toc(loop_start);
        fprintf('Done (%.2f s).\n', loop_duration);

    end

    total_duration = toc(verification_start);
    fprintf('Texture verification complete (Total time: %.2f s).\n', total_duration);
end

%% Save energy profiles

if toggles.save_energy && ~skip_verification

    % Prepare data structure to save
    energy_data = struct();
    energy_data.energy_profiles = all_energy_profiles;
    energy_data.bin_centers = all_bin_centers;
    energy_data.target_centers = all_target_centers;
    energy_data.peak_freqs = all_peak_freqs;
    energy_data.n_textures = n_textures;
    energy_data.ppd = w.ppd;
    energy_data.stimulus_diameter_px = p.stimulus_diameter_px;
    energy_data.stimulus_contrast = p.stimulus_contrast;
    energy_data.noise_filter_count = p.noise_filter_count;
    energy_data.noise_sample_count = p.noise_sample_count;

    % Save energy profiles
    save(energy_filepath, 'energy_data', '-v7.3');

    idx = strfind(energy_filepath, '/measure-pSFT/');
    if ~isempty(idx)
        disp_path = energy_filepath(idx(1):end);
    else
        disp_path = energy_filepath;
    end
    fprintf('Saved energy profiles to: %s\n', disp_path);
end

%% Regenerate texture figures (if loading)

if skip_verification && toggles.save_texture_figures
    % Prepare figure handle
    if ~exist('fg', 'var') || ~isvalid(fg)
        fg = figure('Visible','on','Color',[1 1 1]);
        set(0,'CurrentFigure',fg);
    end

    % Load textures if needed
    if ~exist('textures', 'var')
        load(texture_filepath);
    end

    if ~exist('n_textures', 'var')
        n_textures = size(textures, 3);
    end

    for sf_idx = 1:n_textures
        fprintf('Regenerating texture figure %d of %d...\n', sf_idx, n_textures);

        texture = textures(:,:,sf_idx,1);
        imagesc(texture), colormap(gray), axis square, box off, axis off;

        figure_name = ['Texture ' num2str(sf_idx)];
        figure_path = fullfile(texture_figures_dir, [figure_name '.pdf']);
        saveas(gcf, figure_path);

        idx = strfind(figure_path, '/measure-pSFT/');
        if ~isempty(idx)
            disp_path = figure_path(idx(1):end);
        else
            disp_path = figure_path;
        end
        disp(['Saved ' figure_name ' in ' disp_path]);
        clf
    end
end

%% Regenerate individual profile figures (if loading)

if skip_verification && regenerate_individual
    % Prepare figure handle
    if ~exist('fg', 'var') || ~isvalid(fg)
        fg = figure('Visible','on','Color',[1 1 1]);
        set(0,'CurrentFigure',fg);
    end
    if ~exist('plot_settings', 'var')
        plot_settings = plotSettings();
    end

    % Load textures if needed for texture figures
    if toggles.save_texture_figures
        if ~exist('textures', 'var')
            load(texture_filepath);
        end
    end

    for sf_idx = 1:n_textures
        fprintf('Regenerating figure for texture %d of %d...\n', sf_idx, n_textures);

        % Plot individual energy profile
        energy_profile = all_energy_profiles{sf_idx};
        bin_centers = all_bin_centers{sf_idx};
        target_center = all_target_centers(sf_idx);
        peak_freq = all_peak_freqs(sf_idx);

        percent_offset = 100 * (peak_freq - target_center) / target_center;
        octave_offset = log2(peak_freq / target_center);

        plot(bin_centers, energy_profile, 'k-', 'LineWidth', 1);
        hold on;

        xline(target_center, 'g-', 'LineWidth', 1.5);
        xline(peak_freq, 'r-', 'LineWidth', 1.5);

        % Format figure
        title(sprintf('Target: %.3f cpd | Peak: %.3f cpd | Offset: %+.3f%% (%+.3f oct)', target_center, peak_freq, percent_offset, octave_offset));
        xlabel('Spatial frequency (cpd)');
        ylabel('Normalized power');
        xticks([0.1 0.5 1 5 10 15]);
        xlim([0.1 15]);
        ylim([0 1]);
        set(gca, 'TickDir', 'out', 'XScale', 'log');
        box off;
        legend({'Energy Profile', 'Target', 'Peak'}, 'Location', 'best');

        % Save figure
        figure_name = ['Filter ' num2str(sf_idx) ' Power Spectrum'];
        figure_path = fullfile(energy_figures_dir, [figure_name '.pdf']);
        saveas(gcf, figure_path);

        idx = strfind(figure_path, '/measure-pSFT/');
        if ~isempty(idx)
            disp_path = figure_path(idx(1):end);
        else
            disp_path = figure_path;
        end
        disp(['Saved ' figure_name ' in ' disp_path]);
        clf
    end
end

%% Plot all energy profiles on one axes

if toggles.save_energy_figures || (skip_verification && regenerate_all_profiles)

    % Prepare figure handle if needed
    if ~exist('fg', 'var') || ~isvalid(fg)
        fg = figure('Visible','on','Color',[1 1 1]);
        set(0,'CurrentFigure',fg);
    end
    if ~exist('plot_settings', 'var')
        plot_settings = plotSettings();
    end

    hold on;

    % Plot all energy profiles
    for sf_idx = 1:n_textures
        % Create a linear gradient from blue to red for all profiles
        color1 = plot_settings.colors.blue;
        color2 = plot_settings.colors.red;
        if n_textures > 1
            this_color = color1 + (color2 - color1) * ((sf_idx-1)/(n_textures-1));
        else
            this_color = color1;
        end
        plot(all_bin_centers{sf_idx}, all_energy_profiles{sf_idx}, '-', 'LineWidth', 1, 'Color', this_color);
    end

    % Format figure
    xlabel('Spatial frequency (cpd)');
    ylabel('Power');
    xticks([0.1 0.5 1 5 10 15]);
    xlim([0.1 15]);
    ylim([0 1]);
    set(gca, 'TickDir', 'out', 'XScale', 'log');
    box off;

    % Save figure as vector graphics
    figure_name = 'All Energy Profiles';
    figure_path = fullfile(energy_figures_dir, [figure_name '.pdf']);
    print(gcf, figure_path, '-dpdf', '-painters'); % Ensure vector graphics

    idx = strfind(figure_path, '/measure-pSFT/');
    if ~isempty(idx)
        disp_path = figure_path(idx(1):end);
    else
        disp_path = figure_path;
    end
    disp(['Saved ' figure_name ' in ' disp_path]);
end

% Close figure handle
if toggles.save_texture_figures || toggles.save_energy_figures
    close(fg)
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