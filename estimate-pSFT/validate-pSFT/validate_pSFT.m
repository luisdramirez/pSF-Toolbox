% validate_pSFT
%
% This script validates the pSFT estimation pipeline by simulating stimulus and BOLD time series from pre-defined pSFT parameters.
% It then estimates the pSFT parameters and compares the estimates to the ground truth.
% Results are saved under `validate-pSFT/estimates/validation_pSFT_n<num_subjs>_yyyy-mm-dd_hh-mm-ss/` and `validate-pSFT/figures/validation_pSFT_n<num_subjs>_yyyy-mm-dd_hh-mm-ss/`.

%% Prepare workspace

clear all;
close all;
clc;

curr_time = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));

%% Set directories

project_dir = pwd;
save_dir = 'estimates';
figure_dir = 'figures';

addpath('../functions'); % estimate-pSFT functions
addpath('functions');  % validate-pSFT functions

if ~exist(fullfile(project_dir, save_dir), 'dir')
    mkdir(fullfile(project_dir, save_dir));
end
if ~exist(fullfile(project_dir, figure_dir), 'dir')
    mkdir(fullfile(project_dir, figure_dir));
end

addpath(fullfile(project_dir, save_dir));
addpath(fullfile(project_dir, figure_dir));

%% Toggles

toggles.parallelization = true;
toggles.coarse_grid_search = true;
toggles.fine_grid_search = true;
toggles.disp_on = true; % Display progress

toggles.save_validation_results = true;
toggles.make_summary_plots = true;
toggles.save_summary_plots = true;
toggles.make_voxel_plots = true;
toggles.save_voxel_plots = true;

%% Check for required toolboxes

checkRequiredToolboxes(toggles);

%% Parallelization setup for parfor loop

p.num_cores = 8;
p.num_chunks = p.num_cores - 1;

if toggles.parallelization
    maxNumCompThreads(p.num_cores);
    if p.num_cores > 1
        parpool('local', p.num_chunks)
    end
end

%% Spatial frequency parameters
% Hardcoded SF parameters used for generating pSFT curves.
% Ideally, the range should cover the range of SFs presented in the experimental design.

p.sf_min = 0.5;
p.sf_max = 12;

p.sf_count = 100;
p.sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), p.sf_count);

%% Initialize pSFT model parameters

p.initial_params = [1 1 1 0]; % [mu, sigma, beta, beta_0]

p.pSFT_bounds(1,:) = [6, 4, 25, 10]; % upper bounds
p.pSFT_bounds(2,:) = [0.009, 0.1, -25, -10]; % lower bounds

p.fmincon_options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000, 'display', 'off');

%% Simulation parameters

p.num_subjs = 1;
p.num_ROIs = 1;
p.num_voxels_per_ROI = 1000;

roi_names = {'V1', 'V2', 'V3'};

% Experimental design parameters (matching measure-pSFT)
num_runs = 9;
num_blocks = 6;
TR = 1; % Repetition time (seconds)
blank_duration = 10; % Blank period duration (TRs)

% Stimulus spatial frequencies (logarithmically spaced, matching createTextures.m)
num_center_sfs = 40;
center_sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), num_center_sfs);

% Blank period SF (small value to avoid log(0) in logGauss)
blank_sf = 0.0001;

% Calculate timing
block_duration = num_center_sfs; % 40 TRs per block (1 SF per TR)
num_TRs_per_run = (block_duration * num_blocks) + (blank_duration * (num_blocks + 1)); % 310 TRs
num_TRs = num_TRs_per_run * num_runs; % 2790 total TRs

% Ground truth pSFT parameter ranges
gt_mu_range = [0.5 4];       % pSFT peak range (cpd)
gt_sigma_range = [0.3 1.5];  % pSFT linear bandwidth range (cpd)
gt_beta_range = [1 5];       % BOLD amplitude range
gt_beta0_range = [-0.5 0.5]; % BOLD baseline range

% BOLD noise parameters (SNR in dB)
% SNR_dB = 10 * log10(signal_variance / noise_variance)
% Default: Low noise: SNR = 5.29 dB; Mid noise: SNR = -0.51 dB; High noise: SNR = -4.29 dB (Lerma-Usabiaga et al. 2020 PLOS Computational Biology, Figure 4)
SNR_levels_dB = [5.29, -0.51, -4.29];  % Low, Mid, High noise
num_SNR_levels = length(SNR_levels_dB);
SNR_level_names = cell(1, num_SNR_levels);
for i = 1:num_SNR_levels
    SNR_level_names{i} = sprintf('SNR: %.2f dB', SNR_levels_dB(i));
end

%% Define hemodynamic response function (HRF)

HRFs = cell(p.num_subjs, p.num_ROIs);

for subj = 1:p.num_subjs
    for roi = 1:p.num_ROIs
        HRFs{subj, roi} = repmat(defineHRF(), 1, p.num_voxels_per_ROI);
    end
end

%% Generate stimulus time series for each subject
% Each subject sees a unique random SF order (as in presentStimuli.m)
% Structure per run: blank(10) + [block(40) + blank(10)] x 6

if toggles.disp_on, disp('Generating stimulus time series for each subject...'); end

I_all = cell(p.num_subjs, 1);

for subj = 1:p.num_subjs
    I_subj = [];

    for n_run = 1:num_runs
        I_run = [];

        for n_block = 1:num_blocks
            % Add blank period before each block
            I_run = [I_run; repmat(blank_sf, blank_duration, 1)];

            % Add block: randomly sample SF order (like presentStimuli.m)
            block_sf_order = center_sfs(randperm(num_center_sfs));
            I_run = [I_run; block_sf_order'];
        end

        % Add final blank period
        I_run = [I_run; repmat(blank_sf, blank_duration, 1)];
        I_subj = [I_subj; I_run];
    end

    I_all{subj} = I_subj;
end

if toggles.disp_on
    disp(['  ' num2str(p.num_subjs) ' subjects']);
    disp(['  ' num2str(num_runs) ' runs x ' num2str(num_TRs_per_run) ' TRs/run = ' num2str(num_TRs) ' TRs']);
    disp('Done!');
end

%% Define pSFT model

pSFT_model = @(params, sf) logGauss(params, sf);

%% Initialize simulated data storage

struct_size = cell(p.num_subjs, p.num_ROIs);

simulated_data = struct('params', struct_size, ...
    'SFT', struct_size, ...
    'R', struct_size, ...
    'BOLD', struct_size, ...
    'measured_BOLD', struct_size, ...
    'I', struct_size);

%% Generate ground truth pSFT parameters

if toggles.disp_on, disp('Generating ground truth pSFT parameters...'); end

for subj = 1:p.num_subjs
    for roi = 1:p.num_ROIs
        % Generate random ground truth parameters for each voxel
        gt_mu = gt_mu_range(1) + (gt_mu_range(2) - gt_mu_range(1)) * rand(1, p.num_voxels_per_ROI);
        gt_sigma = gt_sigma_range(1) + (gt_sigma_range(2) - gt_sigma_range(1)) * rand(1, p.num_voxels_per_ROI);
        gt_beta = gt_beta_range(1) + (gt_beta_range(2) - gt_beta_range(1)) * rand(1, p.num_voxels_per_ROI);
        gt_beta0 = gt_beta0_range(1) + (gt_beta0_range(2) - gt_beta0_range(1)) * rand(1, p.num_voxels_per_ROI);

        simulated_data(subj, roi).params = [gt_mu; gt_sigma; gt_beta; gt_beta0];

        % Generate ground truth SFT curves for visualization
        gt_SFT = nan(p.sf_count, p.num_voxels_per_ROI);
        for vox = 1:p.num_voxels_per_ROI
            gt_SFT(:, vox) = pSFT_model([gt_mu(vox), gt_sigma(vox)], p.sfs);
        end
        simulated_data(subj, roi).SFT = gt_SFT;

    end
end

if toggles.disp_on, disp('Done!'); end

%% Generate simulated BOLD time series

if toggles.disp_on, disp('Generating simulated BOLD time series...'); end

for subj = 1:p.num_subjs
    I = I_all{subj};
    for roi = 1:p.num_ROIs
        gt_params = simulated_data(subj, roi).params;
        HRF = HRFs{subj, roi};

        data = generateSimulatedBOLD(gt_params, I, HRF, pSFT_model, SNR_levels_dB);

        simulated_data(subj, roi).I = data.I;
        simulated_data(subj, roi).R = data.R;
        simulated_data(subj, roi).BOLD = data.BOLD;
        simulated_data(subj, roi).measured_BOLD = data.measured_BOLD;
        simulated_data(subj, roi).SNR_dB = data.SNR_dB;
    end
end

if toggles.disp_on, disp('Done!'); end

%% Initialize pSFT struct (for each SNR level)

struct_size = cell(p.num_subjs, p.num_ROIs);

pSFT_template = struct('vox_indices', struct_size, ...
    'measured_BOLD', struct_size, ...
    'HRF', struct_size, ...
    'param_est', struct_size, ...
    'est_SFT', struct_size, ...
    'est_R', struct_size, ...
    'est_BOLD', struct_size, ...
    'r2', struct_size, ...
    'sse', struct_size, ...
    'start_values', struct_size, ...
    'start_sse', struct_size, ...
    'exitflag', struct_size);

% Initialize cell array for each SNR level
all_pSFT = cell(1, num_SNR_levels);
for snr_idx = 1:num_SNR_levels
    all_pSFT{snr_idx} = pSFT_template;
end

%% Estimate pSFT for each SNR level

total_elapsed_time = 0;

for snr_idx = 1:num_SNR_levels

    if toggles.disp_on
        disp(' ');
        disp(['====== SNR Level: ' SNR_level_names{snr_idx} ' ======']);
    end

    for subj = 1:p.num_subjs

        if toggles.disp_on, disp(['+++ S' num2str(subj) ' +++']); end

        for roi = 1:p.num_ROIs

            if toggles.disp_on, disp(['- ' roi_names{roi} ' -']); end

            tic;

            I = simulated_data(subj, roi).I;
            measured_BOLD = simulated_data(subj, roi).measured_BOLD{snr_idx};
            HRF = HRFs{subj, roi};

            % Enter estimatePSFT
            pSFT = estimatePSFT(I, measured_BOLD, HRF, p, toggles);
            all_pSFT{snr_idx}(subj, roi) = pSFT;

            elapsed_time = round(toc / 60, 1);
            if toggles.disp_on, disp(['Elapsed time: ~' num2str(elapsed_time) ' minute(s).']); disp(' '); end
            total_elapsed_time = total_elapsed_time + elapsed_time;

        end
    end
end

if toggles.disp_on, disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']); end
if toggles.parallelization, delete(gcp); end

%% Define validation filename

validation_filename = ['validation_pSFT_n' num2str(p.num_subjs) '_' char(curr_time)];

%% Compute validation metrics for each SNR level

[validation_metrics, pooled_data, simulated_data, all_pSFT] = computeValidationMetrics(simulated_data, all_pSFT, SNR_levels_dB, SNR_level_names, roi_names, p, toggles.disp_on);

%% Save validation results

if toggles.save_validation_results
    save([save_dir '/' validation_filename '.mat'], 'all_pSFT', 'simulated_data', 'validation_metrics', 'pooled_data', 'SNR_levels_dB', 'SNR_level_names', 'p');
    if toggles.disp_on, disp(['Saved ' validation_filename ' in /' save_dir]); end
end

%% Plot settings

plot_settings = plotSettings();

% Generate colors for SNR levels
key_colors = [
    plot_settings.colors.green;      % Green
    plot_settings.colors.yellow;      % Yellow
    plot_settings.colors.orange;      % Orange
    plot_settings.colors.red         % Red
    ];
if num_SNR_levels == 1
    snr_colors = {key_colors(1, :)};
else
    color_indices = linspace(1, size(key_colors, 1), num_SNR_levels);
    snr_colors = cell(1, num_SNR_levels);
    for i = 1:num_SNR_levels
        idx = color_indices(i);
        idx_low = floor(idx);
        idx_high = min(ceil(idx), size(key_colors, 1));
        if idx_low == idx_high
            snr_colors{i} = key_colors(idx_low, :);
        else
            t = idx - idx_low;
            snr_colors{i} = (1 - t) * key_colors(idx_low, :) + t * key_colors(idx_high, :);
        end
    end
end

%% Plot voxel fits

figure_save_dir = fullfile(figure_dir, validation_filename);
if ~exist(figure_save_dir, 'dir'), mkdir(figure_save_dir); end

if toggles.make_voxel_plots

    num_voxels_to_plot = 1;
    snr_idx_for_voxel_plots = 1;
    num_runs_to_plot = 1;

    fg = figure('Visible', 'on', 'Color', 'w');
    set(0, 'CurrentFigure', fg);

    for subj = 1:p.num_subjs
        for roi = 1:p.num_ROIs

            if toggles.save_voxel_plots
                figure_path = fullfile(figure_save_dir, ['S' num2str(subj)], roi_names{roi});
                if ~exist(figure_path, 'dir'), mkdir(figure_path); end
            end

            gt_params = simulated_data(subj, roi).params;
            est_params = all_pSFT{snr_idx_for_voxel_plots}(subj, roi).param_est;

            for vox = 1:num_voxels_to_plot

                %% pSFT curve (simulated and estimated)

                figure_name = ['Vox #' num2str(vox) ' pSFT'];

                % Simulated
                semilogx(p.sfs, simulated_data(subj, roi).SFT(:, vox), ...
                    'Color', plot_settings.colors.black, 'LineWidth', plot_settings.line_width);
                hold on;

                % Estimate
                semilogx(p.sfs, all_pSFT{snr_idx_for_voxel_plots}(subj, roi).est_SFT(:, vox), ...
                    'Color', plot_settings.colors.red, 'LineWidth', plot_settings.line_width);

                % Format figure
                xlabel('Spatial frequency (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                ylabel('Response (au)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                xlim([p.sf_min p.sf_max]); ylim([0 1]);
                xticks([p.sf_min 1 5 p.sf_max]); xticklabels([p.sf_min 1 5 p.sf_max]);
                set(gca, 'TickDir', 'out', 'TickLength', [plot_settings.tick_length plot_settings.tick_length], ...
                    'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                box off;

                title_str = sprintf('Simulated: \\mu=%.2f, \\sigma=%.2f | Estimated: \\mu=%.2f, \\sigma=%.2f', ...
                    gt_params(1, vox), gt_params(2, vox), est_params(1, vox), est_params(2, vox));
                title(title_str, 'FontSize', 10);

                legend({'simulated', 'estimate'}, 'Location', 'northeast', 'Box', 'off');

                hold off;

                if toggles.save_voxel_plots
                    saveas(gcf, [figure_path '/' figure_name '.pdf']); clf;
                    if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_path]); end
                end

                %% BOLD time series (simulated and estimated)

                figure_name = ['Vox #' num2str(vox) ' BOLD Time Series'];

                % Limit to specified number of runs
                num_TRs_to_plot = num_TRs_per_run * num_runs_to_plot;
                time_axis = (1:num_TRs_to_plot) * TR;

                % Simulated (noisy)
                plot(time_axis, all_pSFT{snr_idx_for_voxel_plots}(subj, roi).measured_BOLD(1:num_TRs_to_plot, vox), ...
                    'Color', plot_settings.colors.black, 'LineWidth', plot_settings.line_width);
                hold on;

                % Estimated
                plot(time_axis, all_pSFT{snr_idx_for_voxel_plots}(subj, roi).est_BOLD(1:num_TRs_to_plot, vox), ...
                    'Color', plot_settings.colors.green, 'LineWidth', plot_settings.line_width);

                % Format figure
                xlabel('Time (s)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                ylabel('BOLD (% change)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                xlim([0 num_TRs_to_plot * TR]);
                set(gca, 'TickDir', 'out', 'TickLength', [plot_settings.tick_length plot_settings.tick_length], ...
                    'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                box off;

                title(['R^2 = ' num2str(round(all_pSFT{snr_idx_for_voxel_plots}(subj, roi).r2(vox), 3))]);

                legend({'simulated', 'estimate'}, 'Location', 'best', 'Box', 'off');

                hold off;

                if toggles.save_voxel_plots
                    saveas(gcf, [figure_path '/' figure_name '.pdf']); clf;
                    if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_path]); end
                end

                %% BOLD time series at each noise level

                figure_name = ['Vox #' num2str(vox) ' BOLD by SNR Level'];

                % Compute common y-axis limits across all noise levels
                all_bold_values = simulated_data(subj, roi).BOLD(1:num_TRs_to_plot, vox);
                for snr_idx = 1:num_SNR_levels
                    all_bold_values = [all_bold_values; simulated_data(subj, roi).measured_BOLD{snr_idx}(1:num_TRs_to_plot, vox)];
                end
                bold_ylim = [min(all_bold_values) max(all_bold_values)];

                for snr_idx = 1:num_SNR_levels
                    subplot(num_SNR_levels, 1, snr_idx);

                    % Clean signal
                    plot(time_axis, simulated_data(subj, roi).BOLD(1:num_TRs_to_plot, vox), ...
                        'Color', plot_settings.colors.black, 'LineWidth', plot_settings.line_width);
                    hold on;

                    % Noisy signal at this SNR level
                    plot(time_axis, simulated_data(subj, roi).measured_BOLD{snr_idx}(1:num_TRs_to_plot, vox), ...
                        'Color', snr_colors{snr_idx}, 'LineWidth', plot_settings.line_width);

                    hold off;

                    % Format subplot
                    if snr_idx == num_SNR_levels
                        xlabel('Time (s)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                    end
                    ylabel('BOLD', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                    xlim([0 num_TRs_to_plot * TR]);
                    ylim(bold_ylim);
                    title(SNR_level_names{snr_idx}, 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                    box off;

                    if snr_idx == 1
                        legend({'clean', 'noisy'}, 'Location', 'northeast', 'Box', 'off');
                    end
                end

                if toggles.save_voxel_plots
                    saveas(gcf, [figure_path '/' figure_name '.pdf']); clf;
                    if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_path]); end
                end

            end

        end
    end

    close(fg);

end

%% Summary plots

if toggles.make_summary_plots

    fg = figure('Visible', 'on', 'Color', 'w');
    set(0, 'CurrentFigure', fg);

    %% Effect of noise on parameter recovery (scatter plots for each SNR level)

    figure_name = 'Simulated vs Estimated by SNR Level';

    % Compute common axis limits across all SNR levels for each parameter
    all_mu_values = [];
    all_bandwidth_values = [];
    for snr_idx = 1:num_SNR_levels
        all_mu_values = [all_mu_values, pooled_data{snr_idx}.all_gt_mu, pooled_data{snr_idx}.all_est_mu];
        all_bandwidth_values = [all_bandwidth_values, pooled_data{snr_idx}.all_gt_bandwidth_oct, pooled_data{snr_idx}.all_est_bandwidth_oct];
    end
    max_mu_rounded = ceil(max(all_mu_values));
    max_bandwidth_oct_rounded = ceil(max(all_bandwidth_values));

    % Compute tick values for common axes
    tick_step_mu = max(1, round(max_mu_rounded/4));
    tick_values_mu = 0:tick_step_mu:max_mu_rounded;
    tick_labels_mu = arrayfun(@(x) num2str(x), tick_values_mu, 'UniformOutput', false);

    tick_step_bandwidth_oct = max(1, round(max_bandwidth_oct_rounded/4));
    tick_values_bandwidth_oct = 0:tick_step_bandwidth_oct:max_bandwidth_oct_rounded;
    tick_labels_bandwidth_oct = arrayfun(@(x) num2str(x), tick_values_bandwidth_oct, 'UniformOutput', false);

    for snr_idx = 1:num_SNR_levels

        all_gt_mu = pooled_data{snr_idx}.all_gt_mu;
        all_est_mu = pooled_data{snr_idx}.all_est_mu;
        all_gt_bandwidth_oct = pooled_data{snr_idx}.all_gt_bandwidth_oct;
        all_est_bandwidth_oct = pooled_data{snr_idx}.all_est_bandwidth_oct;
        num_total_voxels = length(all_gt_mu);

        % Subplot row 1: pSFT peak (mu)
        subplot(2, num_SNR_levels, snr_idx);
        scatter(all_gt_mu, all_est_mu, 30, snr_colors{snr_idx}, 'filled', 'MarkerFaceAlpha', 0.7);
        hold on;
        plot([0 max_mu_rounded], [0 max_mu_rounded], 'k--', 'LineWidth', 1);
        hold off;

        % Format figure
        xlabel('\mu_{sim} (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
        if snr_idx == 1
            ylabel('\mu_{est} (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
        end
        xticks(tick_values_mu); yticks(tick_values_mu);
        xticklabels(tick_labels_mu); yticklabels(tick_labels_mu);
        ylim([0 max_mu_rounded]); xlim([0 max_mu_rounded]);
        title({SNR_level_names{snr_idx}, ['r = ' num2str(round(pooled_data{snr_idx}.corr_mu, 3))]}, ...
            'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
        set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
        axis square; box off;

        % Subplot row 2: bandwidth in octaves
        subplot(2, num_SNR_levels, num_SNR_levels + snr_idx);
        scatter(all_gt_bandwidth_oct, all_est_bandwidth_oct, 30, snr_colors{snr_idx}, 'filled', 'MarkerFaceAlpha', 0.7);
        hold on;
        plot([0 max_bandwidth_oct_rounded], [0 max_bandwidth_oct_rounded], 'k--', 'LineWidth', 1);
        hold off;

        % Format figure
        xlabel('\sigma_{sim} (oct)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
        if snr_idx == 1
            ylabel('\sigma_{est} (oct)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
        end
        xticks(tick_values_bandwidth_oct); yticks(tick_values_bandwidth_oct);
        xticklabels(tick_labels_bandwidth_oct); yticklabels(tick_labels_bandwidth_oct);
        ylim([0 max_bandwidth_oct_rounded]); xlim([0 max_bandwidth_oct_rounded]);
        title(['r = ' num2str(round(pooled_data{snr_idx}.corr_bandwidth_oct, 3))], ...
            'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
        set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
        axis square; box off;

    end

    if toggles.save_summary_plots
        saveas(gcf, fullfile(figure_save_dir, [figure_name '.pdf']));
        if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_save_dir]); end
    end
    clf;

    %% Effect of noise on RMSE (line plot with error bars)

    figure_name = 'Effect of Noise on Parameter Recovery';

    % Extract RMSE and correlation values with standard errors for each SNR level
    rmse_mu_by_snr = nan(1, num_SNR_levels);
    rmse_bandwidth_by_snr = nan(1, num_SNR_levels);
    corr_mu_by_snr = nan(1, num_SNR_levels);
    corr_bandwidth_by_snr = nan(1, num_SNR_levels);

    se_rmse_mu_by_snr = nan(1, num_SNR_levels);
    se_rmse_bandwidth_by_snr = nan(1, num_SNR_levels);
    se_corr_mu_by_snr = nan(1, num_SNR_levels);
    se_corr_bandwidth_by_snr = nan(1, num_SNR_levels);

    for snr_idx = 1:num_SNR_levels
        rmse_mu_by_snr(snr_idx) = pooled_data{snr_idx}.rmse_mu;
        rmse_bandwidth_by_snr(snr_idx) = pooled_data{snr_idx}.rmse_bandwidth_oct;
        corr_mu_by_snr(snr_idx) = pooled_data{snr_idx}.corr_mu;
        corr_bandwidth_by_snr(snr_idx) = pooled_data{snr_idx}.corr_bandwidth_oct;

        se_rmse_mu_by_snr(snr_idx) = pooled_data{snr_idx}.se_rmse_mu;
        se_rmse_bandwidth_by_snr(snr_idx) = pooled_data{snr_idx}.se_rmse_bandwidth_oct;
        se_corr_mu_by_snr(snr_idx) = pooled_data{snr_idx}.se_corr_mu;
        se_corr_bandwidth_by_snr(snr_idx) = pooled_data{snr_idx}.se_corr_bandwidth_oct;
    end

    % Subplot 1: RMSE vs SNR with error bars
    subplot(1, 2, 1);
    errorbar(SNR_levels_dB, rmse_mu_by_snr, se_rmse_mu_by_snr, '-o', 'Color', plot_settings.colors.black, ...
        'LineWidth', plot_settings.line_width, 'MarkerFaceColor', plot_settings.colors.black, 'MarkerSize', 8, 'CapSize', 0);
    hold on;
    errorbar(SNR_levels_dB, rmse_bandwidth_by_snr, se_rmse_bandwidth_by_snr, '-s', 'Color', plot_settings.colors.blue, ...
        'LineWidth', plot_settings.line_width, 'MarkerFaceColor', plot_settings.colors.blue, 'MarkerSize', 8, 'CapSize', 0);
    hold off;

    xlabel('SNR (dB)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylabel('RMSE', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    title('Effect of Noise on RMSE', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    legend({'\mu (cpd)', '\sigma (octaves)'}, 'Location', 'northeast', 'Box', 'off');
    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    axis square; box off;

    % Subplot 2: Correlation vs SNR with error bars
    subplot(1, 2, 2);
    errorbar(SNR_levels_dB, corr_mu_by_snr, se_corr_mu_by_snr, '-o', 'Color', plot_settings.colors.black, ...
        'LineWidth', plot_settings.line_width, 'MarkerFaceColor', plot_settings.colors.black, 'MarkerSize', 8, 'CapSize', 0);
    hold on;
    errorbar(SNR_levels_dB, corr_bandwidth_by_snr, se_corr_bandwidth_by_snr, '-s', 'Color', plot_settings.colors.blue, ...
        'LineWidth', plot_settings.line_width, 'MarkerFaceColor', plot_settings.colors.blue, 'MarkerSize', 8, 'CapSize', 0);
    hold off;

    xlabel('SNR (dB)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylabel('Correlation (r)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    title('Effect of Noise on Correlation', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    legend({'\mu (cpd)', '\sigma (octaves)'}, 'Location', 'southwest', 'Box', 'off');
    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    ylim([0 1]);
    axis square; box off;

    if toggles.save_summary_plots
        saveas(gcf, fullfile(figure_save_dir, [figure_name '.pdf']));
        if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_save_dir]); end
    end

    close(fg);

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
