% validate_pSFT
%
% This script validates the pSFT estimation pipeline by simulating stimulus and BOLD time series
% from pre-defined pSFT parameters.
%
% It then estimates the pSFT parameters and compares the estimates to the ground truth.
%
% Results are saved under `validate-pSFT/estimates/validation_pSFT_n<num_subjs>_yyyy-mm-dd_hh-mm-ss/` and `validate-pSFT/figures/validation_pSFT_n<num_subjs>_yyyy-mm-dd_hh-mm-ss/`.

%% Prepare workspace

clear all;
close all;
clc;

curr_time = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));

project_dir = pwd;
save_dir = 'estimates';
figure_dir = 'figures';

addpath('../functions');

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

save_estimates = true;
make_voxel_plots = true;
save_voxel_plots = true;

%% Check for required toolboxes

checkRequiredToolboxes(toggles);

%% Parallelization setup for parfor loop

p.num_cores = 4;
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
% For reviewer demonstration: use multiple subjects, ROIs, and voxels
% to demonstrate validation across different conditions

p.num_subjs = 2;
p.num_ROIs = 3;
p.num_voxels_per_ROI = 20;

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

% BOLD noise parameters
min_BOLD_noise = 0.01;
max_BOLD_noise = 0.1;

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

    % Get subject-specific stimulus time series
    I = I_all{subj};

    for roi = 1:p.num_ROIs

        gt_params = simulated_data(subj, roi).params;
        HRF = HRFs{subj, roi};

        % Store I for this subject-ROI combination
        simulated_data(subj, roi).I = I;

        % Pre-allocate
        R = nan(num_TRs, p.num_voxels_per_ROI);
        clean_BOLD = nan(num_TRs, p.num_voxels_per_ROI);
        noisy_BOLD = nan(num_TRs, p.num_voxels_per_ROI);

        for vox = 1:p.num_voxels_per_ROI

            % Generate neural response from ground truth pSFT
            vox_R = pSFT_model([gt_params(1, vox), gt_params(2, vox)], I);
            R(:, vox) = vox_R;

            % Generate BOLD response using generateBOLD function
            vox_BOLD = generateBOLD(vox_R, HRF(:, vox), gt_params(3, vox), gt_params(4, vox));
            clean_BOLD(:, vox) = vox_BOLD;

            % Add noise to create measured BOLD
            noise_level = min_BOLD_noise + (max_BOLD_noise - min_BOLD_noise) * rand();
            noise = noise_level * randn(num_TRs, 1);
            noisy_BOLD(:, vox) = vox_BOLD + noise;

        end

        simulated_data(subj, roi).R = R;
        simulated_data(subj, roi).BOLD = clean_BOLD;
        simulated_data(subj, roi).measured_BOLD = noisy_BOLD;

    end

end

if toggles.disp_on, disp('Done!'); end

%% Initialize pSFT struct

struct_size = cell(p.num_subjs, p.num_ROIs);

all_pSFT = struct('vox_indices', struct_size, ...
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

%% Estimate pSFT

total_elapsed_time = 0;

for subj = 1:p.num_subjs

    if toggles.disp_on, disp(['+++ S' num2str(subj) ' +++']); end

    for roi = 1:p.num_ROIs

        if toggles.disp_on, disp(['- ' roi_names{roi} ' -']); end

        tic;

        I = simulated_data(subj, roi).I;
        measured_BOLD = simulated_data(subj, roi).measured_BOLD;
        HRF = HRFs{subj, roi};

        % Enter estimatePSFT
        pSFT = estimatePSFT(I, measured_BOLD, HRF, p, toggles);
        all_pSFT(subj, roi) = pSFT;

        elapsed_time = round(toc / 60, 1);
        if toggles.disp_on, disp(['Elapsed time: ~' num2str(elapsed_time) ' minute(s).']); disp(' '); end
        total_elapsed_time = total_elapsed_time + elapsed_time;

    end
end

if toggles.disp_on, disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']); end
if toggles.parallelization, delete(gcp); end

%% Save estimates

estimates_filename = ['validation_pSFT_n' num2str(p.num_subjs) '_' char(curr_time)];

if save_estimates
    save([save_dir '/' estimates_filename '.mat'], 'all_pSFT', 'simulated_data', 'p');
    if toggles.disp_on, disp(['Saved ' estimates_filename ' in /' save_dir]); end
end

%% Compute validation metrics

if toggles.disp_on, disp(' '); disp('=== Validation Metrics ==='); end

validation_metrics = struct();

for subj = 1:p.num_subjs
    for roi = 1:p.num_ROIs

        gt_params = simulated_data(subj, roi).params;
        est_params = all_pSFT(subj, roi).param_est;

        % Extract mu and sigma (the key pSFT parameters)
        gt_mu = gt_params(1, :);
        gt_sigma = gt_params(2, :);
        est_mu = est_params(1, :);
        est_sigma = est_params(2, :);

        % Compute bandwidth in octaves for ground truth and estimated pSFT curves
        num_voxels = size(simulated_data(subj, roi).SFT, 2);
        gt_bandwidth_oct = nan(1, num_voxels);
        est_bandwidth_oct = nan(1, num_voxels);

        for vox = 1:num_voxels
            % Ground truth bandwidth in octaves
            [gt_bandwidth_oct(vox), ~] = cpd2oct(simulated_data(subj, roi).SFT(:, vox), p.sfs);

            % Estimated bandwidth in octaves
            [est_bandwidth_oct(vox), ~] = cpd2oct(all_pSFT(subj, roi).est_SFT(:, vox), p.sfs);
        end

        % Compute RMSE for mu
        rmse_mu = sqrt(mean((gt_mu - est_mu).^2));

        % Compute RMSE for sigma
        rmse_sigma = sqrt(mean((gt_sigma - est_sigma).^2));

        % Compute RMSE for bandwidth in octaves
        rmse_bandwidth_oct = sqrt(mean((gt_bandwidth_oct - est_bandwidth_oct).^2));

        % Compute correlation between ground truth and estimates
        corr_mu = corr(gt_mu', est_mu');
        corr_sigma = corr(gt_sigma', est_sigma');
        corr_bandwidth_oct = corr(gt_bandwidth_oct', est_bandwidth_oct');

        % Mean R^2 across voxels
        mean_r2 = mean(all_pSFT(subj, roi).r2);

        % Store metrics
        validation_metrics(subj, roi).rmse_mu = rmse_mu;
        validation_metrics(subj, roi).rmse_sigma = rmse_sigma;
        validation_metrics(subj, roi).rmse_bandwidth_oct = rmse_bandwidth_oct;
        validation_metrics(subj, roi).corr_mu = corr_mu;
        validation_metrics(subj, roi).corr_sigma = corr_sigma;
        validation_metrics(subj, roi).corr_bandwidth_oct = corr_bandwidth_oct;
        validation_metrics(subj, roi).mean_r2 = mean_r2;

        % Display metrics
        if toggles.disp_on
            disp(['S' num2str(subj) ', ' roi_names{roi} ':']);
            disp(['  RMSE (mu):           ' num2str(round(rmse_mu, 4))]);
            disp(['  RMSE (sigma):        ' num2str(round(rmse_sigma, 4))]);
            disp(['  RMSE (bandwidth oct): ' num2str(round(rmse_bandwidth_oct, 4))]);
            disp(['  Corr (mu):           ' num2str(round(corr_mu, 4))]);
            disp(['  Corr (sigma):        ' num2str(round(corr_sigma, 4))]);
            disp(['  Corr (bandwidth oct): ' num2str(round(corr_bandwidth_oct, 4))]);
            disp(['  Mean R^2:            ' num2str(round(mean_r2, 4))]);
            disp(' ');
        end

    end
end

%% Plot settings

plot_settings = plotSettings();

%% Plot voxel fits

figure_save_dir = fullfile(figure_dir, estimates_filename);
if ~exist(figure_save_dir, 'dir'), mkdir(figure_save_dir); end

if make_voxel_plots

    num_voxels_to_plot = 3;

    fg = figure('Visible', 'on', 'Color', 'w');
    set(0, 'CurrentFigure', fg);

    for subj = 1:p.num_subjs
        for roi = 1:p.num_ROIs

            if save_voxel_plots
                figure_path = fullfile(figure_save_dir, ['S' num2str(subj)], roi_names{roi});
                if ~exist(figure_path, 'dir'), mkdir(figure_path); end
            end

            gt_params = simulated_data(subj, roi).params;
            est_params = all_pSFT(subj, roi).param_est;

            for vox = 1:num_voxels_to_plot

                %% pSFT curve: Ground truth vs Estimate

                figure_name = ['Vox #' num2str(vox) ' pSFT'];

                % Ground truth
                semilogx(p.sfs, simulated_data(subj, roi).SFT(:, vox), ...
                    'Color', plot_settings.colors.black, 'LineWidth', plot_settings.line_width);
                hold on;

                % Estimate
                semilogx(p.sfs, all_pSFT(subj, roi).est_SFT(:, vox), ...
                    'Color', plot_settings.colors.red, 'LineWidth', plot_settings.line_width);

                % Format figure
                xlabel('Spatial Frequency (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                ylabel('Response (R)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                xlim([p.sf_min p.sf_max]); ylim([0 1]);
                xticks([p.sf_min 1 5 p.sf_max]); xticklabels([p.sf_min 1 5 p.sf_max]);
                set(gca, 'TickDir', 'out', 'TickLength', [plot_settings.tick_length plot_settings.tick_length], ...
                    'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                box off;

                title_str = sprintf('GT: \\mu=%.2f, \\sigma=%.2f | Est: \\mu=%.2f, \\sigma=%.2f', ...
                    gt_params(1, vox), gt_params(2, vox), est_params(1, vox), est_params(2, vox));
                title(title_str, 'FontSize', 10);

                legend({'ground truth', 'estimate'}, 'Location', 'northeast', 'Box', 'off');

                hold off;

                if save_voxel_plots
                    saveas(gcf, [figure_path '/' figure_name '.pdf']); clf;
                    if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_path]); end
                end

                %% BOLD time series: Measured vs Estimated vs Ground Truth

                figure_name = ['Vox #' num2str(vox) ' BOLD Time Series'];

                time_axis = (1:num_TRs) * TR;

                % Measured
                plot(time_axis, all_pSFT(subj, roi).measured_BOLD(:, vox), ...
                    'Color', plot_settings.colors.black, 'LineWidth', 0.5);
                hold on;

                % Estimated
                plot(time_axis, all_pSFT(subj, roi).est_BOLD(:, vox), ...
                    'Color', plot_settings.colors.red, 'LineWidth', plot_settings.line_width);

                % Format figure
                xlabel('Time (s)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                ylabel('BOLD (% change)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
                set(gca, 'TickDir', 'out', 'TickLength', [plot_settings.tick_length plot_settings.tick_length], ...
                    'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
                box off;

                title(['R^2 = ' num2str(round(all_pSFT(subj, roi).r2(vox), 3))]);

                legend({'simulated', 'estimate'}, 'Location', 'best', 'Box', 'off');

                hold off;

                if save_voxel_plots
                    saveas(gcf, [figure_path '/' figure_name '.pdf']); clf;
                    if toggles.disp_on, disp(['Saved ' figure_name '.pdf in /' figure_path]); end
                end

            end

        end
    end

    %% Summary scatter plots: Ground truth vs Estimates

    figure_name = 'Validation Summary';

    % Collect all data across subjects and ROIs
    all_gt_mu = [];
    all_est_mu = [];
    all_gt_sigma = [];
    all_est_sigma = [];
    all_gt_bandwidth_oct = [];
    all_est_bandwidth_oct = [];

    for subj = 1:p.num_subjs
        for roi = 1:p.num_ROIs
            all_gt_mu = [all_gt_mu, simulated_data(subj, roi).params(1, :)];
            all_est_mu = [all_est_mu, all_pSFT(subj, roi).param_est(1, :)];
            all_gt_sigma = [all_gt_sigma, simulated_data(subj, roi).params(2, :)];
            all_est_sigma = [all_est_sigma, all_pSFT(subj, roi).param_est(2, :)];

            % Calculate bandwidth in octaves for ground truth and estimated pSFT curves
            num_voxels = size(simulated_data(subj, roi).SFT, 2);
            gt_bandwidth_oct_roi = nan(1, num_voxels);
            est_bandwidth_oct_roi = nan(1, num_voxels);

            for vox = 1:num_voxels
                % Ground truth bandwidth in octaves
                [gt_bandwidth_oct_roi(vox), ~] = cpd2oct(simulated_data(subj, roi).SFT(:, vox), p.sfs);

                % Estimated bandwidth in octaves
                [est_bandwidth_oct_roi(vox), ~] = cpd2oct(all_pSFT(subj, roi).est_SFT(:, vox), p.sfs);
            end

            all_gt_bandwidth_oct = [all_gt_bandwidth_oct, gt_bandwidth_oct_roi];
            all_est_bandwidth_oct = [all_est_bandwidth_oct, est_bandwidth_oct_roi];
        end
    end

    % Calculate total number of voxels
    num_total_voxels = length(all_gt_mu);

    % Subplot 1: pSFT peak
    subplot(1, 3, 1);
    scatter(all_gt_mu, all_est_mu, 50, plot_settings.colors.green, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    max_mu = max([all_gt_mu, all_est_mu]);
    plot([0 max_mu], [0 max_mu], 'k--', 'LineWidth', 1);

    % Format figure
    xlabel('Ground Truth \mu (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylabel('Estimated \mu (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylim([0 max_mu]); xlim([0 max_mu]);
    title(['r = ' num2str(round(corr(all_gt_mu', all_est_mu'), 3)) ', n = ' num2str(num_total_voxels)]);
    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    axis square; box off;
    hold off;

    % Subplot 2: pSFT linear bandwidth
    subplot(1, 3, 2);
    scatter(all_gt_sigma, all_est_sigma, 50, plot_settings.colors.green, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    max_sigma = max([all_gt_sigma, all_est_sigma]);
    plot([0 max_sigma], [0 max_sigma], 'k--', 'LineWidth', 1);

    % Format figure
    xlabel('Ground Truth \sigma (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylabel('Estimated \sigma (cpd)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylim([0 max_sigma]); xlim([0 max_sigma]);
    title(['r = ' num2str(round(corr(all_gt_sigma', all_est_sigma'), 3)) ', n = ' num2str(num_total_voxels)]);
    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    axis square; box off;
    hold off;

    % Subplot 3: pSFT bandwidth in octaves
    subplot(1, 3, 3);
    scatter(all_gt_bandwidth_oct, all_est_bandwidth_oct, 50, plot_settings.colors.green, 'filled', 'MarkerFaceAlpha', 0.7);
    hold on;
    max_bandwidth_oct = max([all_gt_bandwidth_oct, all_est_bandwidth_oct]);
    plot([0 max_bandwidth_oct], [0 max_bandwidth_oct], 'k--', 'LineWidth', 1);

    % Format figure
    xlabel('Ground Truth Bandwidth (octaves)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylabel('Estimated Bandwidth (octaves)', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_label_font_size);
    ylim([0 max_bandwidth_oct]); xlim([0 max_bandwidth_oct]);
    title(['r = ' num2str(round(corr(all_gt_bandwidth_oct', all_est_bandwidth_oct'), 3)) ', n = ' num2str(num_total_voxels)]);
    set(gca, 'TickDir', 'out', 'FontName', plot_settings.font_type, 'FontSize', plot_settings.axes_tick_font_size);
    axis square; box off;
    hold off;

    if save_voxel_plots
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
