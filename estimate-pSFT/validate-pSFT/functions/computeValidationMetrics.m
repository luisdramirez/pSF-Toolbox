% computeValidationMetrics - Compute validation metrics for pSFT estimation
%
% Computes RMSE, correlation, and R^2 metrics for each SNR level, both
% per-subject/ROI and pooled across all voxels.
%
% Inputs:
%   simulated_data   - Struct array with ground truth params and SFT curves
%   all_pSFT         - Cell array of estimated pSFT results per SNR level
%   SNR_levels_dB    - Array of SNR levels in dB
%   SNR_level_names  - Cell array of SNR level names for display
%   roi_names        - Cell array of ROI names
%   p                - Parameter struct with num_subjs, num_ROIs, sfs
%   disp_on          - Boolean for displaying progress
%
% Outputs:
%   validation_metrics - Cell array of per-subject/ROI metrics per SNR level
%   pooled_data        - Cell array of pooled data and metrics per SNR level
%   simulated_data     - Updated with bandwidth_oct field
%   all_pSFT           - Updated with est_bandwidth_oct field

function [validation_metrics, pooled_data, simulated_data, all_pSFT] = computeValidationMetrics(simulated_data, all_pSFT, SNR_levels_dB, SNR_level_names, roi_names, p, disp_on)

if disp_on, disp(' '); disp('=== Validation Metrics ==='); end

num_SNR_levels = length(SNR_levels_dB);

% Initialize validation metrics for each SNR level
validation_metrics = cell(1, num_SNR_levels);

% Initialize pooled data for each SNR level
pooled_data = cell(1, num_SNR_levels);

for snr_idx = 1:num_SNR_levels

    if disp_on
        disp(' ');
        disp(['--- ' SNR_level_names{snr_idx} ' ---']);
    end

    validation_metrics{snr_idx} = struct();

    % Collect all data across subjects and ROIs for pooled analysis
    all_gt_mu = [];
    all_est_mu = [];
    all_gt_sigma = [];
    all_est_sigma = [];
    all_gt_bandwidth_oct = [];
    all_est_bandwidth_oct = [];

    for subj = 1:p.num_subjs
        for roi = 1:p.num_ROIs

            gt_params = simulated_data(subj, roi).params;
            est_params = all_pSFT{snr_idx}(subj, roi).param_est;

            % Extract mu and sigma
            gt_mu = gt_params(1, :);
            gt_sigma = gt_params(2, :);
            est_mu = est_params(1, :);
            est_sigma = est_params(2, :);

            % Compute bandwidth in octaves
            num_voxels = size(simulated_data(subj, roi).SFT, 2);
            gt_bandwidth_oct = nan(1, num_voxels);
            est_bandwidth_oct = nan(1, num_voxels);

            for vox = 1:num_voxels
                [gt_bandwidth_oct(vox), ~] = cpd2oct(simulated_data(subj, roi).SFT(:, vox), p.sfs);
                [est_bandwidth_oct(vox), ~] = cpd2oct(all_pSFT{snr_idx}(subj, roi).est_SFT(:, vox), p.sfs);
            end

            % Store bandwidth in simulated_data for later use in plotting
            if snr_idx == 1
                simulated_data(subj, roi).bandwidth_oct = gt_bandwidth_oct;
            end
            all_pSFT{snr_idx}(subj, roi).est_bandwidth_oct = est_bandwidth_oct;

            % Accumulate for pooled analysis
            all_gt_mu = [all_gt_mu, gt_mu];
            all_est_mu = [all_est_mu, est_mu];
            all_gt_sigma = [all_gt_sigma, gt_sigma];
            all_est_sigma = [all_est_sigma, est_sigma];
            all_gt_bandwidth_oct = [all_gt_bandwidth_oct, gt_bandwidth_oct];
            all_est_bandwidth_oct = [all_est_bandwidth_oct, est_bandwidth_oct];

            % Compute metrics
            errors_mu = gt_mu - est_mu;
            errors_bandwidth_oct = gt_bandwidth_oct - est_bandwidth_oct;

            rmse_mu = sqrt(mean(errors_mu.^2));
            rmse_bandwidth_oct = sqrt(mean(errors_bandwidth_oct.^2));

            corr_mu = corr(gt_mu', est_mu');
            corr_bandwidth_oct = corr(gt_bandwidth_oct', est_bandwidth_oct');

            mean_r2 = mean(all_pSFT{snr_idx}(subj, roi).r2);

            % Store metrics
            validation_metrics{snr_idx}(subj, roi).rmse_mu = rmse_mu;
            validation_metrics{snr_idx}(subj, roi).rmse_bandwidth_oct = rmse_bandwidth_oct;
            validation_metrics{snr_idx}(subj, roi).corr_mu = corr_mu;
            validation_metrics{snr_idx}(subj, roi).corr_bandwidth_oct = corr_bandwidth_oct;
            validation_metrics{snr_idx}(subj, roi).mean_r2 = mean_r2;

            % Display metrics
            if disp_on
                disp(['S' num2str(subj) ', ' roi_names{roi} ':']);
                disp(['  RMSE (mu):           ' num2str(round(rmse_mu, 4))]);
                disp(['  RMSE (bandwidth oct): ' num2str(round(rmse_bandwidth_oct, 4))]);
                disp(['  Corr (mu):           ' num2str(round(corr_mu, 4))]);
                disp(['  Corr (bandwidth oct): ' num2str(round(corr_bandwidth_oct, 4))]);
                disp(['  Mean R^2:            ' num2str(round(mean_r2, 4))]);
                disp(' ');
            end

        end
    end

    % Store pooled data
    pooled_data{snr_idx}.all_gt_mu = all_gt_mu;
    pooled_data{snr_idx}.all_est_mu = all_est_mu;
    pooled_data{snr_idx}.all_gt_sigma = all_gt_sigma;
    pooled_data{snr_idx}.all_est_sigma = all_est_sigma;
    pooled_data{snr_idx}.all_gt_bandwidth_oct = all_gt_bandwidth_oct;
    pooled_data{snr_idx}.all_est_bandwidth_oct = all_est_bandwidth_oct;
    pooled_data{snr_idx}.SNR_dB = SNR_levels_dB(snr_idx);

    % Compute pooled RMSE and correlation
    pooled_data{snr_idx}.rmse_mu = sqrt(mean((all_gt_mu - all_est_mu).^2));
    pooled_data{snr_idx}.rmse_bandwidth_oct = sqrt(mean((all_gt_bandwidth_oct - all_est_bandwidth_oct).^2));
    pooled_data{snr_idx}.corr_mu = corr(all_gt_mu', all_est_mu');
    pooled_data{snr_idx}.corr_bandwidth_oct = corr(all_gt_bandwidth_oct', all_est_bandwidth_oct');

    % Bootstrap confidence intervals for RMSE and correlation
    num_bootstrap = 1000;
    num_voxels_total = length(all_gt_mu);

    boot_rmse_mu = nan(1, num_bootstrap);
    boot_rmse_bandwidth = nan(1, num_bootstrap);
    boot_corr_mu = nan(1, num_bootstrap);
    boot_corr_bandwidth = nan(1, num_bootstrap);

    for boot = 1:num_bootstrap
        % Resample voxel indices with replacement
        boot_idx = randi(num_voxels_total, 1, num_voxels_total);

        % Compute bootstrapped metrics
        boot_rmse_mu(boot) = sqrt(mean((all_gt_mu(boot_idx) - all_est_mu(boot_idx)).^2));
        boot_rmse_bandwidth(boot) = sqrt(mean((all_gt_bandwidth_oct(boot_idx) - all_est_bandwidth_oct(boot_idx)).^2));
        boot_corr_mu(boot) = corr(all_gt_mu(boot_idx)', all_est_mu(boot_idx)');
        boot_corr_bandwidth(boot) = corr(all_gt_bandwidth_oct(boot_idx)', all_est_bandwidth_oct(boot_idx)');
    end

    % Store bootstrap standard errors (for error bars)
    pooled_data{snr_idx}.se_rmse_mu = std(boot_rmse_mu);
    pooled_data{snr_idx}.se_rmse_bandwidth_oct = std(boot_rmse_bandwidth);
    pooled_data{snr_idx}.se_corr_mu = std(boot_corr_mu);
    pooled_data{snr_idx}.se_corr_bandwidth_oct = std(boot_corr_bandwidth);

    % Store 95% confidence intervals
    pooled_data{snr_idx}.ci_rmse_mu = prctile(boot_rmse_mu, [2.5 97.5]);
    pooled_data{snr_idx}.ci_rmse_bandwidth_oct = prctile(boot_rmse_bandwidth, [2.5 97.5]);
    pooled_data{snr_idx}.ci_corr_mu = prctile(boot_corr_mu, [2.5 97.5]);
    pooled_data{snr_idx}.ci_corr_bandwidth_oct = prctile(boot_corr_bandwidth, [2.5 97.5]);

    if disp_on
        disp(['Pooled metrics for ' SNR_level_names{snr_idx} ':']);
        disp(['  RMSE (mu):           ' num2str(round(pooled_data{snr_idx}.rmse_mu, 4)) ' ±' num2str(round(pooled_data{snr_idx}.se_rmse_mu, 4))]);
        disp(['  RMSE (bandwidth oct): ' num2str(round(pooled_data{snr_idx}.rmse_bandwidth_oct, 4)) ' ±' num2str(round(pooled_data{snr_idx}.se_rmse_bandwidth_oct, 4))]);
        disp(['  Corr (mu):           ' num2str(round(pooled_data{snr_idx}.corr_mu, 4)) ' ±' num2str(round(pooled_data{snr_idx}.se_corr_mu, 4))]);
        disp(['  Corr (bandwidth oct): ' num2str(round(pooled_data{snr_idx}.corr_bandwidth_oct, 4)) ' ±' num2str(round(pooled_data{snr_idx}.se_corr_bandwidth_oct, 4))]);
    end

end

end

