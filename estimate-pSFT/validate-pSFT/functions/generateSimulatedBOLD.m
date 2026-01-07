% generateSimulatedBOLD - Generate simulated BOLD time series with noise
%
% Generates clean BOLD responses from ground truth pSFT parameters and adds
% Gaussian noise at multiple SNR levels
%
% Inputs:
%   gt_params       - Ground truth parameters [mu; sigma; beta; beta_0] (4 x num_voxels)
%   I               - Stimulus time series (vector)
%   HRF             - HRF matrix (num_samples x num_voxels)
%   pSFT_model      - Function handle for pSFT model (e.g., @logGauss)
%   SNR_levels_dB   - Array of SNR levels in dB
%
% Outputs:
%   data            - Struct with fields: I, R, BOLD, measured_BOLD, SNR_dB

function data = generateSimulatedBOLD(gt_params, I, HRF, pSFT_model, SNR_levels_dB)

num_TRs = length(I);
num_voxels = size(gt_params, 2);
num_SNR_levels = length(SNR_levels_dB);

% Store stimulus
data.I = I;

% Pre-allocate
R = nan(num_TRs, num_voxels);
clean_BOLD = nan(num_TRs, num_voxels);

for vox = 1:num_voxels
    % Generate neural response from ground truth pSFT
    vox_R = pSFT_model([gt_params(1, vox), gt_params(2, vox)], I);
    R(:, vox) = vox_R;

    % Generate BOLD response
    vox_BOLD = generateBOLD(vox_R, HRF(:, vox), gt_params(3, vox), gt_params(4, vox));
    clean_BOLD(:, vox) = vox_BOLD;
end

data.R = R;
data.BOLD = clean_BOLD;

% Add noise at each SNR level
for snr_idx = 1:num_SNR_levels
    noisy_BOLD = nan(num_TRs, num_voxels);
    SNR_dB = SNR_levels_dB(snr_idx);

    for vox = 1:num_voxels
        signal_std = std(clean_BOLD(:, vox));
        noise_std = signal_std / (10^(SNR_dB / 20));
        noise = noise_std * randn(num_TRs, 1);
        noisy_BOLD(:, vox) = clean_BOLD(:, vox) + noise;
    end

    data.measured_BOLD{snr_idx} = noisy_BOLD;
    data.SNR_dB(snr_idx) = SNR_dB;
end

end
