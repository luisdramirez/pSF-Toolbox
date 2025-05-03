% fitVoxels
%   Fits pSFT parameters to each voxel
%
%   Input Arguments
%       chunk – structure of chunked data (see estimatePSF.m)
%       I – matrix of input spatial frequency time series [time x 1]
%       HIRF – hemodynamic impulse response function [time x 1]
%       p – structure of parameters (see estimatePSF.m)
%       toggles – structure of toggles (see estimatePSF.m)
%
%   Output Arguments
%       chunk – structure with fields:
%           -   start_values – initial pSFT parameters [4 x voxels]
%           -   start_sse – initial SSE values [1 x voxels]
%           -   param_est – estimated pSFT parameters [4 x voxels]
%           -   exitflag – fmincon exit flags [1 x voxels]
%           -   sse – SSE values [1 x voxels]
%           -   r2 – R^2 values [1 x voxels]
%           -   est_SFT – estimated SFT curves [sf_count x voxels]
%           -   est_R – estimated neural time series [time x voxels]
%           -   est_BOLD – estimated BOLD time series [time x voxels]

function chunk = fitVoxels(chunk, I, HIRF, p, toggles)

    num_vox = size(chunk.measured_BOLD,2);

    %% Pre-allocate arrays

    tmp_start_values = nan(size(p.initial_params,2), num_vox);
    tmp_start_sse = nan(1, num_vox); 

    tmp_param_est = nan(size(p.initial_params,2), num_vox);
    tmp_sse = nan(1, num_vox);
    tmp_exitflag = nan(1, num_vox);
    tmp_r2 = nan(1, num_vox);

    tmp_est_SFT = nan(p.sf_count, num_vox);
    tmp_est_R = nan(size(chunk.measured_BOLD,1), num_vox);
    tmp_est_BOLD = nan(size(chunk.measured_BOLD,1), num_vox);

    %% Loop through voxels

    for vox = 1:num_vox
        
        curr_fixed_params_vox = {I, chunk.measured_BOLD(:,vox), HIRF};
        
        %% Coarse-to-fine grid search for starting pSFT parameters
        % pSFT params: [mu, sigma, beta, beta_0]
        
        if toggles.coarse_grid_search
            pSFT_start_vals = gridSearch(p.initial_params, curr_fixed_params_vox, 'coarse', p.pSFT_bounds);
        else
            pSFT_start_vals = p.initial_params;
        end
        
        if toggles.fine_grid_search
            pSFT_start_vals = gridSearch(pSFT_start_vals, curr_fixed_params_vox, 'fine', p.pSFT_bounds);
        end
        
        tmp_start_values(:,vox) = pSFT_start_vals;
        tmp_start_sse(vox) = calcFit(pSFT_start_vals, curr_fixed_params_vox);
                
        %% fmincon to estimate pSFT

        pSFT_model = @(pSFT_free_params) calcFit(pSFT_free_params, curr_fixed_params_vox);
        
        [param_est, sse, exit_flag] = fmincon(pSFT_model, pSFT_start_vals, [],[],[],[], p.pSFT_bounds(2,:), p.pSFT_bounds(1,:), [], p.fmincon_options);
        
        tmp_param_est(:,vox) = param_est;
        tmp_sse(vox) = sse;
        tmp_exitflag(vox) = exit_flag;
        
        %% Evaluate pSFT estimates
        
        % Generate neural response
        R = logGauss(param_est, curr_fixed_params_vox{1});
        R(isnan(R)) = 0;
        tmp_est_R(:,vox) = R;
        
        % Generate BOLD response
        est_BOLD =  conv(R, curr_fixed_params_vox{3});
        est_BOLD = est_BOLD(1:length(R));
        est_BOLD = param_est(3) .* ((est_BOLD./mean(est_BOLD))-1) + param_est(4);
        tmp_est_BOLD(:,vox) = est_BOLD;
        
        % Calculate R^2
        tmp_r2(vox) = 1 - (sum((curr_fixed_params_vox{2} - est_BOLD).^2) / sum((curr_fixed_params_vox{2} - mean(curr_fixed_params_vox{2})).^2));
        
        % Generate SFT curve
        est_SFT = logGauss(param_est, p.sfs);
        tmp_est_SFT(:,vox) = est_SFT;
        
    end

    %% Store chunk

    chunk.start_values = tmp_start_values;
    chunk.start_sse = tmp_start_sse;
    chunk.param_est = tmp_param_est;
    chunk.exitflag = tmp_exitflag;
    chunk.sse = tmp_sse;
    chunk.r2 = tmp_r2;
    chunk.est_SFT = tmp_est_SFT;
    chunk.est_R = tmp_est_R;
    chunk.est_BOLD = tmp_est_BOLD;

end