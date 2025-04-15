% gridSearch - Perform a grid search to find the best mu and sigma for the pSFT model
%   
%   Syntax
%       [best_params, best_sse] = gridSearch(init_params, fixed_params, which_search, param_bounds)
%
%   Input Arguments
%       init_params – initial pSFT parameters
%       fixed_params – fixed parameters
%       which_search – which search to perform ('coarse' or 'fine')
%       param_bounds – bounds for the pSFT parameters
%
%   Output Arguments
%       best_params – best mu and sigma values
%       best_sse – best SSE value

function [best_params, best_sse] = gridSearch(init_params, fixed_params, which_search, param_bounds)

%% Setup mu and sigma values to search through

switch which_search
    
    case 'coarse'
        
        search_steps = 10;

        mu_grid    = logspace(log10(param_bounds(2,1)), log10(param_bounds(1,1)), search_steps); 
        sigma_grid = linspace(param_bounds(2,2), param_bounds(1,2), search_steps); 
        
    case 'fine'
        
        search_steps = 100;
        
        if init_params(1)*1.5 > param_bounds(1,1)
            mu_grid_ub = param_bounds(1,1);
        else
            mu_grid_ub = init_params(1)*1.5;
        end

        if init_params(1)*0.5 < param_bounds(2,1)
            mu_grid_lb = param_bounds(2,1);
        else
            mu_grid_lb = init_params(1)*0.5;
        end
        
        if init_params(2)*1.5 > param_bounds(1,2)
            sigma_grid_ub = param_bounds(1,2);
        else
            sigma_grid_ub = init_params(2)*1.5;
        end
           
        if init_params(2)*0.5 < param_bounds(2,2)
            sigma_grid_lb = param_bounds(2,2);
        else
            sigma_grid_lb = init_params(2)*0.5;
        end
        
        mu_grid    = logspace(log10(mu_grid_lb), log10(mu_grid_ub), search_steps); 
        sigma_grid = linspace(sigma_grid_lb, sigma_grid_ub, search_steps); 
        
end

%% Calculate sse for every combination of mu and sigma

sse_grid = nan(search_steps, search_steps);

for i_mu = 1:search_steps
    for i_sig = 1:search_steps
        
        pSFT_params = [mu_grid(i_mu) sigma_grid(i_sig) init_params(3:end)];
        sse_grid(i_mu, i_sig) = calcFit(pSFT_params, fixed_params);

    end
end

%% Find the best mu and sigma from the grid search

best_sse = min(sse_grid(:));

[row_indices, col_indices] = find(sse_grid == best_sse);

min_mu_values = mean(mu_grid(row_indices));
min_sigma_values = mean(sigma_grid(col_indices));

best_params = [min_mu_values, min_sigma_values, init_params(3:end)];

%% Inspect grid search

%{

[x,y] = meshgrid(sigma_grid, mu_grid);

figure('Color','w', 'Name','Grid Search')
contourf(x,y,sse_grid);

xlabel('sigma'); ylabel('mu'); zlabel('SSE');

% Change color maps
current_colormap = colormap('hot');  % Get the 'hot' colormap
reversed_colormap = flipud(current_colormap);  % Reverse the colormap
colormap(reversed_colormap);  % Set the reversed colormap as the new colormap

colorbar;

% Add cross on best combo
hold on;
plot(best_params(2), best_params(1),'g+', 'MarkerSize', 10, 'LineWidth', 2);
box off;

close all;

%}

end