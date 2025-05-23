% example_pipeline

%% Prepare workspace

clear all; 
close all; 
clc;

project_dir = pwd;
save_dir = 'estimates';
data_dir = 'data';
figure_dir = 'figures';

addpath(fullfile(project_dir, 'functions'));
addpath(fullfile(project_dir, data_dir));    

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

save_pSFT = true;
make_voxel_plots = true;
save_voxel_plots = true;

%% Check for required toolboxes

checkRequiredToolboxes(toggles);

%% Parallelization setup for parfor loop

p.num_cores = 4;
p.num_chunks = p.num_cores-1;

if toggles.parallelization
    maxNumCompThreads(p.num_cores);
    if p.num_cores > 1
        parpool('local', p.num_chunks)
    end
end

%% Spatial frequency parameters
% Hardcoded SF parameters used for generating pSFT curves.
% Ideally, the range should match the range of SFs presented.

p.sf_min = 0.5;
p.sf_max = 12;

p.sf_count = 100;
p.sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), p.sf_count);

%% Initalize pSFT model parameters

p.initial_params = [1 1 1 0]; % [mu, sigma, beta, beta_0]

p.pSFT_bounds(1,:) = [6, 4, 25, 10]; % upper bounds
p.pSFT_bounds(2,:) = [0.009, 0.1, -25, -10]; % lower bounds

p.fmincon_options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000, 'display','off'); 

%% Load data

load(fullfile(data_dir, 'sample_data.mat'));

num_subjs = length(sample_data);
num_ROIs = length(sample_data(1).measured_BOLD);

roi_names = {'V1', 'V2', 'V3'};

%% Define hemodynamic response function (HRF)

HRF = defineHRF();

%% Initialize pSFT struct

struct_size = cell(num_subjs, num_ROIs);    

all_pSFT = struct('vox_indices', struct_size, ...
    'param_est', struct_size, ...    
    'est_SFT', struct_size, ...
    'est_R', struct_size, ...
    'est_BOLD', struct_size, ...
    'r2', struct_size, ...
    'sse', struct_size, ...
    'start_values', struct_size, ...
    'start_sse', struct_size, ...
    'exitflag', struct_size, ...
    'measured_BOLD', struct_size);

%% Estimate pSFT for each subject and ROI

total_elapsed_time = 0;

for subj = 1:num_subjs

    if toggles.disp_on, disp(['+++ S' num2str(subj) ' +++']); end

    I = sample_data(subj).I;
    
    for roi = 1:num_ROIs

        if toggles.disp_on, disp(['- ' roi_names{roi} ' -']); end

        tic;    

        measured_BOLD = sample_data(subj).measured_BOLD{roi};
        
        % Enter estimatePSFT
        pSFT = estimatePSFT(I, measured_BOLD, HRF, p, toggles);
        all_pSFT(subj,roi) = pSFT;

        elapsed_time = round(toc/60,1);
        if toggles.disp_on, disp(['Elapsed time: ~' num2str(elapsed_time) ' minute(s).']); disp(' '); end
        total_elapsed_time = total_elapsed_time + elapsed_time;

    end
end

if toggles.disp_on, disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']); end
if toggles.parallelization, delete(gcp); end

if save_pSFT
    curr_time = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    filename = ['all_pSFT_n' num2str(num_subjs) '_' char(curr_time) '.mat'];
    save([save_dir '/' filename], 'all_pSFT');
    if toggles.disp_on, disp(['Saved ' filename ' in /' save_dir]); end
end

%% Plot voxel fits

if make_voxel_plots

    num_voxels = 1;

    fg = figure('Visible','on','Color','w');
    set(0,'CurrentFigure',fg);

    for subj = 1:num_subjs
        for roi = 1:num_ROIs

            if save_voxel_plots
                figure_path = [figure_dir '/S' num2str(subj) '/V' num2str(roi)];
                if exist(figure_path, 'dir') == 0, mkdir(figure_path); end
            end

            for vox = 1:num_voxels

                %% pSFT curve

                figure_name = ['Vox #' num2str(vox) ' pSFT'];
                
                semilogx(p.sfs, all_pSFT(subj,roi).est_SFT(:, vox), 'k');
            
                % Format figure
                xlabel('log[SF] (cpd)'); ylabel('R');
                xlim([p.sf_min p.sf_max]); ylim([0 1]);
                xticks([p.sf_min 1 5 p.sf_max]); xticklabels([p.sf_min 1 5 p.sf_max]); 
                set(gca,'TickDir','out'); box off;
                title(['\mu = ' num2str(round(all_pSFT(subj,roi).param_est(1,vox),2)) ' | \sigma = ' num2str(round(all_pSFT(subj,roi).param_est(2,vox),2))])
                
                if save_voxel_plots 
                    saveas(gcf,[figure_path '/' figure_name '.png']); clf; 
                    if toggles.disp_on, disp(['Saved ' figure_name '.png in /' figure_path]); end
                end

                %% Voxel time series

                figure_name = ['Vox #' num2str(vox) ' BOLD Time Series'];
                
                plot(all_pSFT(subj,roi).measured_BOLD(:,vox),'k');
                hold on
                plot(all_pSFT(subj,roi).est_BOLD(:,vox),'r');
                
                % Format figure
                xlabel('Time (s)'); ylabel('BOLD (% change)');
                set(gca,'TickDir','out'); box off;
                title(['R^2 = ' num2str(round(all_pSFT(subj,roi).r2(vox),2))])
                legend({'measured','estimate'})
    
                if save_voxel_plots 
                    saveas(gcf,[figure_path '/' figure_name '.png']); clf; 
                    if toggles.disp_on, disp(['Saved ' figure_name '.png in /' figure_path]); end
                end
   
            end
            
        end
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