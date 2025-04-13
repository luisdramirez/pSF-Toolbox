% example_pipeline

%% Prepare workspace

clear all; close all; clc;

project_dir = pwd;

addpath([project_dir '/functions']);
addpath([project_dir '/data']);    

if exist([project_dir '/estimates'], 'dir') == 0, mkdir([project_dir '/estimates']); end
if exist([project_dir '/figures'], 'dir') == 0, mkdir([project_dir '/figures']); end

addpath([project_dir '/estimates']);
addpath([project_dir '/figures']);

%% Toggles

toggles.parallelization = true;
toggles.coarse_grid_search = true;
toggles.fine_grid_search = true;
toggles.disp_on = true; % Display progress

save_pSF = true;
make_voxel_plots = true;
save_voxel_plots = true;

%% Parallelization setup for parfor loop

p.num_cores = 8;
p.num_chunks = p.num_cores-1;

if toggles.parallelization
    maxNumCompThreads(p.num_cores);
    if p.num_cores > 1
        parpool('local', p.num_chunks)
    end
end

%% Spatial frequency parameters
% These are hardcoded SF parameters used for generating pSFT curves.
% Ideally, the range should match the range of SFs presented.

p.sf_min = 0.1;
p.sf_max = 12;

p.sf_count = 100;
p.sfs = 10.^linspace(log10(p.sf_min), log10(p.sf_max), p.sf_count);

%% Initalize pSFT model parameters
% [mu, sigma, beta, beta_0]

p.init_params = [1 1 1 0]; 

p.pSFT_bounds(1,:) = [6, 4, 25, 10];
p.pSFT_bounds(2,:) = [0.009, 0.1, -25, -10]; 

%% Load data

load('data/sample_data.mat');

num_subjs = length(sample_data.measured_BOLD);
num_ROIs = size(sample_data.measured_BOLD{1}, 3);

%% Initialize pSF struct

struct_size = cell(num_subjs, num_ROIs);    

all_pSF = struct('vox_indices', struct_size, ...
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

%% Estimate pSF

HIRF = defineHRF();

total_elapsed_time = 0;

for subj = 1:num_subjs
    for roi = 1:num_ROIs

        if toggles.disp_on, disp(['Estimating pSF for subject ' num2str(subj) ', V' num2str(roi) ' ...']); end

        tic;    

        pSF = estimatePSF(sample_data.measured_BOLD{subj}(:,:,roi), sample_data.I{subj}, HIRF, p, toggles);
        all_pSF(subj,roi) = pSF;

        elapsed_time = round(toc/60,1);

        if toggles.disp_on, disp(['Elapsed time: ~' num2str(elapsed_time) ' minute(s).']); disp(' '); end

        total_elapsed_time = total_elapsed_time + elapsed_time;

    end
end

if toggles.disp_on, disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']); end
if toggles.parallelization, delete(gcp); end
if save_pSF, save('estimates/all_pSF.mat', 'all_pSF'); end

%% Plot voxel fits

if make_voxel_plots

    num_voxels = 5;

    fg = figure('Visible','on','Color','w');
    set(0,'CurrentFigure',fg);

    for subj = 1:num_subjs
        for roi = 1:num_ROIs
            for vox = 1:num_voxels

                %% Spatial frequency tuning curve

                figure_name = ['S' num2str(subj) ' V' num2str(roi) ' Vox #' num2str(vox) ' pSF'];
                
                semilogx(p.sfs, all_pSF(subj,roi).est_SFT(vox,:), 'k');
            
                % Format figure
                xlabel('log[SF] (cpd)'); ylabel('R');
                xlim([p.sf_min p.sf_max]); ylim([0 1]);
                xticks([p.sf_min 0.5 1 5 p.sf_max]); xticklabels([p.sf_min 0.5 1 5 p.sf_max]); 
                set(gca,'TickDir','out'); box off;
                title(['\mu = ' num2str(round(all_pSF(subj,roi).param_est(vox,1),2)) ' | \sigma = ' num2str(round(all_pSF(subj,roi).param_est(vox,2),2))])
                
                if save_voxel_plots, 
                    saveas(gcf,['figures/' figure_name '.png']); clf; 
                    if toggles.disp_on, disp(['Saved ' figure_name '.png in /figures']); end
                end

                %% Estimated voxel time series

                figure_name = ['S' num2str(subj) ' V' num2str(roi) ' Vox #' num2str(vox) ' BOLD Time Series'];
                
                plot(all_pSF(subj,roi).measured_BOLD(vox,:),'k');
                hold on
                plot(all_pSF(subj,roi).est_BOLD(vox,:),'r');
                
                % Format figure
                xlabel('Time (s)'); ylabel('BOLD (% change)');
                set(gca,'TickDir','out'); box off;
                title(['R^2 = ' num2str(round(all_pSF(subj,roi).r2(vox),2))])
                legend({'measured','estimate'})
    
                if save_voxel_plots, 
                    saveas(gcf,['figures/' figure_name '.png']); clf; 
                    if toggles.disp_on, disp(['Saved ' figure_name '.png in /figures']); end
                end
   
            end
        end
    end

    close(fg);

end