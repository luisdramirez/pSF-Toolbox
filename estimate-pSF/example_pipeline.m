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

load(fullfile(data_dir, 'sample_data.mat'));

num_subjs = length(sample_data);
num_ROIs = size(sample_data(1).measured_BOLD, 3); % assumes all subjs have the same num ROIs

%% Initialize pSFT struct

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

%% Estimate pSFT

HIRF = defineHRF();

total_elapsed_time = 0;

for subj = 1:num_subjs
    for roi = 1:num_ROIs

        if toggles.disp_on, disp(['++++ S' num2str(subj) ' V' num2str(roi) ' ++++']); end

        tic;    

        pSFT = estimatePSF(sample_data(subj).measured_BOLD(:,:,roi), sample_data(subj).I, HIRF, p, toggles);
        all_pSF(subj,roi) = pSFT;

        elapsed_time = round(toc/60,1);

        if toggles.disp_on, disp(['Elapsed time: ~' num2str(elapsed_time) ' minute(s).']); disp(' '); end

        total_elapsed_time = total_elapsed_time + elapsed_time;

    end
end

if toggles.disp_on, disp(['Total elapsed time: ~' num2str(total_elapsed_time) ' minutes']); end
if toggles.parallelization, delete(gcp); end

if save_pSF
    curr_time = string(datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss'));
    save([save_dir '/all_pSF_n' num2str(num_subjs) '_' char(curr_time) '.mat'], 'all_pSF');
    if toggles.disp_on, disp(['Saved all_pSF.mat in /' save_dir]); end
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
                
                semilogx(p.sfs, all_pSF(subj,roi).est_SFT(:, vox), 'k');
            
                % Format figure
                xlabel('log[SF] (cpd)'); ylabel('R');
                xlim([p.sf_min p.sf_max]); ylim([0 1]);
                xticks([p.sf_min 0.5 1 5 p.sf_max]); xticklabels([p.sf_min 0.5 1 5 p.sf_max]); 
                set(gca,'TickDir','out'); box off;
                title(['\mu = ' num2str(round(all_pSF(subj,roi).param_est(1,vox),2)) ' | \sigma = ' num2str(round(all_pSF(subj,roi).param_est(2,vox),2))])
                
                if save_voxel_plots 
                    saveas(gcf,[figure_path '/' figure_name '.png']); clf; 
                    if toggles.disp_on, disp(['Saved ' figure_name '.png in /' figure_path]); end
                end

                %% Voxel time series

                figure_name = ['Vox #' num2str(vox) ' BOLD Time Series'];
                
                plot(all_pSF(subj,roi).measured_BOLD(:,vox),'k');
                hold on
                plot(all_pSF(subj,roi).est_BOLD(:,vox),'r');
                
                % Format figure
                xlabel('Time (s)'); ylabel('BOLD (% change)');
                set(gca,'TickDir','out'); box off;
                title(['R^2 = ' num2str(round(all_pSF(subj,roi).r2(vox),2))])
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