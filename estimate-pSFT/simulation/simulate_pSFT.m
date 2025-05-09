% simulate_pSFT

close all; 
clear variables; 
clc;

addpath('../functions');

%% Settings

% Plot settings
font_type = 'Helvetica';
axes_label_font_size = 14;
axes_tick_font_size = 13;
tick_length = 0.020;
line_width = 1;

% Colors (RGB)
red = [204 0 0]/255;
green = [0 153 0]/255;
blue = [0 76 152]/255;
purple = [102 51 204]/255;
black = [0 0 0];
white = [1 1 1];
gray = white/2;

%% Define pSFT model

pSFT_model = @(params,sf) logGauss(params, sf);

%% Spatial frequency parameters

sf_min = 0.1;
sf_max = 12;
sf_count = 1000;
log_sf = logspace(log10(sf_min), log10(sf_max), sf_count); 

%% Baseline pSFT

% pSFT parameters
mu = 1;
sigma = 0.5;

% Generate neural response with log gaussian function handle 
params = [mu, sigma];
R = pSFT_model(params, log_sf);

% Plot pSFT
figure('Name','Baseline pSFT','Color','w')

plot(log_sf, R, 'Color', green);

% Format figure
xlabel('Spatial frequency (cpd)','FontName',font_type,'FontSize',axes_label_font_size); 
ylabel('Response (au)','FontName',font_type,'FontSize',axes_label_font_size);
set(gca, 'TickDir','out', 'XScale','log','TickLength',[tick_length tick_length],'FontName',font_type,'FontSize',axes_tick_font_size);
xticks([sf_min 1 5 sf_max]); xticklabels([sf_min 1 5 sf_max]); xlim([sf_min sf_max])
yticks([0 0.5 1])
box off;

%% Change in peak

mu_inc = 2;
params_inc_mu = [mu * mu_inc, sigma];
R_mu(1,:) = pSFT_model(params_inc_mu, log_sf);

mu_dec = 2;
params_dec_mu = [mu / mu_dec, sigma];
R_mu(2,:) = pSFT_model(params_dec_mu, log_sf);

% Plot responses 
figure('Name', 'Change in peak','Color','w')

plot(log_sf,R,'Color',green); hold on;

plot(log_sf,R_mu(1,:),'Color',red) % increase
plot(log_sf,R_mu(2,:),'Color',blue) % decrease

line([params(1) params(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', green);
line([params_inc_mu(1) params_inc_mu(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', red);
line([params_dec_mu(1) params_dec_mu(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', blue);

% Format figure
xlabel('Spatial frequency (cpd)'); ylabel('Response (au)');
set(gca, 'TickDir','out', 'XScale','log','TickLength',[tick_length tick_length]);
if mu ~= 1
    xticks([sf_min 1 mu sf_max]); xticklabels([sf_min 1 mu sf_max]);
else
    xticks([sf_min mu sf_max]); xticklabels([sf_min mu sf_max]);
end
xlim([sf_min sf_max])
xlim([sf_min sf_max])
ylim([0 1])
box off;
legend({'Baseline pSFT', [num2str(mu_inc) ' \times \mu'], ['\mu / ' num2str(mu_dec)]})

%% Change in bandwidth

sigma_inc = 2;
params_inc_sigma = [mu, sigma * sigma_inc];
R_sigma(1,:) = pSFT_model(params_inc_sigma, log_sf);

sigma_dec = 2;
params_dec_sigma = [mu, sigma / sigma_dec];
R_sigma(2,:) = pSFT_model(params_dec_sigma, log_sf);

% Plot responses 
figure('Name', 'Change in bandwidth','Color','w')

plot(log_sf,R,'Color',green); hold on;

plot(log_sf, R_sigma(1,:),'Color',red) % increase
plot(log_sf, R_sigma(2,:),'Color',blue) % decrease

line([params(1) params(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', green);
line([params_inc_sigma(1) params_inc_sigma(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', red);
line([params_dec_sigma(1) params_dec_sigma(1)], [min(ylim) max(ylim)],'LineStyle','--', 'Color', blue);

% Format figure
xlabel('Spatial frequency (cpd)'); ylabel('Response (au)');
set(gca, 'TickDir','out', 'XScale','log','TickLength',[tick_length tick_length]);
if mu ~= 1
    xticks([sf_min 1 mu sf_max]); xticklabels([sf_min 1 mu sf_max]);
else
    xticks([sf_min mu sf_max]); xticklabels([sf_min mu sf_max]);
end
xlim([sf_min sf_max])
xlim([sf_min sf_max])
ylim([0 1])
box off;
legend({'Baseline pSFT', [num2str(sigma_inc) ' \times \sigma'], ['\sigma / ' num2str(sigma_dec)]})

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