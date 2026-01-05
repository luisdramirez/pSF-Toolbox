% createTextures - Generate stimuli for acquiring pSF
%   Generates a set of bandpass filtered noise textures for a given image size and ppd.
%
% Syntax
%   [textures, filters] = createTextures(radius_px, contrast, noise_filter_count, noise_sample_count, ppd, save_textures, texture_filepath)
%
% Input Arguments
%   radius_px – radius of the stimulus in pixels
%   contrast – contrast of the noise
%   noise_filter_count – number of noise filters
%   noise_sample_count – number of noise samples
%   ppd – pixels per degree
%   save_textures – boolean to save the textures
%   texture_filepath – path to save the textures
%   verify_energy – boolean to verify the energy of the textures
%
% Output Arguments
%   textures – [height, width, filter_count, noise_sample] of bandpass filtered noise textures
%   filters – structure containing filter parameters (count, min, max, width, gauss_smoothening_sd, centers, lower_bound, upper_bound, f_Nyquist, masks)

function [textures, filters] = createTextures(radius_px, contrast, noise_filter_count, noise_sample_count, ppd, save_textures, texture_filepath)

%% Bandpass filter parameters

% Define spatial frequency range
sf_min = 0.5; % cycles per degree (cpd) (default = 0.5)
sf_max = 12; % cpd (default = 12)
width = 0.1; % cpd (default = 0.1)
gauss_smoothening_sd = 1; % pixels, Small amount of smoothing to avoid ringing but preserve bandwidth

% Logarithmically space the center frequencies
centers = 10.^linspace(log10(sf_min), log10(sf_max), noise_filter_count); % cpd

% Define upper and lower bounds of each SF band
lower_bound = centers - width/2; % cpd
upper_bound = centers + width/2; % cpd

% Calculate Nyquist frequency
% In this case, the Nyquist frequency represents the minimum number of pixels needed for half a cycle from light to dark.
f_Nyquist = ppd / 2; % angular

% Calculate lower and upper bounds
f_low = lower_bound / f_Nyquist; % Normalized frequency (0 to 1, where 1 = Nyquist)
f_high = upper_bound / f_Nyquist;

%% Create bandpass filtered noise textures

% Initialize texture and mask arrays
textures = nan(radius_px, radius_px, noise_filter_count, noise_sample_count);
masks = nan(radius_px, radius_px, noise_filter_count);

% Loop through each filter (spatial frequency band)
for n_filter = 1:noise_filter_count

    % Create 2D bandpass filter mask
    bandpass_filter = Bandpass2(radius_px, f_low(n_filter), f_high(n_filter));

    % Smooth the filter to avoid ringing artifacts
    bandpass_filter = imgaussfilt(bandpass_filter, gauss_smoothening_sd);

    % Clip filter values to [0, 1] range (smoothing may cause minor overshoots)
    bandpass_filter = max(0, min(1, bandpass_filter));

    % figure, imagesc(bandpass_filter), colormap("gray"), axis square, axis off, box off, title("1: bandpass filter")

    % Store the filter mask
    masks(:, :, n_filter) = bandpass_filter;

    % Generate noise samples for this frequency band
    for n_noise_sample = 1:noise_sample_count

        % Generate white noise
        noise_texture = 2 * rand(radius_px) - 1;

        % figure, imagesc(noise_texture), colormap("gray"), axis square, axis off, box off, title("2a: noise sample") 

        % Convert to frequency domain
        fft_texture = fftshift(fft2(noise_texture));

        % figure, imagesc(abs(fft_texture)), colormap("gray"), axis square, axis off, box off, title("2b: fft noise")

        % Apply the bandpass filter
        bp_fft_texture = fft_texture .* bandpass_filter;

        % figure, imagesc(abs(bp_fft_texture)), colormap("gray"), axis square, axis off, box off, title("3: bandpass fft noise")

        % Convert back to spatial domain
        filtered_texture = ifft2(ifftshift(bp_fft_texture));
        filtered_texture = real(filtered_texture);
        filtered_texture = filtered_texture - mean(filtered_texture(:)); % Remove DC component

        % Normalize amplitude
        max_value = max(abs(filtered_texture(:)));

        if max_value > 0
            filtered_texture = filtered_texture / (2 * max_value);
        end

        % Apply contrast scaling and shift to mean gray level
        filtered_texture = 127 * (1 + filtered_texture * 2 * contrast);

        % figure, imagesc(filtered_texture), colormap("gray"), axis square, axis off, box off, title("4: bandpass-filtered noise") 

        % Store the final texture
        textures(:, :, n_filter, n_noise_sample) = filtered_texture;

    end

end

%% Store filters

% Package filter parameters into a structure
filters.sf_min = sf_min;
filters.sf_max = sf_max;
filters.width = width;
filters.gauss_smoothening_sd = gauss_smoothening_sd;
filters.centers = centers;
filters.lower_bound = lower_bound;
filters.upper_bound = upper_bound;
filters.f_Nyquist = f_Nyquist;
filters.masks = masks;

%% Save textures

% Save the generated textures and filter info to a file if requested
if nargin > 6
    if save_textures
        save(texture_filepath, 'textures', 'filters');
    end
end

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