% createTextures - Generate stimuli for acquiring pSF
%   Generates a set of bandpass filtered noise textures for a given image size, pixel size, and ppd.
%
%   Syntax
%       textures = createTextures(radius_px, contrast, noise_sample_count, ppd, px_size, save_textures, texture_filepath)
%
%   Input Arguments
%       radius_px – radius of the stimulus in pixels
%       contrast – contrast of the noise
%       noise_sample_count – number of noise samples
%       ppd – pixels per degree
%       px_size – size of each pixel in cm
%       save_textures – boolean to save the textures
%       texture_filepath – path to save the textures
%
%   Output Arguments
%       textures – [height, width, filter_count, noise_sample] of bandpass filtered noise textures
%       filters – structure containing filter parameters (count, min, max, width, gauss_smoothening_sd, centers, lower_bound, upper_bound, f_Nyquist, masks)

function [textures, filters] = createTextures(radius_px, contrast, noise_filter_count, noise_sample_count, ppd, px_size, save_textures, texture_filepath)
    
    %% Bandpass filter parameters

    sf_min = 0.5;
    sf_max = 12;
    width = 0.1;
    gauss_smoothening_sd = 0.1 * ppd;

    centers = 10.^linspace(log10(sf_min), log10(sf_max), noise_filter_count);
    lower_bound = centers - width/2;
    upper_bound = centers + width/2;
    f_Nyquist = 1 / (2 * px_size);

    %% Create bandpass filtered noise textures
    
    textures = nan(radius_px, radius_px, noise_filter_count, noise_sample_count);
    masks = nan(radius_px, radius_px, noise_filter_count);

    for n_filter = 1:noise_filter_count

        bandpass_filter = Bandpass2(radius_px, lower_bound(n_filter)/f_Nyquist, upper_bound(n_filter)/f_Nyquist);
        bandpass_filter = imgaussfilt(bandpass_filter, gauss_smoothening_sd);
        bandpass_filter = (bandpass_filter - min(bandpass_filter(:))) ./ (max(bandpass_filter(:)) - min(bandpass_filter(:)));
        
        masks(:, :, n_filter) = bandpass_filter;
        
        for n_noise_sample = 1:noise_sample_count

            noise_texture = 2 * rand(radius_px) - 1; 
            fft_texture = fftshift(fft2(noise_texture));

            filtered_texture = fft_texture .* bandpass_filter; 
            filtered_texture = ifft2(ifftshift(filtered_texture));
            filtered_texture = real(filtered_texture);
            filtered_texture = filtered_texture - mean(filtered_texture(:));

            max_value = max(abs(filtered_texture(:)));

            if max_value > 0
                filtered_texture = filtered_texture / (2 * max_value);
            end

            filtered_texture = 127 * (1 + filtered_texture * 2 * contrast);

            textures(:, :, n_filter, n_noise_sample) = filtered_texture;

        end

    end

    %% Store filters

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

    if nargin > 6
        if save_textures
            save(texture_filepath, 'textures', 'filters'); 
        end
    end

end