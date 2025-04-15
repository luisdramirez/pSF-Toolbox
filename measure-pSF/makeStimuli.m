% makeStimuli - Generate stimuli for acquiring pSF
%   Generates a set of bandpass filtered noise textures for a given image size, pixel size, and ppd.
%
%   Syntax
%       textures = makeStimuli(image_size_px, ppd, px_size)
%
%   Input Arguments
%       image_size_px – [height, width] of the image in pixels
%       ppd – pixels per degree
%       px_size – size of each pixel in degrees
%
%   Output Arguments
%       textures – [height, width, filter_count, noise_sample] of bandpass filtered noise textures
%       stimulus – structure containing stimulus parameters
%       filter – structure containing filter parameters and masks

function [textures, stimulus, filter] = makeStimuli(image_size_px, ppd, px_size)
    
    %% Stimulus parameters

    stimulus.contrast = 1;
    stimulus.noise_sample_count = 10;
    stimulus.size_px = image_size_px;

    gray = 127;

    %% Bandpass filter parameters

    filter.count = 40;
    filter.min = 0.1;
    filter.max = 12;
    filter.width = 0.1;
    filter.gauss_sd = 0.1*ppd;

    filter.centers = 10.^linspace(log10(filter.min), log10(filter.max), filter.count);
    filter.lower_bound = filter.centers - filter.width;
    filter.upper_bound = filter.centers + filter.width;
    filter.f_Nyquist = 1/(2 * px_size);

    filter.masks = nan(image_size_px(1), image_size_px(2), filter.count);

    %% Create bandpass filtered noise textures
    
    textures = nan(image_size_px(1), image_size_px(2), filter.count, stimulus.noise_sample_count);
    
    for n_filter = 1:filter.count

        sf_bp_filter = Bandpass2(image_size_px, filter.lower_bound(n_filter)/filter.f_Nyquist, filter.upper_bound(n_filter)/filter.f_Nyquist);
        sf_bp_filter = imgaussfilt(sf_bp_filter, filter.gauss_sd);
        sf_bp_filter_norm = (sf_bp_filter - min(sf_bp_filter(:))) ./ (max(sf_bp_filter(:)) - min(sf_bp_filter(:)));
        
        filter.masks(:, :, n_filter) = sf_bp_filter_norm;
        
        for n_noise_sample = 1:stimulus.noise_sample_count

            init_noise = 2*rand(image_size_px)-1; 
            fft_noise = fftshift(fft2(init_noise));

            filtered_noise = filter.masks(:, :, n_filter) .* fft_noise; 
            filtered_noise = ifft2(ifftshift(filtered_noise));
            filtered_noise = (filtered_noise - min(filtered_noise(:)))./max(filtered_noise(:));
            filtered_noise = abs(filtered_noise * stimulus.contrast * gray + gray);

            textures(:, :, n_filter, n_noise_sample) = filtered_noise;

        end

    end

end