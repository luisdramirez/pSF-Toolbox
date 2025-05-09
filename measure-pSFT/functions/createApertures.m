% createApertures
%   Creates a stimulus and fixation apertures.
%
% Syntax
%   [stimulus_aperture, fixation_aperture] = createApertures(p, w)
%
%   Input Arguments
%       p – structure containing scan parameters
%       w – structure containing window parameters
%
%   Output Arguments
%       stimulus_aperture – stimulus aperture
%       fixation_aperture – fixation aperture

function [stimulus_aperture, fixation_aperture] = createApertures(p, w)

    %% Stimulus aperture

    % Create cartesian coordinates for texture: -x to x, -y to y
    [x,y] = meshgrid(-w.screen_width_px/2:w.screen_width_px/2-1,-w.screen_height_px/2:w.screen_height_px/2-1);
    
    % Calculate eccentricity (e.g., radius) of each pixel coordinate in the cartesian grid
    r = sqrt((x).^2 + (y).^2); % solving for r in: x^2 + y^2 = r^2
        
    % Set values inside desired eccentricity to be fully transparent
    aperture_mask = ones(w.screen_height_px, w.screen_width_px);
    aperture_mask(r <= p.aperture_radius_px) = 0;
    
    % Create aperture texture
    stimulus_aperture = nan(w.screen_height_px, w.screen_width_px, 2);
    stimulus_aperture(:,:,1) = ones(size(aperture_mask)) * w.gray;
    stimulus_aperture(:,:,2) = imgaussfilt(aperture_mask,  0.1*w.ppd) * 255;

    %% Fixation aperture

    [x, y] = meshgrid(-(p.fixation_aperture_px/2):(p.fixation_aperture_px/2) -1, -(p.fixation_aperture_px/2):(p.fixation_aperture_px/2)-1);
    r = sqrt(x.^2+y.^2); 

    aperture_mask = ones(p.fixation_aperture_px); 
    aperture_mask(r >= (p.fixation_aperture_radius_px/2)) = 0; 
    
    fixation_aperture = nan(p.fixation_aperture_px, p.fixation_aperture_px, 2); 
    fixation_aperture(:,:,1) = ones(p.fixation_aperture_px) * w.gray; 
    fixation_aperture(:,:,2) = imgaussfilt(aperture_mask, 1) * 255; 

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