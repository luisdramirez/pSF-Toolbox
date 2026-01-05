% stimulusParams - Define stimulus parameters
%   Define stimulus radius, contrast, noise filter count, and noise sample count
%
% Syntax
%   p = stimulusParams(p, w)
%
% Input Arguments
%   p – structure containing scan parameters
%   w – structure containing window parameters
%
% Output Arguments
%   p – updated parameters structure

function p = stimulusParams(p, w)

p.aperture_radius_deg = 9; % defines the apparent stimulus size
p.aperture_radius_px = round(p.aperture_radius_deg * w.ppd);
p.stimulus_radius_px = round(p.aperture_radius_px * 1.25); % default = 1.1

p.stimulus_contrast = 0.9; % default = 0.9
p.noise_filter_count = 40; % default = 40
p.noise_sample_count = 10; % default = 10

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