% generateBOLD - Generate BOLD response from neural response
%   Convolves neural response with HRF and applies amplitude and baseline scaling.
%
%   Syntax
%       BOLD = generateBOLD(R, HRF, beta, beta_0)
%
%   Input Arguments
%       R      - Neural response time series [time x 1]
%       HRF    - Hemodynamic response function [time x 1]
%       beta   - BOLD amplitude
%       beta_0 - BOLD baseline offset
%
%   Output Arguments
%       BOLD   - BOLD response time series [time x 1]

function BOLD = generateBOLD(R, HRF, beta, beta_0)

% Convolve neural response with HRF
BOLD = conv(R, HRF);

% Truncate to original length
BOLD = BOLD(1:length(R));

% Apply amplitude and baseline scaling
BOLD = beta .* ((BOLD ./ mean(BOLD)) - 1) + beta_0;

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

