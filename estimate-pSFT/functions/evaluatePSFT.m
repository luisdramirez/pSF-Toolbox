% evaluatePSFT - Evaluate goodness of fit of pSFT parameters
%   Evaluates the sum of squares error (SSE) between the measured and estimated BOLD.
%   1. Neural response, R, is generated from logGauss()
%   2. R is convolved with the HRF to generate BOLD time series
%   3. Generated BOLD is converted to percent signal change [ ( y/mean(y) ) - 1 ]
%   4. Generated BOLD is then compared to measured BOLD
%
%   Syntax
%       sse = evaluatePSFT(free_params, fixed_params)
%
%   Input Arguments
%       free_params – [mu, sigma, beta, beta_0]
%       fixed_params – {I, measured_BOLD, HRF}
%
%   Output Arguments
%       sse – sum of squares error

function sse = evaluatePSFT(free_params, fixed_params)

%% Generate neural response

I = fixed_params{1};

R = logGauss(free_params, I);
R(isnan(R)) = 0;

%% Generate BOLD response

HRF = fixed_params{3};
beta = free_params(3);
beta_0 = free_params(4);

est_BOLD = generateBOLD(R, HRF, beta, beta_0);

%% Get goodness of fit

measured_BOLD = fixed_params{2};
sse = SSE(measured_BOLD, est_BOLD);

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