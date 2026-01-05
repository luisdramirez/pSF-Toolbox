% generateSFTimeSeries - Generate stimulus time series from run info
%   Converts SF indices into a vector of SF values.
%   Blank periods are filled with a small value (0.0001) to avoid log(0) in logGauss.
%
%   Syntax
%       I = generateSFTimeSeries(p, t)
%
%   Input Arguments
%       p       - structure containing scan parameters (including p.I_idx with SF indices)
%       t       - structure containing timing parameters
%
%   Output Arguments
%       I - [num_TRs x 1] stimulus time series of SF values (cpd)

function I = generateSFTimeSeries(p, t, filters)

%% Parameters

blank_sf = 0.0001; % Small value to avoid log(0) in logGauss

% Calculate durations in TRs
blank_TRs = round(t.blank_period_dur / t.TR);
block_TRs = round(t.block_dur / t.TR);
num_TRs = round(t.scan_dur / t.TR);

%% Build stimulus time series

I = nan(num_TRs, 1);
TR_idx = 1;

for n_block = 1:p.num_blocks

    % Add blank period before each block
    I(TR_idx : TR_idx + blank_TRs - 1) = blank_sf;
    TR_idx = TR_idx + blank_TRs;

    % Grab presented SFs
    block_sf = filters.centers(p.I_idx(:, n_block));
    I(TR_idx : TR_idx + block_TRs - 1) = block_sf;
    TR_idx = TR_idx + block_TRs;

end

% Add final blank period
I(TR_idx : TR_idx + blank_TRs - 1) = blank_sf;

%% Validate

if any(isnan(I))
    warning('Stimulus time series contains NaN values. Check timing parameters.');
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
