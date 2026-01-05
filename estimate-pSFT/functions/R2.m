% R2 - Calculate coefficient of determination (R^2)
%   Computes R^2 between measured and estimated time series.
%
%   Syntax
%       r2 = R2(measured, estimated)
%
%   Input Arguments
%       measured  - Measured time series [time x 1]
%       estimated - Estimated time series [time x 1]
%
%   Output Arguments
%       r2 - Coefficient of determination (R^2)

function r2 = R2(measured, estimated)

SS_res = SSE(measured, estimated);
SS_tot = SSE(measured, mean(measured));

r2 = 1 - (SS_res / SS_tot);

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

