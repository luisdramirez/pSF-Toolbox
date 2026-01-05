% SSE - Calculate sum of squared errors (SSE)
%   Computes SSE between measured and estimated time series.
%
%   Syntax
%       sse = SSE(measured, estimated)
%
%   Input Arguments
%       measured  - Measured time series [time x 1]
%       estimated - Estimated time series [time x 1]
%
%   Output Arguments
%       sse - Sum of squared errors

function sse = SSE(measured, estimated)

sse = sum((measured - estimated).^2);

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


