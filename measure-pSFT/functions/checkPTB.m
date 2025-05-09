% checkPTB - Check if Psychtoolbox is installed
%
% Syntax
%   checkPTB()
%

function checkPTB()
    if ~exist('PsychtoolboxVersion', 'file')
        error(['Psychtoolbox is not installed or not on the MATLAB path.\n' ...
               'Please install it from http://psychtoolbox.org/download.']);
    else
        v = PsychtoolboxVersion;
        fprintf('Psychtoolbox is installed. Version: %s\n', v);
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