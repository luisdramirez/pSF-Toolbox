% checkRequiredToolboxes - Checks for required MATLAB toolboxes.
%   If any are missing, an error is returned.

function checkRequiredToolboxes(toggles)
    
    if toggles.parallelization
        required_toolboxes = {'Optimization Toolbox', ...
            'Parallel Computing Toolbox'};
    else
        required_toolboxes = {'Optimization Toolbox'};
    end
    
    % Get list of installed toolboxes
    installed_toolboxes = arrayfun(@(x) x.Name, ver, 'UniformOutput', false);
    
    % Identify missing toolboxes
    missing = setdiff(required_toolboxes, installed_toolboxes);
    
    if ~isempty(missing)
        error('The following required MATLAB toolboxes are missing:\n- %s\n\nPlease install them using the Add-On Explorer or contact your administrator.', ...
            strjoin(missing, '\n- '));
    else
        disp('All required MATLAB toolboxes are installed.');
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