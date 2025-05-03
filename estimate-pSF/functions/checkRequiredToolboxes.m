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
