function checkPTB()
    if ~exist('PsychtoolboxVersion', 'file')
        error(['Psychtoolbox-3 is not installed or not on the MATLAB path.\n' ...
               'Please install it from http://psychtoolbox.org/download.']);
    else
        v = PsychtoolboxVersion;
        fprintf('Psychtoolbox-3 is installed. Version: %s\n', v);
    end
end
