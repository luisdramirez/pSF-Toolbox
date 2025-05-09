% defineHRF - Define the canonical hemodynamic impulse response function (HRF)
%   Boynton & Heeger, 1996 Journal of Neuroscience 
%
%   Syntax
%       HRF = defineHRF()
%
%   Output Arguments
%       HRF – hemodynamic impulse response function

function HRF = defineHRF()

    %% HRF parameters
    
    delta = 2.05; % delay between stimulus onset and fMRI response (Boynton & Heeger, 1998)
    tau = 1.08; % time constant (Boynton & Heeger, 1998)
    n_HRF = 3;  % phase delay (Boynton & Heeger, 1998)
    B = 0;  % baseline BOLD activity, offset from 0 (in reality, this is would be some random big number)
    tmp_t = 0:20;
    t_shift = max(tmp_t-delta,0);
    
    %% Calculate HRF

    HRF = (((t_shift/tau).^(n_HRF-1)) .* exp(-(t_shift/tau))) / (tau*(factorial(n_HRF-1))) + B; % gamma function
    HRF = HRF./max(HRF); % normalize the hypothetical HRF so that it peaks at 1
    
    %% Plot HRF

    %{
    figure('Color', [1 1 1], 'Name', 'HRF')
    plot(tmp_t, HRF,'k'); 
    xlim([0 25]); 
    xlabel('Time (s)');  ylabel('Arbitrary M.R. units');
    box off; axis square; set(gca, 'TickDir','out'); 
    %}

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