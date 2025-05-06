function [bandwidth_oct, sf_halfmax] = cpd2oct(pSFT, sf)

%% Finds the lowest and highest SF in cycles/degree (cpd) that produced half the maximum response.

max_pSFT = max(pSFT);
sf_halfmax(1) = sf(find(pSFT >= max_pSFT/2, 1, 'first'));
sf_halfmax(2) = sf(find(pSFT >= max_pSFT/2, 1, 'last'));

%% Convert cpd to octaves
% To convert cpd to octaves, a ratio is created between the higher and lower SF. Then,
% this ratio is log transformed at log base 2, because an increase in SF in octave units doubles the frequency in cpd.
% log2(SF_H/SF_L) = log2(SF_H) - log2(SF_L) = Full-width at half maxx

bandwidth_oct = log2(sf_halfmax(2)/sf_halfmax(1));

%{

figure
plot(sf,pSFT); hold on;
line([sf_halfmax(1) sf_halfmax(1)],[0 max(pSFT)],'Color',[0 1 0],'LineStyle','--');
line([sf_halfmax(2) sf_halfmax(2)],[0 max(pSFT)],'Color',[0 1 0],'LineStyle','-');
set(gca,'XScale','log');

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