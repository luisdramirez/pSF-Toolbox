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
