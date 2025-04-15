% logGauss - Generate neural response from log-Gaussian function
%
%   Syntax
%       R = logGauss(params, x)
%
%   Input Arguments
%       params – [mu, sigma]
%       x – time-series of spatial frequencies
%
%   Output Arguments
%       R – neural response

function R = logGauss(params, x)

    mu = params(1);
    sigma = params(2);

    R = exp( -( (log(x) - log(mu)).^2 ) ./ (2 * sigma.^2) );

end

