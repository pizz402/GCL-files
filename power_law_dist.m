function [ki] = power_law_dist(k0,n,gamma)

% p(k) = k^(-gamma)

kf = n^(1/(gamma-1)); % natural cutoff

ki = round( k0 * (1 + ((kf/k0).^(1-gamma)-1)*rand(1,n)).^(1/(1-gamma)) );


end

