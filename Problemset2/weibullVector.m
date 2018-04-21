%Returns the value of the Weibull distribution for a certain value of u or
%a vector containing multiple u values.
function y = weibullVector(u,labda,k)

y = k/labda .* (u ./ labda) .^(k-1).*exp( -(u./labda).^k);


end

