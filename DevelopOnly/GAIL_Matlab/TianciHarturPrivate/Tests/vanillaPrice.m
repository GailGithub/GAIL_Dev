function [ vanillaCall, vanillaPut ] = vanillaPrice( S0, K, r, sigma, T, weight )
% vanillaPrice is a function that calculates the price of vanilla option
d1 = (log(S0/K)+T.*(r+(sigma.^2)/2))./(sigma.*sqrt(T));
d2 = (log(S0/K)+T.*(r+(sigma.^2)/2))./(sigma.*sqrt(T))-(sigma.*sqrt(T));
vanillaCall = weight*(S0.*normcdf(d1) - K.*exp(-r.*T).*normcdf(d2));
vanillaPut = weight*(K*exp(-r.*T)-S0+(S0.*normcdf(d1) - K*exp(-r.*T).*normcdf(d2)));
end

