function [ vanillaCall, vanillaPut ] = vanillaPrice( S0, K, r, sigma, T )
% vanillaPrice is a function that calculates the price of vanilla option
d1 = log(S0/K)+(r+(sigma^2)/2)*T;
d2 = d1-sigma*sqrt(T);
vanillaCall = S0*normcdf(d1) - K*exp(-r*T)*normcdf(d2);
vanillaPut = K*exp(-r*T) - S0 + vanillaCall;
end

