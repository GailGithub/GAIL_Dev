function [ p ] = stdnormcdf( z )
% this function is to define cumulative distribution function (CDF) of the
% standard normal distribution.
p = 0.5*erfc(-z./sqrt(2));
% Use the complementary error function, rather than .5*(1+erf(z/sqrt(2))),
% to produce accurate near-zero results for large negative x.
end

