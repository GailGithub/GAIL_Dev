function [ x ] = stdnorminv( p )
% this function is the inverse function of CDF of standard normal distribution
x = -sqrt(2).*erfcinv(2*p);
end

