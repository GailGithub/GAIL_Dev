function [fmin,out_param] = sinetest;
%
% >> [fmin,out_param] = sinetest;
% >> abs(fmin-0) < 1e-6, abs(mean(out_param.intervals)) 
%    ans = 1 ans = 0

format compact


[fmin,out_param]=funmin_g_CSC(@(x) sin(x));
funmin_g_demo(fmin, out_param)