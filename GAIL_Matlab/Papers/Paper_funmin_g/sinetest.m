function [fmin,out_param] = sinetest(a,b);
%
% >> [fmin,out_param] = sinetest(0,1);
% >> abs(fmin-0) < 1e-6, abs(mean(out_param.intervals)) 
%    ans = 1 ans = 0
%
%
% >> [fmin,out_param] = sinetest(0,pi);
% >> abs(fmin-0) < 1e-6
%    ans = 1 
% >> x = abs(mean(out_param.intervals))
%    x = 0  3.1416
%
%
format compact


[fmin,out_param]=funmin_g_CSC(@(x) sin(x),a,b);
funmin_g_demo(fmin, out_param)