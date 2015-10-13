function [fmin,out_param] = quadratic(a,b);
%
% >> [fmin,out_param] = quadratic(-2,3);
% >> abs(fmin+8) < 1e-6, abs(mean(out_param.intervals)+2) 
%    ans = 1 ans = 0
% >> out_param.intervals
%  ans =
%    -2     
%    -2     
%
% >> [fmin,out_param] = quadratic(-2,4);
% >> abs(fmin+8) < 1e-6, abs(mean(out_param.intervals)+2) 
%    ans = 1 ans = 0
% >> out_param.intervals
%  ans =
%    -2     4
%    -2     4
%

format compact

f=@(x) -(x-1).^2+1; 
[fmin,out_param] = funmin_g_CSC(f,a,b,1e-7,1e-4,10,10,1000000);
%funmin_g_demo(fmin, out_param)
