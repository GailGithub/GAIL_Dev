function [fmin,out_param] = quadratic;
%
% >> [fmin,out_param] = quadratic;
% >> abs(fmin+4.289999999999999) < 1e-6, abs(mean(out_param.intervals)+2) 
%    ans = 1 ans = 0

format compact

f=@(x) -(x-0.3).^2+1; 
[fmin,out_param] = funmin_g_CSC(f,-2,2,1e-7,1e-4,10,10,1000000);
%funmin_g_demo(fmin, out_param)
