function [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(m);
% The Fejer-Jackson inequality says that the following partial sum is 
% positive for a given positive integer n and for x in (0, pi): 
%     fnx = @(n,x)((1./(1:n)) * sin((1:n)' * x));
%
% Example 1:
% >> [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(32);
% >>  abs(fmin) < 1e-6, abs(fmax-1.804097) < 1e-6
%  ans = 1  ans = 1
% >>  funmin_g_demo(fmin, out_min); hold on; funmin_g_demo(fmax, out_max); hold off;
% >>  out_min.intervals
%    ans =
%         0    3.1401
%         0    3.1416
%
%
% Example 2:
% >> [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(128);
% >>  abs(fmin) < 1e-6, abs(fmax-1.839745) < 1e-6
%  ans = 1  ans = 1
% >>  funmin_g_demo(fmin, out_min); hold on; funmin_g_demo(fmax, out_max); hold off;
%
%
% Example 3:
% >> [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(512);
%     Warning***
% >>  abs(fmin) < 1e-6, abs(fmax-1.848874) < 1e-6
%  ans = 1  ans = 1
% >>  funmin_g_demo(fmin, out_min); hold on; funmin_g_demo(fmax, out_max); hold off;
%
%
% Reference: Fejer-Jackson inequality, Nick Trefethen
%         http://www.chebfun.org/examples/fourier/FejerJackson.html


fnx = @(n,x)((1./(1:n)) * sin((1:n)' * x));
f = @(x) fnx(m,x);
a = 0; b = pi;
in_param.abstol = 1e-6;
in_param.TolX = 1e-6;
[fmin,out_min]=funmin_g(f,a,b,in_param); 

g = @(x) fnx(m,x)*-1;
[fmax,out_max]=funmin_g(g,a,b,in_param);
fmax=-1*fmax;
out_max.f = f;

