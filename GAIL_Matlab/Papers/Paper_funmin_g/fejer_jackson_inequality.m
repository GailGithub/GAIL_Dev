function [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(m)
% The Fejer-Jackson inequality says that the following partial sum is 
% positive for a given positive integer n and for x in (0, pi): 
%     fnx = @(n,x)((1./(1:n)) * sin((1:n)' * x));
%
% Example 1:
% >> [fmin, fmax, out_min, out_max] =  fejer_jackson_inequality(32);
% >> abs(fmin) < 1e-6, abs(fmax-1.804097) < 1e-6
%  ans = 1  ans = 1
% >>  funmin_g_demo(fmin, out_min); hold on; funmin_g_demo(fmax, out_max); hold off;
%    ***intlen = 0.00***
%
%
% Source: Fejer-Jackson inequality, Nick Trefethen
%         http://www.chebfun.org/examples/fourier/FejerJackson.html


fnx = @(n,x)((1./(1:n)) * sin((1:n)' * x));
f = @(x) fnx(m,x);
[fmin,out_min]=funmin_g(f,0,pi); % Why is 0 not in out_min.intervals?

g = @(x) fnx(m,x)*-1;
[fmax,out_max]=funmin_g(g,0,pi);
fmax=-1*fmax;
out_max.f = f;



