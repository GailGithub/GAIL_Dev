function [fmin, fmax, out_min, out_max] = fejer_jackson_inequality(m)
%  
% Example 1:
% >> [fmin, fmax] =  fejer_jackson_inequality(32);
% >> abs(fmin) < 1e-6, abs(fmax-1.804097) < 1e-6
%  ans = 1  ans = 1
%
%
%
% Source: Fejer-Jackson inequality, Nick Trefethen
%         http://www.chebfun.org/examples/fourier/FejerJackson.html


fnx = @(n,x)(sin(x'*(1:n))*(1./(1:n))')';
f = @(x) fnx(m,x);
[fmin,out_min]=funmin_g(f,0,pi);

g = @(x) fnx(m,x)*-1;
[fmax,out_max]=funmin_g(g,0,pi);
fmax=-1*fmax;

h = 0.000001;
x = 0:h:pi;
plot(x,f(x))