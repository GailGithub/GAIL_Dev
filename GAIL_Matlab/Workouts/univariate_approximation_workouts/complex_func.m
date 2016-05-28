% complex_func: funappx can approximate a complex-valued function defined on
% finite real valued interval [a,b]. 
%
% The example is taken from MATLAB documentation "doc INTERP1"
clear all; close all;
a = 1;
b = 10;
f = @(x) (5*x)+(x.^2*1i);
xq = a:0.001:b;
[fappx,out] = funappx_g(f,a,b,1e-13,100,100,10^8)
gail.funappx_g_check(fappx,out);
vq = fappx(xq);
v  = f(xq);
errest =  max(abs(v-vq))

figure
plot(xq,real(v),'*r',xq,real(vq),'-r');
hold on
plot(xq,imag(v),'*b',xq,imag(vq),'-b');
hold off