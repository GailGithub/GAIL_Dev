% funappxNoPenalty can approximate a complex-valued function defined on
% finite real valued interval [a,b]. 
%
% The example is taken from MATLAB documentation "doc INTERP1"

a = 1;
b = 10;
x = a:b;
f = @(x) (5*x)+(x.^2*1i);
xq = a:0.001:b;
[fappx,out] = funappxNoPenalty_g(f,a,b) 
errest =  max(abs(f(xq)-fappx(xq)))