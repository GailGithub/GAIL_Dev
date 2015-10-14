function [fmin,out_min] = flat_bottom;
% Widen  output interval of flat bottom function
%
% Example:
% >> [fmin,out_min] = flat_bottom; 
%    intlen***
% >> x = out_min.intervals
%     x = 0.2311   0.7689
%
% 
format compact
f=@(x) exp(-1./(x-.5).^2);
[fmin,out_min] = funmin_g_CSC(f,0,1);
funmin_g_demo(fmin,out_min);
intervals = out_min.intervals; % funmin_g gives [ 0.2340 0.7659 ]
x = out_min.x;

y = f(x);
[~, index] = find(abs(y - fmin) < out_min.abstol);
leftint = find([1 diff(index)~=1]);
rightint = find([diff(index)~=1 1]);
[x(index(leftint)), x(index(rightint))]; 
% Shouldn't out_min.intervals be [0.2311   0.7689]
 