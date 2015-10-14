function [fmin, out_param] = flatbottom2(a,b);
%
% a = -1; b = 1;   %singularities
%
% Example 1:
% >> [fmin, out_param] = flatbottom2(-.5,0.5);
% >> abs(fmin-1) < out_param.abstol
%    ans = 1
%
%
% Example 2:
% >> [fmin, out_param] = flatbottom2(-1,1);
%      *** Function f(x) = Inf at x = -1 1


f = @(x) 1./sqrt(1-x.^4);  

[fmin, out_param] = funmin_g_CSC(f,a,b);

funmin_g_demo(fmin, out_param);
 