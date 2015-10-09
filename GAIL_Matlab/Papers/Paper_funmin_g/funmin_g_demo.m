function [fmin, out] = funmin_g_demo(f,a,b)
% 
% Example 1: f is a flat bottom function. Cf. Wolfram Alpha's result at 
% tinyurl.com/ox63f2c
%
% >> f = @(x) 1./sqrt(1-x.^4);  a = -0.5; b = 0.5; 
% >> [fmin, out] = funmin_g_demo(f,a,b);
%    ***
% >> abs(fmin-1) < 1e-6  
%    ans = 1 
%
%
% 
  
[fmin, out] = funmin_g(f,a,b);


h = 0.00001;
x = a:h:b;
plot(x,f(x));

if ~isempty(out.intervals)
    hold on;
    plot(mean(out.intervals), fmin, 'bo')
    hold off;
end



