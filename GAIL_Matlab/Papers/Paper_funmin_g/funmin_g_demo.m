function funmin_g_demo(fmin, out_param)
% 
% Example 1: f is a flat bottom function. Cf. Wolfram Alpha's result at 
% tinyurl.com/ox63f2c
%
% >> f = @(x) 1./sqrt(1-x.^4);  a = -0.2; b = 0.2; tol = 1e-15;
% >> [fmin, out] = funmin_g(f,a,b,tol,tol);
% >> abs(fmin-1) < 1e-6
%    ans = 1 
% >> funmin_g_demo(fmin, out)
%    intlen = ***e-04 fval1 = ***e-16 fval2 = ***e-16
%

  
a = out_param.a;
b = out_param.b;
f = out_param.f;
h = 1.0/out_param.nmax;

h = 0.00001;
x = a:h:b;
plot(x,f(x));

if ~isempty(out_param.intervals)
    hold on;
    if out_param.volumeX > out_param.TolX,
        intlen = out_param.intervals(2) - out_param.intervals(1)
        fval1 = f(out_param.intervals(2)) - fmin
        fval2 = f(out_param.intervals(1)) - fmin
        plot(out_param.intervals(1), fmin, 'r<', 'MarkerSize', 12,'LineWidth',2)
        plot(out_param.intervals(2), fmin, 'r>', 'MarkerSize', 12,'LineWidth',2)
    end
    plot(mean(out_param.intervals),  fmin, 'ro', 'MarkerSize', 12,'LineWidth',2)
    hold off;
end



