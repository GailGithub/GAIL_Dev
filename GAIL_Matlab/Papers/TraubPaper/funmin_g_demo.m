function [intlen,fval1,fval2] = funmin_g_demo(fmin, out_param)
% FUNMIN_G_DEMO visualize f and minimum from funmin_g
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
LatexInterpreter;
a = out_param.a;
b = out_param.b;
f = out_param.f;
h = 1.0/out_param.nmax;

%h = 0.00001;
x = a:h:b;
figure;
plot(x,f(x));
xlabel('$x$')

hold on;
%plot(out_param.output_x, fmin, 'ro', 'MarkerSize', 12,'LineWidth',2)
[m,n]=size(out_param.intervals);
if ~isempty(out_param.intervals)
   if (n>1)
     plot([mean(out_param.intervals(:,1)), mean(out_param.intervals(:,n))], [fmin, fmin], 'ro', 'MarkerSize', 12,'LineWidth',2)
   else
    plot(mean(out_param.intervals(:,1)),  fmin, 'ro', 'MarkerSize', 12,'LineWidth',2)
   end
end
hold off;

legend('$f$', 'minimum from \texttt{funmin\_g}')


