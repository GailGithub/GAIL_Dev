function wingo
% Examples 1 and 2 from: Wingo, Dallas R. "Globally minimizing polynomials
% without evaluating derivatives." International Journal of Computer
% Mathematics 17.3-4 (1985): 287-294.
%
% >> wingo
%    check = 1 
%    check = 1

%% Example 1
f = @(x) 1.0/6 * x.^6 - 52.0 / 25 * x.^5 + 39.0 / 80 * x.^4 ...
       + 71.0/10 * x.^3 - 79.0/20 * x.^2 - x + 1.0/10;
a = -1.5; b = 11; tol = 1e-7; Xtol = 1e-4;
[fmin,out_param] = funmin_g(f,a,b,tol,Xtol,10,10,1000000);
funmin_g_demo(fmin, out_param);

fmin_true = -29763.233; xmin_true = 10;
ferror = abs(fmin - fmin_true);
xerror = abs(mean(out_param.intervals) - xmin_true);
check = ferror < tol || xerror < Xtol

%% Example 2
f = @(x) (x-1).^10;
a = 0; b = 2; tol = 1e-7; Xtol = 1e-4;
[fmin,out_param] = funmin_g(f,a,b,tol,Xtol,100,100,1000000);
funmin_g_demo(fmin, out_param);

fmin_true = 0; xmin_true = 1;
ferror = abs(fmin - fmin_true);
xerror = abs(mean(out_param.intervals) - xmin_true);
check = ferror < tol || xerror < Xtol

 
