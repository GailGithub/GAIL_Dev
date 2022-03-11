f = @(x) -5*exp(-100*(x-0.15).^2) - exp(-80*(x-0.65).^2);
[fmin,out] = funmin_g(f);