syms gamma x t;
integrand = @(x,t) exp(-gamma^2*(x-t)^2-x^2-t^2);
result = int(int(integrand(x,t),x,-inf,inf),t,-inf,inf);

