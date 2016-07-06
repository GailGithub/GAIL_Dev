function YX = expr_cv(n) 
    u = rand(n, 1);
    Y = exp(u);
    YX = [Y u];
end