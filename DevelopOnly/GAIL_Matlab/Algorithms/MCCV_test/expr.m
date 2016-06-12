function YX = expr(n) 
    p = rand(n, 1);
    Y = exp(p);
    YX = [Y p];
end