function YY = expr_av(n) 
    u = rand(n, 1);
    u_a = 1 - u;
    Y = exp(u);
    Y_a = exp(u_a);
    YY = [Y Y_a];
end