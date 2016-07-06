function YX = sinr_cv(n) 
    u = rand(n, 1);
    Y = sin(u);
    YX = [Y u];
end