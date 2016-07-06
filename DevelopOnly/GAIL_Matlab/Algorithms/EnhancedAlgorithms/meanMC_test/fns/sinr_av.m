function YY = sinr_av(n) 
    u = rand(n, 1);
    u_a = 1 - u;
    Y = sin(u);
    Y_a = sin(u_a);
    YY = [Y Y_a];
end