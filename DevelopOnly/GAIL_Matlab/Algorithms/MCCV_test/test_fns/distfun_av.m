function YX = distfun_av(n) 
    p1 = rand(n,2);
    p2 = rand(n,2);

    Y = sqrt(sum((p1 - p2).^2,2));
    Y_a = sqrt(sum(((1-p1) - (1-p2)).^2,2));
    YX = [Y Y_a];
end