function YX = distfun_cv(n) 
    p1 = rand(n,2);
    p2 = rand(n,2);
    Y = sqrt(sum((p1 - p2).^2,2));
    YX = [Y p1 p2];
end