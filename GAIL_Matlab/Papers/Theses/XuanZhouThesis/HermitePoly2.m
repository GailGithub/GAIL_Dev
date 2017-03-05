function hk2 = HermitePoly2(n)
H = HermitePoly(n);
hk2 = zeros(1,2*n-1);
H2 = H'*H;
for j = 1:n
    hk2(2*j-1) = hk(j)^2;
    for k = j:-1:1
        hk2(2*j-1) = hk2(j)+2*hk(k)*hk(2*j