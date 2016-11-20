function YX = Test_CV_integral(n)
u = rand(n,1);
YX = [2.*exp(-4.*u.^2) exp(-2.*u)];
end

