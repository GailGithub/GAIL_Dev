function YX = Test_CV_integral_Two(n)
u = rand(n,2);
YX = [4.*exp(sum(u.^2,2)) sum(u.^2,2)];
end

