function fval = fctVec(fct,dsites)
% Evaluate the function value vector
N = size(dsites,1);
fval = zeros(N,1);
for j = 1:N
    fval(j) = fct(dsites(j,:)');
end
end