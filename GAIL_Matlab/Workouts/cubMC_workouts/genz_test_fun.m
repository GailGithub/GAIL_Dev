%This is a suite of test functions from Allan Genz, which contains six
%multidimensional integrands. the test function accepts six inputs and give
%the output of the function values.
% the first API is x, which is the random variable, x should be a matrix
% of size n x d. n is the sample size, d the dimension.
% the second API is index, which identify which test function to use.
% the third API is the dimension of the integrand.
% the fourth API is the parameter alpha.
% the fifth API is the parameter beta.
% the sixth API is the paramter r.
function f_val = genz_test_fun (x,index,dim,alpha,beta,r)
if index == 1 % Genz "Oscillatory"
    f_val = cos (2.0 * pi * r + sum (x(:,1:dim).*...
        repmat(alpha(1:dim),size(x,1),1),2));
end
if index == 2; % Genz "Product Peak";
    f_val = 1/prod(repmat(alpha(1:dim),size(x,1),1).^2+...
        (x(:,1:dim)-repmat(beta(1:dim),size(x,1),1)).^2,2);
end
if index == 3 % Genz "Corner Peak"
    f_val = 1/(1+sum(repmat(alpha(1:dim),size(x,1),1).*x(:,1:dim),2))...
        .^(dim+r);
end
 if index == 4 % Genz "Gaussian"
     f_val = exp(-sum(repmat(alpha(1:dim),size(x,1),1).^2.*...
         (x(:,1:dim)-repmat(beta(1:dim),size(x,1),1)).^2,2));
 end
 if index == 5 % Genz "Continuous"
     f_val = exp(-sum(repmat(alpha(1:dim),size(x,1),1).*...
         abs(x(:,1:dim)-repmat(beta(1:dim),size(x,1),1)),2));
 end
 if index == 6 % Genz "Discontinuous"
     if ( any ( beta(1:dim) < x(1:dim) ))
         f_val = 0;
     else
         f_val = exp(sum(repmat(alpha(1:dim),size(x,1),1).*x(:,1:dim),2));
     end
 end
end