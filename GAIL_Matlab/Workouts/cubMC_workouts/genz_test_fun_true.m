%GENZ_TEST_FUN_TRUE Provides true integrals of seven test functions
% from Allan Genz and Keister. the test function accepts six
%inputs and give the output of the integration.
% the first API is interval, which is a matrix of size 2 x dim.
% the second API is index, which identify which test function to use.
% the third API is the dimension of the integrand.
% the fourth API is the parameter alpha.
% the fifth API is the parameter beta.
function f_true = genz_test_fun_true (interval,index,dim,alpha,beta)
switch index
    case 2; % Genz "Product Peak";
        f_true = prod(atan((interval(2,:)-beta(1:dim))./alpha(1:dim))./alpha(1:dim)...
            - atan(interval(1,:)-beta(1:dim)./alpha(1:dim))./alpha(1:dim));
    case 4 % Genz "Gaussian"
        f_true = (sqrt(pi)/2)^dim*prod(1./alpha(1:dim))*...
            prod(erf(interval(2,:).*alpha(1:dim)-alpha(1:dim).*beta(1:dim))-...
            erf(interval(1,:).*alpha(1:dim)-alpha(1:dim).*beta(1:dim)));
    case 5 % Genz "Continuous"
        f_true_iter = ones(dim,1);
        for i =1:dim
            if beta(i) <= interval(1,i)
                f_true_iter(i) = -1/alpha(i)*(exp(alpha(i)*(beta(i)-interval(2,i)))...
                    -exp(alpha(i)*(beta(i)-interval(1,i))));
            elseif beta(i) >= interval(2,i)
                f_true_iter(i) = 1/alpha(i)*(exp(alpha(i)*(interval(2,i)-beta(i)))...
                    -exp(alpha(i)*(interval(1,i)-beta(i))));
            else
                f_true_iter(i) = 1/alpha(i)*(2-exp(alpha(i)*...
                    (beta(i)-interval(2,i)))-exp(alpha(i)*...
                    (interval(1,i)-beta(i))));
            end
        end
        f_true = prod(f_true_iter);
    case 6 % Genz "Discontinuous"
        if ( any ( beta(1:dim) < interval(1,:) ))
            f_true = 0;
        else
            f_true_iter = ones(dim,1);
            for i = 1:dim        
            f_true_iter(i) = (exp(alpha(i).*min(interval(2,i),beta(i)))-...
                exp(alpha(i).*interval(1,i)))./alpha(i);
            end
            f_true = prod(f_true_iter);
        end
    case 7 % Keister test function
        if dim == 1
            f_true = 1.2;
        end
        if dim == 2
            f_true = 0.9;
        end
        if dim == 3
            f_true = 0.6;
        end
        if dim == 4
            f_true = 0.3;
        end
end
end
