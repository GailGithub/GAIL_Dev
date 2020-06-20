%GENZ_TEST_FUN_TRUE  Provides true integrals of seven test functions
%from Alan Genz and Keister. The test function accepts six inputs and
%give the output of the true solution of the integral.
%
% The first argument is hyperbox, which is a matrix of size 2 x dim.
%
% The second argument is index, which identifies a test function to use.
% the value is 1,...,7.
%
% The third argument is the dimension of the integrand.
%
% The fourth argument is the parameter alpha, which should be a row vector of
% size 1 x dim.
%
% The fifth argument is the parameter beta, which should be a row vector of
% size 1 x dim.
%
% The sixth argument is the paramter r, which sould be a constant.
%
function f_true = genz_test_fun_true(hyperbox,index,dim,alpha,beta,r)
switch index
    case 1 % Genz "Oscillatory"
        s = zeros(2^dim,1);
        sign = zeros(2^dim,1);
        s(1) = 2*pi*r+hyperbox(2,:)*alpha(1:dim)';
        sign(1) = 1;
        for i = 1:dim
            s(2^(i-1)+1:2^i) = s(1:2^(i-1))-alpha(dim-i+1)*(hyperbox(2,dim-i+1)-hyperbox(1,dim-i+1));
            sign(2^(i-1)+1:2^i) = -sign(1:2^(i-1));
        end
        switch mod(dim,4)
            case 1
                f_true = sum(sign.*sin(s))/prod(alpha(1:dim));
            case 2
                f_true = sum(sign.*(-cos(s)))/prod(alpha(1:dim));
            case 3
                f_true = sum(sign.*(-sin(s)))/prod(alpha(1:dim));
            case 0
                f_true = sum(sign.*cos(s))/prod(alpha(1:dim));
        end
        
    case 2 % Genz "Product Peak";
        f_true = prod(atan((hyperbox(2,:)-beta(1:dim))./alpha(1:dim))./alpha(1:dim)...
            - atan(hyperbox(1,:)-beta(1:dim)./alpha(1:dim))./alpha(1:dim));
        
    case 3 % Genz "Corner Peak"
        s = zeros(2^dim,1);
        sign = zeros(2^dim,1);
        s(1) = 1+hyperbox(2,:)*alpha(1:dim)';
        sign(1) = 1;
        for i = 1:dim
            s(2^(i-1)+1:2^i) = s(1:2^(i-1))-alpha(dim-i+1)*(hyperbox(2,dim-i+1)-hyperbox(1,dim-i+1));
            sign(2^(i-1)+1:2^i) = -sign(1:2^(i-1));
        end
        f_true = (-1)^dim*sum(sign.*s.^(-r))/prod(alpha(1:dim))/prod(r:r+dim-1);
        
    case 4 % Genz "Gaussian"
        f_true = (sqrt(pi)/2)^dim*prod(1./alpha(1:dim))*...
            prod(erf(hyperbox(2,:).*alpha(1:dim)-alpha(1:dim).*beta(1:dim))-...
            erf(hyperbox(1,:).*alpha(1:dim)-alpha(1:dim).*beta(1:dim)));
        
    case 5 % Genz "Continuous"
        f_true_iter = ones(dim,1);
        for i =1:dim
            if beta(i) <= hyperbox(1,i)
                f_true_iter(i) = -1/alpha(i)*(exp(alpha(i)*(beta(i)-hyperbox(2,i)))...
                    -exp(alpha(i)*(beta(i)-hyperbox(1,i))));
            elseif beta(i) >= hyperbox(2,i)
                f_true_iter(i) = 1/alpha(i)*(exp(alpha(i)*(hyperbox(2,i)-beta(i)))...
                    -exp(alpha(i)*(hyperbox(1,i)-beta(i))));
            else
                f_true_iter(i) = 1/alpha(i)*(2-exp(alpha(i)*...
                    (beta(i)-hyperbox(2,i)))-exp(alpha(i)*...
                    (hyperbox(1,i)-beta(i))));
            end
        end
        f_true = prod(f_true_iter);
        
    case 6 % Genz "Discontinuous"
        if ( any ( beta(1:dim) < hyperbox(1,:) ))
            f_true = 0;
        else
            f_true = prod((exp(alpha(1:dim).*min(hyperbox(2,1:dim),beta(1:dim)))-...
                exp(alpha(1:dim).*hyperbox(1,1:dim)))./alpha(1:dim));
        end
        
    case 7 % Keister test function
        if hyperbox == [zeros(1,dim);ones(1,dim)];
            if dim == 1
                f_true = 1.3803884470431429744;
            end
            if dim == 2
                f_true = 1.808186634594926;
            end
            if dim == 3
                f_true = pi^(3/2)/(2*exp(1/4));
            end
            if dim == 4
                f_true = 2.1659293025745063474;
            end
            if dim == 5
                f_true = pi^(5/2)/(12*exp(1/4));
            end
            if dim == 6
                f_true = -2.3273037292979391292;
            end
            if dim == 7
                f_true = -31*pi^(7/2)/(120*exp(1/4));
            end
            if dim == 8
                f_true = -30.609075003558562675;
            end
        end
end
f_true = f_true ./ prod(hyperbox(2,:)-hyperbox(1,:)); % pdf of uniform
end
