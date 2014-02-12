%GENZ_TEST_FUN_TRUE  Provides true integrals of seven test functions
%from Alan Genz and Keister. The test function accepts six inputs and
%give the output of the true solution of the integral.
%
% The first argument is hyperbox, which is a matrix of size 2 x dim.
% The second argument is index, which identifies a test function to use.
% The third argument is the dimension of the integrand.
% The fourth argument is the parameter alpha, which is a dim x 1 vector.
% The fifth argument is the parameter beta.
% The sixth argument is the parameter r.
function f_true = genz_test_fun_true(hyperbox,index,dim,alpha,beta,r)
switch index
    case 1 % Genz "Oscillatory"
        if dim == 1;
            f_true = sin(2*pi*r + alpha(1)*hyperbox(2))/alpha(1) - ...
                sin(2*pi*r + alpha(1)*hyperbox(1))/alpha(1);
        elseif dim == 2
            f_true = ((-cos(2*pi*r + alpha(1)* hyperbox(2,1)+alpha(2)*hyperbox(2,2))+...
                cos(2*pi*r + alpha(1)* hyperbox(1,1)+alpha(2)*hyperbox(2,2)))...
                -(-cos(2*pi*r + alpha(1)* hyperbox(2,1)+alpha(2)*hyperbox(1,2))+...
                cos(2*pi*r + alpha(1)* hyperbox(1,1)+alpha(2)*hyperbox(1,2))))...
                /(alpha(1)*alpha(2));
        else
            s = zeros(dim,1);
            sign = zeros(dim,1);
            s(1) = 2*pi*r+hyperbox(2,:)*alpha;
            for i = 1:dim
                s(2^(i-1)+1:2^i) = s(1:2^(i-1))-alpha(dim-i+1)*(hyperbox(2,dim-i+1)-hyperbox(1,dim-i+1));
            end
            sign(1) = 1;
            for i = 1:dim
                sign(2^(i-1)+1:2^i) = -sign(1:2^(i-1));
            end
            switch mod(dim,4)
                case 1
                    f_true = sum(sign.*sin(s))/prod(alpha);
                case 2
                    f_true = sum(sign.*(-cos(s)))/prod(alpha);
                case 3
                    f_true = sum(sign.*(-sin(s)))/prod(alpha);
                case 0
                    f_true = sum(sign.*cos(s))/prod(alpha);
            end
            %f_true = nan;
        end
        
    case 2 % Genz "Product Peak";
        f_true = prod(atan((hyperbox(2,:)-beta(1:dim))./alpha(1:dim))./alpha(1:dim)...
            - atan(hyperbox(1,:)-beta(1:dim)./alpha(1:dim))./alpha(1:dim));
        
    case 3 % Genz "Corner Peak"
        if dim == 1
            f_true = -1/(alpha(1)*r*(1+alpha(1)*hyperbox(2))^r)...
                +1/(alpha(1)*r*(1+alpha(1)*hyperbox(1))^r);
        elseif dim == 2
            f_true = 1/(r*(r+1)*alpha(1)*alpha(2))...
                *(((1+alpha(1)*hyperbox(2,1)+alpha(2)*hyperbox(2,2))^(-r)...
                -(1+alpha(1)*hyperbox(1,1)+alpha(2)*hyperbox(2,2))^(-r))-...
                ((1+alpha(1)*hyperbox(2,1)+alpha(2)*hyperbox(1,2))^(-r)...
                -(1+alpha(1)*hyperbox(1,1)+alpha(2)*hyperbox(1,2))^(-r)));
        else
            s = zeros(dim,1);
            sign = zeros(dim,1);
            s(1) = 1+hyperbox(2,:)*alpha;
            for i = 1:dim
                s(2^(i-1)+1:2^i) = s(1:2^(i-1))-alpha(dim-i+1)*(hyperbox(2,dim-i+1)-hyperbox(1,dim-i+1));
            end
            sign(1) = 1;
            for i = 1:dim
                sign(2^(i-1)+1:2^i) = -sign(1:2^(i-1));
            end
            f_true = (-1)^dim*sum(sign.*s.^(-r))/prod(alpha)/prod(r:r+dim-1);
            %f_true = nan;
        end
        
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
            f_true_iter = ones(dim,1);
            for i = 1:dim
                f_true_iter(i) = (exp(alpha(i).*min(hyperbox(2,i),beta(i)))-...
                    exp(alpha(i).*hyperbox(1,i)))./alpha(i);
            end
            f_true = prod(f_true_iter);
        end
        
    case 7 % Keister test function
        if hyperbox == [zeros(1,dim);ones(1,dim)];
            if dim == 1
                f_true = 1.2238;
            end
            if dim == 2
                f_true = 0.9040;
            end
            if dim == 3
                f_true = 0.6113;
            end
        else
            f_true =  nan;
        end
end
end
