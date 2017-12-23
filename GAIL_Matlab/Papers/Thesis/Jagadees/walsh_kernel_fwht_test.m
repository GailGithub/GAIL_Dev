function walsh_kernel_fwht_test()

n = 2^3;
dim = 1;
sobstr=sobolset(dim); %generate a Sobol' sequence
x = sobstr(1:n,1:dim); %grab Sobol' points

xun = x; %bitrevorder(x);
order = 2;
a = 0.01;
[K, Ktilde] = kernel(xun,order,a)
fprintf('done')
end

% x : input points of size [n,d]
function dm = diffMatrix(x)
[n, d] = size(x);
A = reshape(x, n,1,d);
dm = repmat(A*n, [1, n, 1]);
A1 = reshape(x, 1,n,d);
dm1 = repmat(A1*n, [n, 1, 1]);
dm = bitxor(dm, dm1)/n;  % bitwise subtraction
end

function [K, Ktilde] = kernel(xun,order,a)

%a1 = @(x)(-floor(log2(x)));
    function out = a1(x)
        out = -floor(log2(x));
        out(x==0) = 0;  % a1 is zero when x is zero
    end

%t1 = @(x)(2.^(-a1(x)));
    function out = t1(x)
        out = (2.^(-a1(x)));
        out(x==0) = 0;  % t1 is zero when x is zero
    end
s1 = @(x)(1-2*x);
ts2 = @(x)((1-5*t1(x))/2 - (a1(x)-2).*x);

%omega2 = @(x, a)prod(1 + a*(s1(x) + ts2(x) - 1), ndims(x));
% to avoid subtracting "1", s1 is used directly
omega2_1D = @(x, a)(1 + a*(-2*x + ts2(x)));
[n, dim] = size(xun);
if dim > 1
    omega2 = @(x, a)prod(omega2_1D(x,a), ndims(x));
else
    omega2 = omega2_1D;
end
% figure(2); x=[0:0.001:1]; plot(x, omega2(x), '.'); grid on; axis([0 1 -1 2])

if a==1
    fprintf('debug')
end

K = omega2(xun, a);
Ktilde = abs(fwht(K, n, 'hadamard'));

if 1
    %figure() ; plot(K)
    dm = diffMatrix(xun);
    km = omega2(dm,a);
    lambda = eig(km)/length(xun);
    temp = lambda./Ktilde(n:-1:1);
    max(temp)
    figure; plot(temp)

    [V,I] = sort(diag(fwht(eye(n))*km*fwht(eye(n))), 'desc');
    [diag(hadamard(n)*km*hadamard(n)) V*n*n]
end

if ~isreal(Ktilde)
    fprintf('Ktilde has complex vals \n');
end

end

function g = bin2gray(b)
g(1) = b(1);
for i = 2 : length(b)
    x = xor(str2num(b(i-1)), str2num(b(i)));
    g(i) = num2str(x);
end
end
