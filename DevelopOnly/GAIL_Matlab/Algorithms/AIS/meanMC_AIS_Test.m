function tests = meanMC_AIS_Test
tests = functiontests(localfunctions);
end

% Defining the Keistertrue function to compare the results of the tests 
% that uses Keister integrals.

function I = Keistertrue(d)
%KEISTERTRUE computes the true value of the Keister integral in dimension d
%  accuracy might degrade as d increases due to round-off error
cosinteg=zeros(1,d);
cosinteg(1)=sqrt(pi)/(2*exp(1/4));
sininteg=zeros(1,d);
%sininteg(1)=integral(@(x) exp(-x.*x).*sin(x),0,inf);
sininteg(1)=4.244363835020225e-01;
cosinteg(2)=(1-sininteg(1))/2;
sininteg(2)=cosinteg(1)/2;
for j=3:d
   cosinteg(j)=((j-2)*cosinteg(j-2)-sininteg(j-1))/2;
   sininteg(j)=((j-2)*sininteg(j-2)+cosinteg(j-1))/2;
end
I=(2*(pi.^(d/2))/gamma(d/2))*cosinteg(d);
end

% Tests for the program meanMC_CLT_AIS:

function testMax_CLT(testCase)
Y = @(t,b_value)max(t+b_value,0).*exp(-t.*b_value-b_value.^2/2);
b = [-2 2];
d = 1;
abstol = 0.005;
alpha = 0.01;
nSig = 1e4;
fudge = 1.2;
actual = meanMC_CLT_AIS(Y, b, d, abstol, alpha, nSig, fudge);
expected = 1/sqrt(2*pi);
verifyEqual(testCase,actual,expected,'AbsTol',abstol)
end

function testKeister_CLT(testCase)
d = 1;
Y = @(x,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(x.*x,2))).*exp((1/2-b.^2).*sum(x.*x,2));
b = [0.5 3];
abstol = 0.005;
alpha = 0.01;
nSig = 1e4;
fudge = 1.2;
actual = meanMC_CLT_AIS(Y, b, d, abstol, alpha, nSig, fudge);
expected = Keistertrue(1);
verifyEqual(testCase,actual,expected,'AbsTol',abstol)
end

% Tests for the program meanMC_AIS_g:

function testMax_g(testCase)
Y = @(t,b_value)max(t+b_value,0).*exp(-t.*b_value-b_value.^2/2);
b = [-2 2];
d = 1;
abstol = 0.005;
alpha = 0.01;
nSig = 1e4;
fudge = 1.2;
actual = meanMC_AIS_g(Y, b, d, abstol, alpha, nSig, fudge);
expected = 1/sqrt(2*pi);
verifyEqual(testCase,actual,expected,'AbsTol',abstol)
end

function testKeisterDimen_g(testCase)
n = 3;
b = [0.5 3];
abstol = 0.005;
alpha = 0.01;
nSig = 1e4;
fudge = 1.2;
for d = 1:n
    Y = @(x,b) ((sqrt(2.*pi).*b).^d).*cos(b.*sqrt(sum(x.*x,2))).*exp((1/2-b.^2).*sum(x.*x,2));
    actual = meanMC_AIS_g(Y, b, d, abstol, alpha, nSig, fudge);
    expected = Keistertrue(d);
  verifyEqual(testCase,actual,expected,'AbsTol',abstol)      
end
end
