wKernel = @(x)12*( (1/4) - 2.^(floor(log2(x))-1) );
%figure(1); x=[0:0.001:1]; plot(x, wKernel(x), '.'); grid on; axis([0 1 -1 2])

a1 = @(x)(-floor(log2(x )));
t1 = @(x)(2.^(-a1(x)));
t2 = @(x)(2.^(-2*a1(x)));
s1 = @(x)(1-2*x);
s2 = @(x)(1/3 - 2*(1-x).*x);
ts2 = @(x)((1-5*t1(x))/2 - (a1(x)-2).*x);
ts3 = @(x)((1-43*t2(x))/18 +(5*t1(x)-1).*x + (a1(x)-2).*(x.^2));

omega2 = @(x)(s1(x) + ts2(x));
figure(2); x=[0:0.001:1]; plot(x, omega2(x), '.'); grid on; axis([0 1 -1 2])
figure(21); x=[0:0.001:1]; plot(x, omega2(0.2*x), '.'); grid on; axis([0 1 -1 2])
figure(22); x=[0:0.001:1]; plot(x, omega2(2*x), '.'); grid on; axis([0 1 -1 2])

omega3 = @(x)(s1(x) + s2(x) + ts3(x));
%figure(3); x=[0:0.001:1]; plot(x, omega3(x), '.'); grid on; axis([0 1 -1 2])

% based on the Omega plots
% reordered with bitrevorder
M1 = [2    1/2  0   -1/2; ...
      1/2  2   -1/2  0; ...
      0   -1/2  2    1/2; ...
     -1/2  0    1/2  2];

diffMat = [0 1/2 1/4 3/4; 1/2 0 3/4 1/4; 1/4 3/4 0 1/2; 3/4 1/4 1/2 0];
h4 = hadamard(4);
h42 = h4*h4
M2 = wKernel(diffMat)

trans = h4*M2*h4
[e,v] = eig(M2,'vector')

return


[e,v] = eig(M1);

h4'*M1

h4*M1

h4