
r = 3;  % order

for r=1:4
n = 256;
N = n;
x = (0:1/n:(n-1)/n);


tilde_g_0 = 0;
m = 1:(-1 + N/2);
tilde_g_h1 = N./abs(m).^r;
m = (N/2):(-1 + N);
tilde_g_h2 = N./abs(N-m).^r;
tilde_g = [tilde_g_0 tilde_g_h1 tilde_g_h2];
g = ifft(tilde_g);
plot(g)
hold on
end



const = -(-1)^(r/2)*((2*pi)^r)/factorial(r);
if r==2
  bernPoly = @(x)const*(-x.*(1-x) + 1/6);
else
  bernPoly = @(x)const*( ( (x.*(1-x)).^2 ) - 1/30);
end

figure; plot(bernPoly(x)); hold on


error = sum(abs(bernPoly(x) - g))

fprintf('')