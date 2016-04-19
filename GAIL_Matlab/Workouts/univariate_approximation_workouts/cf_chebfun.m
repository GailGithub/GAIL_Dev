function cf_chebfun(f, a, b)
% Examples:
% f1 = @(x) x.^4 .* sin(1./x); a = -1; b = 1; cf_chebfun(f1, a, b)
%
% f2 = @(x) f1(x) + 10.*x.^2;  cf_chebfun(f2, a, b)
%
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
% f3 = @(x) B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%    - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; cf_chebfun(f3, a, b)
%
% f = @(x) sign(x);  a = -1; b = 1; cf_chebfun(f, a, b)
% 

tic, [fappx, fout] = funappxNoPenalty_g(f,a,b,1e-13,5,5,10^8); toc
%tic, [fappx2, fout] = funappx_g(f,a,b,1e-15,'nmax',10^8), toc

splitting on;
tic, c = chebfun(f, [a,b]); toc

x=a:0.00001:b;
subplot(2,3,1), plot(x,f(x)); title(['f(x) = ',func2str(f)]); axis tight
subplot(2,3,2), plot(x,fappx(x)); title(['funappx\_g approx.']); axis tight
subplot(2,3,3), plot(x,c(x)); title(['Chebfun approx.']); axis tight

subplot(2,3,5), semilogy( x, abs( fappx(x) - f(x)), 'k' );  title('funappx\_g error'); axis tight
err = abs(c(x) - f(x));
subplot(2,3,6), semilogy( x, err, 'k' );   title ('chebfun error'); axis tight; hold on;
[~,ind] = find(err > 1e-14);
semilogy( x(ind), err(ind), 'ro' );   hold off;
   
keyboard
