function cf_chebfun(f, a, b)
% Examples:
% f1 = @(x) x.^4 .* sin( 3.371568833721756 ./x); a = -1; b = 1; cf_chebfun(f1, a, b)
%
% f2 = @(x) f1(x) + 10.*x.^2;  cf_chebfun(f2, a, b)
%
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
% f3 = @(x) B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%    - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; cf_chebfun(f3, a, b)
%
% f = @(x) sign(x);  a = -1; b = 1; cf_chebfun(f, a, b)
% 
format compact
format long
abstol = 1e-8;

%% funappx_g
tic, [fappx, fout] = funappxNoPenalty_g(f,a,b,abstol,'nmax',10^8), toc
gail.funappx_g_check(fappx,fout)
%% chebfun
splitting on;
tic, c = chebfun(f, [a,b]), toc

x=a:0.00001:b;
subplot(2,3,1), plot(x,f(x)); title(['f(x) = ',func2str(f)]); axis tight
subplot(2,3,2), plot(x,fappx(x)); title(['funappx\_g approx.']); axis tight
subplot(2,3,3), plot(x,c(x)); title(['Chebfun approx.']); axis tight

err =abs( fappx(x) - f(x));
subplot(2,3,5), semilogy( x, err, 'k' );  title('funappx\_g error'); axis tight; hold on
[~,ind] = find(err > abstol);
semilogy( x(ind), err(ind), 'ro' );   hold off;
   
err = abs(c(x) - f(x));
subplot(2,3,6), semilogy( x, err, 'k' );   title ('chebfun error'); axis tight; hold on;
[~,ind] = find(err > abstol);
semilogy( x(ind), err(ind), 'ro' );   hold off;
   
%% funmin_g
tic, [fmin,out2] = funminNoPenalty_g(f, a, b, abstol, abstol), toc
figure;
funmin_g_demo(fmin, out2);

%% fminbnd
tic, x = fminbnd(@(x) f(x),a,b), toc

%% chebfun
tic, [cy,cx] = min(c), toc

%% integral_g
tic, [fint,out3] = integralNoPenalty_g(f, a, b, abstol), toc
tic, cint = sum(c),toc
keyboard
