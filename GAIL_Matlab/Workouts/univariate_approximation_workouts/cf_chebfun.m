function cf_chebfun(f, a, b, abstol)
% Examples:
% f1 = @(x) x.^4 .* sin(1./x); a = -1; b = 1; abstol = 1e-6; cf_chebfun(f1, a, b, abstol)
%
% f2 = @(x) f1(x) + 10.*x.^2; abstol = 1e-6;   cf_chebfun(f2, a, b, abstol) 
%
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
% f3 = @(x) B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%    - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; abstol = 1e-14;  
% cf_chebfun(f3, a, b, abstol)
%
% f4 = @(x)sin(10*pi*x.^4)-x, a = 1; b = 2; abstol = 1e-14; cf_chebfun(f4, a, b, abstol)
%
% f5 = @(x) sign(x);  a = -1; b = 1; cf_chebfun(f5, a, b, abstol)
%  

set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
format compact
format long

%% funappx_g
tic, [fappx, fout] = funappxNoPenalty_g(f,a,b,abstol,'nmax',10^8), toc
% gail.funappx_g_check(fappx,fout)
%% chebfun
tic, c = chebfun(f, [a,b],'splitting','on'), toc

x=a:0.00001:b;
figure(1)
subplot(2,3,1), plot(x,f(x)); title(['f(x) = ',func2str(f)]); axis tight
subplot(2,3,2), plot(x,fappx(x)); title(['funappx\_g approx.']); axis tight
subplot(2,3,3), plot(x,c(x)); title(['Chebfun approx.']); axis tight

err = abs( fappx(x) - f(x));
subplot(2,3,5), semilogy( x, err, 'k' );  title('funappx\_g errors'); axis tight; hold on
[~,ind] = find(err > abstol*10);
semilogy( x(ind), err(ind), 'ro' );   hold off;
figure(2);
semilogy( x, err, 'k' );  title('funappx\_g errors'); axis tight;  hold on;
semilogy( x(ind), err(ind), 'ro' );  
gail.save_eps('TraubPaperOutput', 'funappx_g_errors');
   
chebfuntol=1e-14;
err = abs(c(x) - f(x));
figure(1); subplot(2,3,6), semilogy( x, err, 'k' );   title ('Chebfun errors'); axis tight; hold on;
[~,ind] = find(err > chebfuntol*10);
semilogy( x(ind), err(ind), 'ro' );   hold off;
figure(3);
semilogy( x, err, 'k' );   title ('Chebfun errors'); axis tight;  hold on;
semilogy( x(ind), err(ind), 'ro' ); 
gail.save_eps('TraubPaperOutput', 'chebfun_errors');


