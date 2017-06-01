function cf_chebfun_min(f, a, b, abstol, truey)
% CF_CHEBFUN_MIN Compares funmin_g with Chebfun
%
% Example 1:
% f1 = @(x) x.^4 .* sin(1./((x==0)+x)); a = -1; b = 1; abstol = 1e-6; cf_chebfun_min(f1, a, b, abstol,f1(-1))
%
% Example 2:
% f2 = @(x) f1(x) + 10.*x.^2;  cf_chebfun_min(f2, a, b, abstol, 0) 
%
% Example 3:
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
% f3 = @(x) -B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%    - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; abstol = 1e-6;  
% cf_chebfun_min(f3, a, b, abstol, -1)
%
% Example 4:
% f4 = @(x)sin(10*pi*x.^4)-x, a = 0; b = 2; abstol = 1e-6; truey= -3; cf_chebfun_min(f4, a, b, abstol, truey)
%
% Example 5:
% f5 = @(x) sign(x);  a = -1; b = 1; abstol = 1e-6; truey= -1; cf_chebfun_min(f5, a, b, abstol, truey)
%  

set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
format compact
format long

%% funmin_g
tic, [fmin,out] = funmin_g(f, a, b, abstol), toc
%out.intervals
funmin_g_demo(fmin, out);
trueerr = abs(truey - fmin)
hold on; 

%% fminbnd
tic, [fx,fy] = fminbnd(@(x) f(x), a, b), toc
plot(fx,fy,'bs')
trueerr=abs(truey - fy)

%% chebfun
tic, c = chebfun(f, [a,b],'chebfuneps', abstol,'splitting','on');
tic, [cy,cx] = min(c), toc
trueerr = abs(truey - cy)
plot(cx,cy,'g*')
legend('f', 'funmin\_g', 'fminbnd', 'chebfun')
hold off

%% symbolic toolbox
% syms x
% g = sym(f);
% assume(a <= x <= b) 
% sym_x_min = solve(diff(g,'x'))
% sym_y_min = subs(g,x,sym_x_min)
% syms x clear

%% integral_g
%tic, [fint,out3] = integralNoPenalty_g(f, a, b, abstol,'nmax',10^8), toc
%tic, cint = sum(c),toc
%keyboard
