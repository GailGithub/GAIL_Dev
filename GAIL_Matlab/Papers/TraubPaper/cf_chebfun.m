function cf_chebfun(f, a, b, abstol, varargin)
% CF_CHEBFUN compares funappx_g with Chebfun
%
% Example 1:
% f1 = @(x) x.^4 .* sin(1./((x==0)+x)); a = -1; b = 1; abstol = 1e-6; cf_chebfun(f1, a, b, abstol)
%
% Example 2:
% f2 = @(x) f1(x) + 10.*x.^2; abstol = 1e-6;   cf_chebfun(f2, a, b, abstol) 
%
% Example 3:
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
% f3 = @(x) B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%    - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; abstol = 1e-14;  
% cf_chebfun(f3, a, b, abstol)
%
% Example 4:
% f4 = @(x)sin(10*pi*x.^4)-x, a = 1; b = 2; abstol = 1e-14; cf_chebfun(f4, a, b, abstol)
%
% Example 5:
% f5 = @(x) sign(x);  a = -1; b = 1; cf_chebfun(f5, a, b, abstol)
%  
%  

gail.InitializeDisplay
%set(0,'defaultaxesfontsize',24,'defaulttextfontsize',24) %make font larger
% format compact
format long

% MATLABBlue = [0, 0.447, 0.741];
% MATLABOrange = [0.85,  0.325, 0.098];
% MATLABPurple = [0.494,  0.184, 0.556];
% MATLABGreen = [0.466,  0.674, 0.188];
% MATLABDkOrange = [0.85,  0.325, 0.098]*0.6;
% MATLABLtOrange = 0.5*[0.85,  0.325, 0.098] + 0.5*[1 1 1];
% if nargin < 5,
%     testparallel = false;
% else
%     testparallel = varargin{1};
% end

%% funappx_g
t1 = 0;
%if ~testparallel,
    tic, [fappx, fout] = funappx_g(f,a,b,abstol,'nmax',10^8), t1=toc
    disp('---------------------');
    % gail.funappx_g_check(fappx,fout)
% else
%     tic, [fappx, fout] = par_funappx_g(4, f,a,b,abstol,'nmax',10^8), t2=toc
%     time_ratio1 = t1/t2
%     disp('---------------------');
% end

%% chebfun
tic, c = chebfun(f,[a,b],'chebfuneps', abstol,'splitting','on'), t3=toc
time_ratio2 = t1/t3
disp('---------------------');

x=a:0.00001:b;
figure(1)
subplot(2,3,1), plot(x,f(x)); 
%title(fstr); axis tight
subplot(2,3,2), plot(x,fappx(x)); title(['funappx\_g approx.']); axis tight
subplot(2,3,3), plot(x,c(x)); title(['Chebfun approx.']); axis tight

err1 = abs( fappx(x) - f(x));
subplot(2,3,5), semilogy( x, err1, 'k' );  title('funappx\_g errors'); axis tight; hold on
[~,ind1] = find(err1 > abstol*10);
semilogy( x(ind1), err1(ind1), '.' );   hold off;

chebfuntol=1e-14;
err = abs(c(x) - f(x));
figure(1); subplot(2,3,6), semilogy( x, err, 'k' );   title ('Chebfun errors'); axis tight; hold on;

[~,ind] = find(err > chebfuntol*10);
semilogy( x(ind), err(ind), 'r.' );   hold off;


figure(2);
h=semilogy( x, err1, '-', x(ind1), err1(ind1), '.', 'color', MATLABOrange);   hold on;
% semilogy( x, err, '-', 'color', MATLABBlue); hold on
% semilogy( x(ind), err(ind), '.' , 'color', MATLABOrange);  
small1 = max(-20,log10(0.1*min(err1)));
large1 = log10(10*max(err1));
axis([a b 10^small1 10^large1])
xlabel('\(x\)')
%ylabel('{\tt funappx\_g} error')
%gail.save_eps('TraubPaperOutput', 'funappx_g_errors');

h=semilogy( x, err, '-', x(ind), err(ind), '.');  hold off;
% semilogy( x, err, '-', 'color', MATLABBlue); 
% %axis tight;  
% hold on;
% semilogy( x(ind), err(ind), '.', 'color', MATLABOrange ); 
small = max(-20,log10(0.1*min(err)));
large = log10(10*max(err));
axis([a b 10^small 10^large])
xlabel('\(x\)')
%ylabel('Chebfun error')
legend(h,{'{\tt funappx\_g} error', 'Chebfun error'},'location', 'northwest','box','off')
set(gca,'ytick',10.^(5*ceil(small/5):5:5*floor(large/5)))
gail.save_eps('TraubPaperOutput', 'chebfun_errors');

% Example output:
%
% delta = 0.2; B = 1./(2*delta.^2); c = -0.2;
%   f3 = @(x) B*(4*delta.^2 + (x-c).^2 + (x-c-delta).*abs(x-c-delta) ...
%      - (x-c+delta).*abs(x-c+delta)).*(abs(x-c) <= 2*delta); a = - 1; b = 1; abstol = 1e-14;  
%   cf_chebfun(f3, a, b, abstol)
% fappx = 
%   griddedInterpolant with properties:
% 
%             GridVectors: {[1x45088972 double]}
%                  Values: [1x45088972 double]
%                  Method: 'linear'
%     ExtrapolationMethod: 'linear'
% fout = 
%            f: @(x)B*(4*delta.^2+(x-c).^2+(x-c-delta).*abs(x-c-delta)-(x-c+delta).*abs(x-c+delta)).*(abs(x-c)<=2*delta)
%            a: -1
%            b: 1
%       abstol: 1.000000000000000e-14
%          nlo: 10
%          nhi: 1000
%        ninit: 216
%         nmax: 100000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 20
%      npoints: 45088972
%       errest: 4.943971292550902e-15
% t1 =
%   49.791136166000001
% ---------------------
% fappx = 
%   griddedInterpolant with properties:
% 
%             GridVectors: {[1x30199129 double]}
%                  Values: [1x30199129 double]
%                  Method: 'linear'
%     ExtrapolationMethod: 'linear'
% fout = 
%            f: @(x)B*(4*delta.^2+(x-c).^2+(x-c-delta).*abs(x-c-delta)-(x-c+delta).*abs(x-c+delta)).*(abs(x-c)<=2*delta)
%            a: -1
%            b: 1
%       abstol: 1.000000000000000e-14
%          nlo: 20
%          nhi: 1000
%        ninit: 76
%         nmax: 100000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 21
%      npoints: 30199129
%       errest: 8.066471428701743e-15
% t2 =
%   25.345884486999999
% time_ratio1 =
%    1.964466309768675
% ---------------------
% c =
%    chebfun column (6 smooth pieces)
%        interval       length     endpoint values  
% [      -1,    -0.6]        1         0        0 
% [    -0.6,    -0.4]        3   8.7e-15      0.5 
% [    -0.4,-0.00078]        3       0.5      0.5 
% [-0.00078,  0.0047]       20       0.5     0.48 
% [  0.0047,     0.2]        3      0.48  3.1e-05 
% [     0.2,       1]        1         0        0 
% vertical scale =   1    Total length = 31
% t3 =
%    1.827106135000000
% time_ratio2 =
%   27.251364993090565
% --------------------
