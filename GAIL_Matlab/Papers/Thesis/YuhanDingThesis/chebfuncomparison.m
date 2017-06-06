%% Chebfun stuff
clearvars, close all
 
% 
% xplot=-1:0.001:1;
% f0cheb=chebfun(@sin,[-1 1])
%  
% 
% f1=@(x) sin(x.^2); 
% f1cheb=chebfun(f1,[-1 1])
% figure, plot(xplot,[f1(xplot); f1cheb(xplot)],'linewidth',2)
%  
% 
% f2=@(x) sin(1000.*x.^2); 
% f2cheb=chebfun(f2,[-1 1])
% figure, plot(xplot,[f2(xplot); f2cheb(xplot)],'linewidth',2)
%  
% 
% f3=@(x) exp(-1e5*(x-0.1).^2)
% f3cheb=chebfun(f3,[-1 1])
% figure, plot(xplot,[f3(xplot); f3cheb(xplot)],'linewidth',2)
 
xplot=-1:0.00001:1;
f5=@(x) exp(1-0.095.^2./(x.*(0.19-x))).*(x>=0).*(x<=0.19);
tic;
f5cheb=chebfun(f5,[-1 1]);
toc;
tic;
[fappx,out_param]=funappxglobal_g(f5,-1,1,1e-5,10,100);
toc;
figure, plot(xplot,f5(xplot),'linewidth',2)
gail.save_eps('WorkoutFunappxOutput', 'chebcomp1');
figure, plot(xplot,f5cheb(xplot),'r+',xplot, ppval(fappx,xplot),'b-','linewidth',2)
legend('chebfun','funappxglobal\_g')
ylim([-0.1 1.1]);
gail.save_eps('WorkoutFunappxOutput', 'chebcomp2');
max(abs(f5cheb(xplot)-f5(xplot)))
max(abs(ppval(fappx,xplot)-f5(xplot)))