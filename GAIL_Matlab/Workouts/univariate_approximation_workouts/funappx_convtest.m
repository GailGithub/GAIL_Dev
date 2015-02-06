function [npoints,errest,t,npointsglobal,errestglobal,tglobal]=funappx_convtest()
format compact
tol = zeros(1,15);
errest = tol;
npoints = tol;
t = tol;
errestglobal = tol;
npointsglobal = tol;
tglobal = tol;
j=1;
%f = @(x) x.^2; 
f = @(x) exp(-100*(x-sqrt(2)/2).^2);
warning('off','MATLAB:funappx_g:peaky')
warning('off','MATLAB:funappx_g:exceedbudget')
warning('off','MATLAB:funappxglobal_g:peaky')
warning('off','MATLAB:funappxglobal_g:exceedbudget')

for i=-15:-1,
  tol(j) = 10^(i);
  in_param.abstol = 10^(i);
  tic,
  [~, out_param] = funappx_g(f, in_param);
  t(j) = toc;
  errest(j) = out_param.errest;
  npoints(j) = out_param.npoints;
  
  tic,
  [~, out_param] = funappxglobal_g(f, in_param);
  tglobal(j) = toc;
  errestglobal(j) = out_param.errorbound;
  npointsglobal(j) = out_param.npoints;
  
  j=j+1;
end

warning('on','MATLAB:funappxglobal_g:exceedbudget')
warning('on','MATLAB:funappx_g:peaky')
warning('on','MATLAB:funappx_g:exceedbudget')
warning('on','MATLAB:funappxglobal_g:peaky')

figure(1)
subplot(2,1,1)
loglog(npoints, errest)
title('Time and Computational Cost of funappx\_g VS error tolerance')
ylabel('error estimation')
xlabel('# of points')

subplot(2,1,2)
loglog(t, errest)
ylabel('error estimation')
xlabel('time cost')

gail.save_eps('WorkoutFunappxOutput', 'WorkoutFunAppxConvTest1');

figure(2)
subplot(2,1,1)
loglog(npointsglobal, errestglobal)
title('Time and Computational Cost of funappxglobal\_g VS error tolerance')
ylabel('error estimation')
xlabel('# of points')

subplot(2,1,2)
loglog(tglobal, errestglobal)
gail.save_eps('WorkoutFunappxOutput', 'WorkoutFunAppxConvTest2');
ylabel('error estimation')
xlabel('time cost')
end
