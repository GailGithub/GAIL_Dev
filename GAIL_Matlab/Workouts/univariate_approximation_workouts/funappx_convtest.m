clear all; clf; format compact
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
for i=-15:-1,
  tol(j) = 10^(i)
  in_param.abstol = 10^(i);
  tic,
  [~, out_param] = funappx_g(f, in_param);
  t(j) = toc
  errest(j) = out_param.errest   
  npoints(j) = out_param.npoints
  
  tic,
  [~, out_param] = funappxglobal_g(f, in_param);
  tglobal(j) = toc
  errestglobal(j) = out_param.errorbound   
  npointsglobal(j) = out_param.npoints
  
  j=j+1;
end

figure(1)
subplot(2,1,1)
loglog(npoints, errest)

subplot(2,1,2)
loglog(t, errest)

figure(2)
subplot(2,1,1)
loglog(npointsglobal, errestglobal)

subplot(2,1,2)
loglog(tglobal, errestglobal)