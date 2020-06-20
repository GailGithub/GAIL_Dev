%DT_CUBMC_G  lengthy doctests for cubMC_g
%
% Example 1: 
% Estimate the integral with integrand Genz "Oscillatory" in
% genz_test_fun, index 1.
% >>  in_param.dim = 3;alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 1;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
% f_true = 0.06***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-2)+0.001
% Q = 0.06***
% 
%
% Example 2: 
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 2.
% >> index = 2;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
% f_true = 0.66***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-3)+0.001
% Q = 0.66***
% 
%
% Example 3: 
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 3.
% >> index = 3;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
% f_true = 0.04***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-2)+0.001
% Q = 0.04***
% 
%
% Example 4: 
% Estimate the integral with integrand Genz "Gaussian" in
% genz_test_fun, index 4.
% >> index = 4;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.001
% f_true = 0.62***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-3)+0.001
% Q = 0.62***
%
%
% Example 5: 
% Estimate the integral with integrand Genz "Continuous" in
% genz_test_fun, index 5.
% >> index = 5;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.38***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-3)
% Q = 0.38***
%
%
% Example 6: 
% Estimate the integral with integrand Genz "Discontinuous" in
% genz_test_fun, index 6.
% >> index = 6;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.44***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-3)
% Q = 0.44***
%
%
% Example 7: 
% Estimate the integral with integrand Keister test function in
% genz_test_fun, index 7.
% >> index = 7; in_param.dim = 3;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-4)
% Q = 2.16***
%
%
% Example 8: 
% Estimate the integral with integrand Genz "Oscillatory" in
% genz_test_fun, index 1.
% >> in_param.dim = 3; alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
% >> hyperbox = [-ones(1,in_param.dim);ones(1,in_param.dim)];index = 1;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.59***
% >> abstol = 1e-2; reltol = 1e-13;
% >> Q = cubMC_g(f,hyperbox,'uniform',abstol,reltol);
% >> abs(f_true-Q) <= max(abstol, reltol * abs(f_true))
%    1
%
%
% Example 9: 
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 2.
% >> index = 2;hyperbox = [-ones(1,in_param.dim);ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.30***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-2,1e-13)
% Q = 0.30***
% 
%
% Example 10: 
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 3.
% >> index = 3;hyperbox = [ones(1,in_param.dim)/2;2*ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.00***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-2,1e-3)
% Q = 0.00***
% 
%
% Example 11: 
% Estimate the integral with integrand Genz "Gaussian" in
% genz_test_fun, index 4.
% >> index = 4;hyperbox = [-ones(1,in_param.dim)/2;ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.005
% f_true = 0.38***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-2,1e-3)+0.005;
% >> check = abs(Q - f_true)< max(1e-2, 1e-3*f_true)
%   1
%
%
% Example 12: 
% Estimate the integral with integrand Genz "Continuous" in
% genz_test_fun, index 5.
% >> index = 5;hyperbox = [-ones(1,in_param.dim)/2;ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)+0.005
% f_true = 0.25***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-2,1e-2)+0.005
% Q = 0.25***
%
%
% Example 13: 
% Estimate the integral with integrand Genz "Discontinuous" in
% genz_test_fun, index 6.
% >> index = 6;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.51***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,1e-3)
% Q = 0.51***
%
%