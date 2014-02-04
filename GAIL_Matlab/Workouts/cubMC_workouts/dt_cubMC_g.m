%DT_CUBMC_G  lengthy doctests for cubMC_g
%
% Example 1: 
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 2.
% >>  in_param.dim = 3;alpha = ones(1:in_param.dim); beta = [1/3 1/4 2];r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 2;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta)+0.001
% f_true = 0.26***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3)+0.001
% Q = 0.26***
% 
%
% Example 2: 
% Estimate the integral with integrand Genz "Gaussian" in
% genz_test_fun, index 4.
% >>  in_param.dim = 3;alpha = ones(1:in_param.dim); beta = [1/3 1/4 2];r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 4;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta)+0.001
% f_true = 0.10***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3)+0.001
% Q = 0.10***
%
%
% Example 3: 
% Estimate the integral with integrand Genz "Continuous" in
% genz_test_fun, index 5.
% >>  in_param.dim = 3;alpha = ones(1:in_param.dim); beta = [1/3 1/4 2];r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 5;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta)
% f_true = 0.13***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3)
% Q = 0.13***
%
%
% Example 4: 
% Estimate the integral with integrand Genz "Discontinuous" in
% genz_test_fun, index 6.
% >>  in_param.dim = 3;alpha = ones(1:in_param.dim); beta = [1/3 1/4 2];r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 6;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta)
% f_true = 0.19***
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3)
% Q = 0.19***
%
%
% Example 5: 
% Estimate the integral with integrand Kesiter test function in
% genz_test_fun, index 7.
% >>  in_param.dim = 3;alpha = ones(1:in_param.dim); beta = [1/3 1/4 2];r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 7;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3)
% Q = 0.61***
%
%