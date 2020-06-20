%DT_CUBLATTICE_G  lengthy doctests for cubLattice_g
%
% Example 1:
% Estimate the integral with integrand Genz "Oscillatory" in
% genz_test_fun, index 1.
% >> in_param.dim = 3;alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
% >> hyperbox = [zeros(1,in_param.dim);ones(1,in_param.dim)];index = 1;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.06***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-2,1,f_true,'max'))
% check = 1
%
%
% Example 2:
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 2.
% >> index = 2;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.66***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-3);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
% check = 1
%
%
% Example 3:
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 3.
% >> index = 3;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.04***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-2,1,f_true,'max'))
% check = 1
%
%
% Example 4:
% Estimate the integral with integrand Genz "Gaussian" in
% genz_test_fun, index 4.
% >> index = 4;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.62***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-3);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
% check = 1
%
%
% Example 5:
% Estimate the integral with integrand Genz "Continuous" in
% genz_test_fun, index 5.
% >> index = 5;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.38***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-3);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-3,1e-3,1,f_true,'max'))
% check = 1
%
%
% Example 6:
% Estimate the integral with integrand Genz "Oscillatory" in
% genz_test_fun, index 1.
% >> in_param.dim = 3;alpha = ones(1,in_param.dim); beta = 1./ (1:in_param.dim);r=1;
% >> hyperbox = [-ones(1,in_param.dim);ones(1,in_param.dim)];index = 1;
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.59***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-2,1e-13);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-2,1e-13,1,f_true,'max'))
% check = 1
%
%
% Example 7:
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 2.
% >> index = 2;hyperbox = [-ones(1,in_param.dim);ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.30***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-2,1e-13);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-2,1e-13,1,f_true,'max'))
% check = 1
%
%
% Example 8:
% Estimate the integral with integrand Genz "Product Peak" in
% genz_test_fun, index 3.
% >> index = 3;hyperbox = [ones(1,in_param.dim)/2;2*ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.00***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-2,1e-3);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-2,1e-3,1,f_true,'max'))
% check = 1
%
%
% Example 9:
% Estimate the integral with integrand Genz "Gaussian" in
% genz_test_fun, index 4.
% >> index = 4;hyperbox = [-ones(1,in_param.dim)/2;ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.38***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-2,1e-3);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-2,1e-3,1,f_true,'max'))
% check = 1
%
%
% Example 10:
% Estimate the integral with integrand Genz "Continuous" in
% genz_test_fun, index 5.
% >> index = 5;hyperbox = [-ones(1,in_param.dim)/2;ones(1,in_param.dim)];
% >> f = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
% >> f_true = genz_test_fun_true(hyperbox,index,in_param.dim,alpha,beta,r)
% f_true = 0.24***
% >> Q = cubLattice_g(f,hyperbox,'uniform',1e-2,1e-2);
% >> check = double(abs(f_true-Q) < gail.tolfun(1e-2,1e-2,1,f_true,'max'))
% check = 1
%
%
