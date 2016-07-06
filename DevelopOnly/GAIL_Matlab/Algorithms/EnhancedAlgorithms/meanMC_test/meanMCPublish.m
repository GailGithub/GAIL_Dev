%% Documentation for meanMC
% meanMC is a class that uses the IID Monte Carlo method to estimate the
% mean of a random variable.
%
% meanMC offers the option for users to use the following variance reduction methods:
% 1. control variate with coefficients estimated by linear regression
% 2. control variate with coefficients estimated by ridge regression
% and 3. antithetic variate
%
% meanMC is a vectorized class, which means you can create an array of
% meanMC 1ects and estimate their means at one go.

%% No variance reduction
% In this example, we are interested in estimating E[sin(U)], where U is
% uniformly distributed over [0, 1]. No variance reduction is applied.
obj1.in_param.abstol = 2e-2; % set absolute error tolerance
obj1.in_param.reltol = 2e-1; % set relative error tolerance
obj1.method = {'plain'}; % no variance reduction 
obj1.Yrand = @(n) sin(rand(n, 1)); 
test1 = meanMC(obj1);
[mu param] = genMu(test1)
%% Control variate
% In this example, we are interested in estimating E[sin(U)], where U is
% uniformly distributed over [0, 1]. We use U as the control variate.
obj2.method = {'cv'};
obj2.cv_param.YXrand = @sinr_cv;
obj2.cv_param.muX = 0.5;
obj2.cv_param.ncv = 1.1e3; 
obj2.cv_param.ridge = 1e-3;
test2 = meanMC(obj2);
[mu param] = genMu(test2)
%% Antithetic variate
% In this example, we are interested in estimating E[sin(U)], where U is
% uniformly distributed over [0, 1]. We use sin(1-U) as the antithetic variate.
obj3.method = {'av'};
obj3.av_param.YYrand = @sinr_av;
test3 = meanMC(obj3);
[mu param] = genMu(test3)
%% Combined methods
% In this example, we are interested in estimating E[sin(U)], where U is
% uniformly distributed over [0, 1].
%
% This method will automatically choose the optimal method from the list of
% methods the user suggests. 
obj4.method = {'cv','av','plain'};
obj4.nc = 1e3; % number of samples (per method) used to compare different methods
obj4.Yrand = @(n) sin(rand(n, 1));
obj4.cv_param.YXrand = @sinr_cv;
obj4.cv_param.muX = 0.5;
obj4.cv_param.ridge = [1e-3 2e-3]; 
% use ridge regression to estimate control variate coefficients
% if multiple ridge coefficients are suggested, this will be interpreted as
% two different methods (since control variate coefficients are estimated
% in 2 different ways)
obj4.av_param.YYrand = @sinr_av;
test4 = meanMC(obj4);
[mu param] = genMu(test4)
%% Multiple problems at one go
% In this example, we are interested in estimating E[sin(U)] and E[exp(U)] at the same time, 
% where U is uniformly distributed over [0, 1]
% you can specify different values for all numerical parameters in two
% ways: give an array or give a cell array.
%
% Exception: You can only use cell arrays to indicate different values 
% for cv_param.muX and cv_param.ridge to avoid possible confusion.
%
% All other paramters need to be given in the form of a cell array when
% you want to specify them seperately for each problem.
obj5.in_param.abstol = [1e-2 2e-2]; % different abstol for the two problems
% alternatively obj5.in_param.abstol = {1e-2 2e-2};
obj5.in_param.reltol = 0;
obj5.in_param.alpha = [1e-2 2e-2]; % different alpha for the two problems
obj5.method = {'cv','av','plain'};
obj5.nc = 1e3;
obj5.Yrand = {@(n) sin(rand(n, 1)), @(n) exp(rand(n, 1))};
obj5.cv_param.YXrand = {@sinr_cv, @expr_cv};
obj5.cv_param.muX = 0.5;
obj5.cv_param.ncv = 1e3;
obj5.cv_param.ridge = {[0 1e-3] 0}; % different ridge coefficients
obj5.av_param.YYrand = {@sinr_av, @expr_av};
test = meanMC(obj5);
test5 = meanMC(obj5);
[mu param] = genMu(test5);
mu
param{1}
param{2}