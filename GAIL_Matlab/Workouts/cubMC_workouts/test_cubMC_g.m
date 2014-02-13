% TEST_CUBMC_G This is the driver script to test cubMC_g algorithm
%using integrands of dimensions up to 3 
%clear all;close all;clc;
in_param.measure  = 'uniform';
disp(horzcat('Dim ', 'FcnIdx ',   ' Error')); 
disp(        '---------------------------');
for dim=1:7
  in_param.dim =dim;%the function dimension
  startingpoint = zeros(1,in_param.dim);
  endingpoint = ones(1,in_param.dim);
  hyperbox = [startingpoint;endingpoint];% the integration interval
  in_param.abstol = 1e-2;% the absolute tolerance
  in_param.alpha = 0.01;% the uncertainty
  in_param.n_sigma = 1e4;% the sample size to estimate sigma
  in_param.fudge =1.1;% standard deviation inflation factor
  in_param.timebudget = 100;% time budget
  in_param.nbudget = 1e9;% sample budget
  alpha = ones(1,in_param.dim); 
  beta = 1./ (1:in_param.dim); 
  r=2;
  % three coefficients in genz_test_fun and genz_test_fun_true
  for index=[1:6]
    test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
    f_true = genz_test_fun_true (hyperbox,index,in_param.dim,alpha,beta,r);
    % true solution
    [Q,out_param]=cubMC_g(test_function,hyperbox,in_param);% the results using cubMC_g
    error = abs(Q-f_true);
    numstr=horzcat(num2str(dim), '   ', num2str(index), '       ', num2str(error));
    if error > in_param.abstol,% if error does not meet tolerance, mark it
      disp([numstr,'   ****']);
    else
      disp(numstr);
    end
  end
end

%% The following output was obtained on 2014-02-11  by Lan Jiang
%  
% Test_cubMC_g
% Dim FcnIdx  Error
% ---------------------------
% 1   1       0.00066017
% 1   2       0.0029365
% 1   3       0.0013206
% 1   4       0.00026806
% 1   5       0.0027263
% 1   6       0.00039949
% 1   7       0.00056086
% 2   1       0.00093275
% 2   2       0.00036589
% 2   3       0.00072042
% 2   4       0.0013418
% 2   5       0.00067463
% 2   6       0.00062401
% 2   7       0.001117
% 3   1       0.0019781
% 3   2       0.00050505
% 3   3       0.00037568
% 3   4       0.0012184
% 3   5       4.7978e-05
% 3   6       0.00010979
% 3   7       0.0016297
