%TEST_CUBMC_G   This is the driver script to test cubMC_g algorithm
%using integrands of dimensions up to 3 
%clear all;close all;clc;
in_param.measure  = 'uniform';
disp(horzcat('Dim ', 'FcnIdx ',   'Error')); 
disp(        '---------------------------');
%f=@(x) exp(-x(1).^2-x(2).^2);% the test function
%f=@(x) x(:,1)+x(:,2);
for dim=1:3
  in_param.dim =dim;%the function dimension
  startingpoint = zero(1,in_param.dim);
  endingpoint = ones(1,in_param.dim);
  hyperbox = [startingpoint;endingpoint];% the integration interval
  in_param.abstol = 1e-3;% the absolute tolerance
  in_param.alpha = 0.01;% the uncertainty
  in_param.n_sigma = 1e4;% the sample size to estimate sigma
  in_param.fudge =1.1;% standard deviation inflation factor
  in_param.timebudget = 100;% time budget
  in_param.nbudget = 1e9;% sample budget
  alpha = ones(1:in_param.dim); % one coefficient
  beta = [1/3 1/4 2];
  %beta = ones(1:in_param.dim);% the other coefficent
  r=1;
  for index= 7;%[2,4:7]
    test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
    f_true = genz_test_fun_true (hyperbox,index,in_param.dim,alpha,beta);
    % true solution
    [Q,out_param]=cubMC_g(test_function,hyperbox,in_param);% the results using cubMC_g
    error = abs(Q-f_true);
    numstr=horzcat(num2str(dim), '   ', num2str(index), '       ', num2str(error));
    if error > in_param.abstol,
      disp([numstr,'   ****']);
    else
      disp(numstr);
    end
  end
end
