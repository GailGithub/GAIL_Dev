% TEST_CUBSobol This is the driver script to test cubSobol algorithm
%using integrands of various dimensions
clear all;close all;clc;
in_param.measure  = 'uniform';
disp(horzcat('Dim FcnIdx     Q             Error        time        n')); 
disp(        '--------------------------------------------------------------------');
for dim=1:8
  abstol = 1e-4;% the absolute tolerance
  %alpha = ones(1,dim); 
  alpha = 2*rand(1,dim); 
  beta = 1./ (1:dim).*(2*rand(1,dim)); 
  r=2;
  % three coefficients in genz_test_fun and genz_test_fun_true
  for index=1:7
    true_integ = genz_test_fun_true ([zeros(1,dim);ones(1,dim)],index,dim,alpha,beta,r);
    % true solution
    test_function = @(x)genz_test_fun(x,index,dim,alpha,beta,r)./true_integ;
    [Q,err,time,n,overbudget]=cubSobol(test_function,dim,abstol);% the results using cubMC_g
    error = abs(Q-1);
    numstr=horzcat(int2str(dim), '   ', int2str(index), '    ', ...
       num2str(Q,'%+13.6e'), '    ', num2str(error,'%13.4e'),...
       '    ', num2str(time,'%13.6f'),'    ', int2str(n));
    if overbudget,% overbudget
      numstr=horzcat(numstr,'XXX');
    end
    if error > abstol,% if error does not meet tolerance, mark it
      numstr=horzcat(numstr,'   failed');
    end
    disp(numstr)
  end
end