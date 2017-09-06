% TEST_CUBSobol This is the driver script to test cubSobol algorithm
%using integrands of dimensions up to 3 
clear all;close all;clc;
in_param.measure  = 'uniform';
disp(horzcat('Dim ', 'FcnIdx ',  ' Q          f_true       Error    time    n')); 
disp(        '--------------------------------------------------------------------');
for dim=1:8
  abstol = 1e-4;% the absolute tolerance
  alpha = ones(1,dim); 
  beta = 1./ (1:dim); 
  r=2;
  % three coefficients in genz_test_fun and genz_test_fun_true
  for index=[1:7]
    test_function = @(x)genz_test_fun(x,index,dim,alpha,beta,r);
    f_true = genz_test_fun_true ([zeros(1,dim);ones(1,dim)],index,dim,alpha,beta,r);
    % true solution
    [Q,err,time,n]=cubSobol(test_function,dim,abstol);% the results using cubMC_g
    error = abs(Q-f_true);
    numstr=horzcat(int2str(dim), '   ', int2str(index), '    ', ...
       num2str(Q,'%+13.6e'), '    ', num2str(f_true,'%+13.6e'),'    ', ...
       num2str(error,'%13.4e'),...
       '    ', num2str(time,'%13.6f'),'    ', int2str(n));
    if error > abstol,% if error does not meet tolerance, mark it
      disp([numstr,'   ****']);
    else
      disp(numstr);
    end
  end
end