% TEST_CUBMC_G This is the driver script to test cubMC_g algorithm
%using seven integrands of dimensions up to 8
%clear all;close all;clc;
format long
in_param.measure  = 'uniform';
disp(horzcat('Dim  ', ' FcnIdx ',  ' Q            f_true         Error')); 
disp(        '------------------------------------------------------------');
for dim=1:4
  in_param.dim =dim;%the function dimension
  startingpoint = zeros(1,in_param.dim);
  endingpoint = ones(1,in_param.dim);
  hyperbox = [startingpoint;endingpoint];% the integration interval
  in_param.abstol = 1e-2;% the absolute tolerance
  in_param.alpha = 1e-2;% the uncertainty
  in_param.n_sigma = 1e4;% the sample size to estimate sigma
  in_param.fudge =1.1;% standard deviation inflation factor
  in_param.timebudget = 300;% time budget
  in_param.nbudget = 1e10;% sample budget
  alpha = ones(1,in_param.dim); 
  beta = 1./ (1:in_param.dim); 
  r=2;
  % three coefficients in genz_test_fun and genz_test_fun_true
  for index=[1:7]
    test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
    f_true = genz_test_fun_true (hyperbox,index,in_param.dim,alpha,beta,r);
    % true solution
    [Q,out_param]=cubMC_g(test_function,hyperbox,in_param);% the results using cubMC_g
    error = abs(Q-f_true);
    numstr=horzcat(num2str(dim), '     ', num2str(index), '       ',...
        num2str(Q), '       ', num2str(f_true),'       ', num2str(error));
    if error > in_param.abstol,% if error does not meet tolerance, mark it
      disp([numstr,'     ****']);
    else
      disp(numstr);
    end
  end
end

%% The following output was obtained on 2014-02-28  by Lan Jiang
%  
% Dim FcnIdx  Q            f_true         Error
% ---------------------------------------------------
% 1   1       0.84115       0.84147       0.00031717
% 1   2       0.78502       0.7854       0.00037834
% 1   3       0.37513       0.375       0.00012598
% 1   4       0.74974       0.74682       0.002917
% 1   5       0.63327       0.63212       0.0011525
% 1   6       1.7159       1.7183       0.0023409
% 1   7       1.3782       1.3804       0.0021664
% 2   1       0.49586       0.49675       0.00089481
% 2   2       0.7273       0.7283       0.00099608
% 2   3       0.10216       0.10185       0.00030964
% 2   4       0.68868       0.68899       0.00030784
% 2   5       0.49672       0.49744       0.00072303
% 2   6       1.114       1.1147       0.00070263
% 2   7       1.8063       1.8082       0.0018937
% 3   1       0.062137       0.062359       0.00022213
% 3   2       0.66103       0.66257       0.001542
% 3   3       0.021925       0.021701       0.00022398
% 3   4       0.62161       0.6209       0.00070482
% 3   5       0.38442       0.38305       0.0013682
% 3   6       0.43699       0.44098       0.0039984
% 3   7       2.1689       2.1683       0.00055366
% 4   1       -0.35004       -0.35176       0.0017239
% 4   2       0.5874       0.58868       0.0012771
% 4   3       0.0036658       0.0038056       0.00013974
% 4   4       0.54645       0.54337       0.0030758
% 4   5       0.2864       0.28684       0.00044835
% 4   6       0.1271       0.12525       0.0018505
% 4   7       2.1698       2.1659       0.0039153
% 5   1       -0.64743       -0.64933       0.0018976
% 5   2       0.51271       0.51341       0.00069956
% 5   3       0.00057889       0.00056713       1.1762e-05
% 5   4       0.46137       0.4646       0.0032351
% 5   5       0.20932       0.20995       0.0006327
% 5   6       0.026101       0.027731       0.0016302
% 5   7       1.1374       1.1353       0.0020326
% 6   1       -0.7684       -0.76938       0.00097805
% 6   2       0.44013       0.44147       0.0013459
% 6   3       7.3767e-05       7.3494e-05       2.7338e-07
% 6   4       0.39163       0.39023       0.0014033
% 6   5       0.15102       0.15094       8.411e-05
% 6   6       0.0050143       0.0050293       1.4935e-05
% 6   7       -2.3197       -2.3273       0.007637
% 7   1       -0.69845       -0.69782       0.0006228
% 7   2       0.37672       0.37548       0.0012361
% 7   3       8.1634e-06       8.4259e-06       2.6245e-07
% 7   4       0.324       0.32324       0.00076769
% 7   5       0.10689       0.10698       9.0749e-05
% 7   6       0       0.00077232       0.00077232
% 7   7       -11.0614       -11.0568       0.0045222
% 8   1       -0.46467       -0.46704       0.0023611
% 8   2       0.31432       0.3166       0.0022834
% 8   3       1.0008e-06       8.6621e-07       1.3457e-07
% 8   4       0.26724       0.2648       0.0024348
% 8   5       0.075203       0.074953       0.00024963
% 8   6       0       0.00010283       0.00010283
% Warning: tried to evaluate at 147596348
% samples, which is more than the allowed
% maximum of 111232448 samples. Just use
% the maximum sample budget. 
% > In meanMC_g>meanMC_g_err at 433
%   In meanMC_g at 246
%   In cubMC_g at 186
%   In test_cubMC_g at 26 
% 8   7       -30.612       -30.6091       0.0029272