% TEST_CUBMC_G This is the driver script to test cubMC_g algorithm
%using seven integrands of dimensions up to 8
%clear all;close all;clc;
format long
in_param.measure  = 'uniform';
disp(horzcat('Dim  ', ' FcnIdx ',  '      Q    ','        f_true     ',...
    '        Error      ','      Sample Used    ', '    status  ')); 
disp(        '----------------------------------------------------------------------------------------------');
for dim=1:8
  in_param.dim =dim;%the function dimension
  startingpoint = zeros(1,in_param.dim);%the lower limits of the integral
  endingpoint = ones(1,in_param.dim);%the upper limits of the integral
  hyperbox = [startingpoint;endingpoint];% the integration interval
  in_param.abstol = 1e-3;% the absolute tolerance
  in_param.reltol = 1e-3;% the relative tolerance
  in_param.alpha = 1e-2;% the uncertainty
  in_param.nSig = 1e4;% the sample size to estimate sigma
  in_param.n1 = 1e4;% the initial sample size to estimate Q
  in_param.fudge =1.2;% standard deviation inflation factor
  in_param.timebudget = 300;% time budget
  in_param.nbudget = 1e10;% sample budget
  alpha = ones(1,in_param.dim);
  beta = 1./ (1:in_param.dim); 
  r=2; % three coefficients in genz_test_fun and genz_test_fun_true
  for index=1:7 % index refers to different integrands in genz_test_fun
    test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
    % the test function
    f_true = genz_test_fun_true (hyperbox,index,in_param.dim,alpha,beta,r);
    % true integral of the test function
    [Q,out_param]=cubMC_g(test_function,hyperbox,in_param);
    % the results by using cubMC_g
    abserr = abs(Q-f_true);% the absolute error
    relerr = abs((Q-f_true)/f_true);% the relative error
    numstr=horzcat(num2str(dim), '     ', num2str(index), '       ',...
        num2str(Q,'%10.5e'), '       ', num2str(f_true,'%10.5e'),...
        '       ', num2str(abserr,'%10.5e'),...
        '         ', num2str(out_param.ntot));
    % print the results
    if abserr > in_param.abstol && relerr > in_param.reltol,
    %if both absolute error and relative error does not meet tolerance
      disp([numstr,'     NoErrMet']);% mark it as "both err exceed"
    elseif abserr < in_param.abstol && relerr > in_param.reltol,
        % if only relative error does not meet the tolerance
        disp([numstr,'     AbsErrMet']);% mark it as "rel err exceed"
    elseif abserr > in_param.abstol && relerr < in_param.reltol,
        %if only the absolute error does not meet the tolerance
        disp([numstr,'     RelErrMet']);% mark it as "abs err exceed"
    else
        disp([numstr,'             BothErrMet']);% otherwise disp "OK"
    end
  end
end

%% The following output was obtained on 2014-10-04  by Lan Jiang
%  
% Dim   FcnIdx       Q            f_true             Error            Sample Used        status  
% ----------------------------------------------------------------------------------------------
% 1     1       8.41343e-01       8.41471e-01       1.28378e-04         987270             OK
% 1     2       7.85419e-01       7.85398e-01       2.10586e-05         1240479             OK
% 1     3       3.75075e-01       3.75000e-01       7.54451e-05         2123305             OK
% 1     4       7.46753e-01       7.46824e-01       7.12502e-05         1717254             OK
% 1     5       6.32131e-01       6.32121e-01       1.05773e-05         1470892             OK
% 1     6       1.71827e+00       1.71828e+00       9.18923e-06         2926898             OK
% 1     7       1.38024e+00       1.38039e+00       1.44901e-04         4288514             OK
% 2     1       4.96534e-01       4.96751e-01       2.17433e-04         3674991             OK
% 2     2       7.28306e-01       7.28296e-01       9.94141e-06         1198115             OK
% 2     3       1.02018e-01       1.01852e-01       1.65905e-04         735294     RelErrExceed
% 2     4       6.88943e-01       6.88992e-01       4.82295e-05         1619323             OK
% 2     5       4.97712e-01       4.97440e-01       2.71680e-04         1224954             OK
% Warning: tried to evaluate at 26284672 samples, which is more than the allowed maximum of 19195890
% samples. Just use the maximum sample budget. 
% > In meanMC_g>meanMC_g_err at 476
%   In meanMC_g>meanmctolfun at 241
%   In meanMC_g at 206
%   In cubMC_g at 207
%   In Test_cubMC_g at 30 
% 2     6       1.11448e+00       1.11469e+00       2.07804e-04         19215912             OK
% 2     7       1.80842e+00       1.80819e+00       2.33968e-04         11062714             OK
% 3     1       6.25176e-02       6.23593e-02       1.58278e-04         6125806     RelErrExceed
% 3     2       6.62386e-01       6.62570e-01       1.83902e-04         1206712             OK
% 3     3       2.18969e-02       2.17014e-02       1.95497e-04         156042     RelErrExceed
% 3     4       6.20775e-01       6.20903e-01       1.28494e-04         1567173             OK
% 3     5       3.82914e-01       3.83055e-01       1.40058e-04         1042065             OK
% Warning: tried to evaluate at 29202842 samples, which is more than the allowed maximum of 23732368
% samples. Just use the maximum sample budget. 
% > In meanMC_g>meanMC_g_err at 476
%   In meanMC_g>meanmctolfun at 241
%   In meanMC_g at 206
%   In cubMC_g at 207
%   In Test_cubMC_g at 30 
% 3     6       4.41391e-01       4.40984e-01       4.07293e-04         23752390             OK
% 3     7       2.16808e+00       2.16831e+00       2.26099e-04         27881870             OK
% 4     1       -3.51893e-01       -3.51764e-01       1.28932e-04         6565204             OK
% 4     2       5.88906e-01       5.88680e-01       2.26208e-04         1201652             OK
% 4     3       3.81270e-03       3.80556e-03       7.14627e-06         46144     RelErrExceed
% 4     4       5.43207e-01       5.43373e-01       1.66097e-04         1516906             OK
% 4     5       2.86845e-01       2.86844e-01       1.55998e-06         838936             OK
% 4     6       1.25211e-01       1.25251e-01       3.97478e-05         10934471             OK
% 4     7       2.16650e+00       2.16593e+00       5.70649e-04         87430732             OK
% 5     1       -6.49099e-01       -6.49331e-01       2.31641e-04         4551358             OK
% 5     2       5.13471e-01       5.13409e-01       6.19648e-05         1212369             OK
% 5     3       5.54533e-04       5.67130e-04       1.25967e-05         20022     RelErrExceed
% 5     4       4.64712e-01       4.64603e-01       1.09414e-04         1445830             OK
% 5     5       2.10056e-01       2.09952e-01       1.03099e-04         651357             OK
% 5     6       2.77204e-02       2.77308e-02       1.04365e-05         3558925             OK
% Warning: tried to evaluate at 351256589 samples, which is more than the allowed maximum of
% 207329933 samples. Just use the maximum sample budget. 
% > In meanMC_g>meanMC_g_err at 476
%   In meanMC_g>meanmctolfun at 241
%   In meanMC_g at 206
%   In cubMC_g at 207
%   In Test_cubMC_g at 30 
% 5     7       1.13577e+00       1.13532e+00       4.49289e-04         207349955             OK
% 6     1       -7.69107e-01       -7.69376e-01       2.69869e-04         2809603             OK
% 6     2       4.41452e-01       4.41474e-01       2.14539e-05         1154627             OK
% 6     3       6.71034e-05       7.34937e-05       6.39027e-06         20022     RelErrExceed
% 6     4       3.90184e-01       3.90227e-01       4.34789e-05         1356084             OK
% 6     5       1.50864e-01       1.50939e-01       7.51571e-05         507202             OK
% 6     6       4.98953e-03       5.02927e-03       3.97454e-05         1292676     RelErrExceed
% Warning: tried to evaluate at 412039854 samples, which is more than the allowed maximum of
% 169339672 samples. Just use the maximum sample budget. 
% > In meanMC_g>meanMC_g_err at 476
%   In meanMC_g>meanmctolfun at 241
%   In meanMC_g at 206
%   In cubMC_g at 207
%   In Test_cubMC_g at 30 
% 6     7       -2.32639e+00       -2.32730e+00       9.18060e-04         169359694             OK
% 7     1       -6.97792e-01       -6.97824e-01       3.16745e-05         4254664             OK
% 7     2       3.75514e-01       3.75484e-01       3.04815e-05         1097426             OK
% 7     3       7.87780e-06       8.42590e-06       5.48093e-07         20022     RelErrExceed
% 7     4       3.23409e-01       3.23235e-01       1.73487e-04         1211068             OK
% 7     5       1.07006e-01       1.06978e-01       2.81095e-05         360319             OK
% 7     6       9.11923e-04       7.72320e-04       1.39602e-04         480750     RelErrExceed
% 7     7       -1.10596e+01       -1.10568e+01       2.71187e-03         109936324     AbsErrExceed
% 8     1       -4.66671e-01       -4.67036e-01       3.65306e-04         7898153             OK
% 8     2       3.16850e-01       3.16602e-01       2.48055e-04         964833             OK
% 8     3       7.98513e-07       8.66209e-07       6.76959e-08         20022     RelErrExceed
% 8     4       2.64834e-01       2.64801e-01       3.27813e-05         1030309             OK
% 8     5       7.48311e-02       7.49531e-02       1.22054e-04         251055     RelErrExceed
% 8     6       0.00000e+00       1.02833e-04       1.02833e-04         20022     RelErrExceed
% 8     7       -3.06156e+01       -3.06091e+01       6.56798e-03         41933007     AbsErrExceed
