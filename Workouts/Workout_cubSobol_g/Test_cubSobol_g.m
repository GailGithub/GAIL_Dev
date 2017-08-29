% TEST_CUBSOBOL_G This is the driver script to test cubSobol_g algorithm
%using seven integrands of dimensions up to 8
%clear all;close all;clc;
function [ut_abserr,ut_relerr,abstol,reltol] = Test_cubSobol_g
%[
dimsize = 3;
indexsize = 2;
abstol = 1e-3;
reltol = abstol;
format long
in_param.measure  = 'uniform';
disp('');
disp(horzcat('Dim  ', ' FcnIdx ',  '      Q    ','         f_true     ',...
    '          Err      ','      Sample Used    ', '         Stats  ')); 
disp(        '-----------------------------------------------------------------------------------------------------');
ut_abserr = nan(dimsize,indexsize);
ut_relerr = nan(dimsize,indexsize);
for dim=1:dimsize
  in_param.dim =dim;%the function dimension
  startingpoint = zeros(1,in_param.dim);%the lower limits of the integral
  endingpoint = ones(1,in_param.dim);%the upper limits of the integral
  hyperbox = [startingpoint;endingpoint];% the integration interval
  in_param.abstol = abstol;% the absolute tolerance
  in_param.reltol = reltol;% the relative tolerance
  in_param.alpha = 1e-2;% the uncertainty
  in_param.nSig = 1e4;% the sample size to estimate sigma
  in_param.n1 = 1e4;% the initial sample size to estimate Q
  %in_param.fudge =1.2;% standard deviation inflation factor
  in_param.timebudget = 300;% time budget
  in_param.nbudget = 1e10;% sample budget
  alpha = ones(1,in_param.dim);
  beta = 1./ (1:in_param.dim); 
  r=2; % three coefficients in genz_test_fun and genz_test_fun_true
  for index=1:indexsize % index refers to different integrands in genz_test_fun
    test_function = @(x)genz_test_fun(x,index,in_param.dim,alpha,beta,r);
    % the test function
    f_true = genz_test_fun_true (hyperbox,index,in_param.dim,alpha,beta,r);
    % true integral of the test function
    [Q,out_param]=cubSobol_g(test_function,hyperbox,in_param);
    % the results by using cubLattice_g
    abserr = abs(Q-f_true);% the absolute error
    relerr = abs((Q-f_true)/f_true);% the relative error
    numstr=horzcat(num2str(dim), '     ', num2str(index), '       ',...
        num2str(Q,'%10.5e'), '       ', num2str(f_true,'%10.5e'),...
        '       ', num2str(abserr,'%10.5e'),...
        '         ', num2str(out_param.n));
    % print the results
    if abserr > in_param.abstol && relerr > in_param.reltol,
    %if both absolute error and relative error does not meet tolerance
      disp([numstr,'            NoErrMet']);% mark it as "both err exceed"
    elseif abserr < in_param.abstol && relerr > in_param.reltol,
        % if only relative error does not meet the tolerance
        disp([numstr,'             AbsErrMet']);% mark it as "rel err exceed"
    elseif abserr > in_param.abstol && relerr < in_param.reltol,
        %if only the absolute error does not meet the tolerance
        disp([numstr,'            RelErrMet']);% mark it as "abs err exceed"
    else
        disp([numstr,'             BothErrMet']);% otherwise disp "OK"
    end
    ut_abserr(dim,index) = abserr;
    ut_relerr(dim,index) = relerr;
  end
end
end



%% The following output was obtained on 2014-10-13  by Lan Jiang
%  
% Dim   FcnIdx       Q            f_true             Error            Sample Used        status  
% ----------------------------------------------------------------------------------------------
% 1     1       8.41429e-01       8.41471e-01       4.16988e-05         975445             BothErrMet
% 1     2       7.85471e-01       7.85398e-01       7.25717e-05         1236277             BothErrMet
% 1     3       3.75266e-01       3.75000e-01       2.65681e-04         2129675             BothErrMet
% 1     4       7.46756e-01       7.46824e-01       6.76680e-05         1728700             BothErrMet
% 1     5       6.32161e-01       6.32121e-01       4.03169e-05         1462142             BothErrMet
% 1     6       1.71795e+00       1.71828e+00       3.29683e-04         2997628             BothErrMet
% 1     7       1.38083e+00       1.38039e+00       4.39163e-04         4339991             BothErrMet
% 2     1       4.96824e-01       4.96751e-01       7.24267e-05         3649273             BothErrMet
% 2     2       7.28163e-01       7.28296e-01       1.33307e-04         1179211             BothErrMet
% 2     3       1.01966e-01       1.01852e-01       1.14565e-04         707289     AbsErrMet
% 2     4       6.89122e-01       6.88992e-01       1.30139e-04         1604731             BothErrMet
% 2     5       4.97274e-01       4.97440e-01       1.65640e-04         1241193             BothErrMet
% Warning: At step 2, tried to evaluate at 27050414 samples, which is more than the remaining
% 16945370 samples. We will use all the sample left to estimate the mean. 
% > In meanMC_g>meanMC_g_err at 493
%   In meanMC_g>meanmctolfun at 243
%   In meanMC_g at 208
%   In cubMC_g at 211
%   In Test_cubMC_g at 30 
% 2     6       1.11489e+00       1.11469e+00       2.08576e-04         16965382             BothErrMet
% 2     7       1.80862e+00       1.80819e+00       4.32649e-04         10925304             BothErrMet
% 3     1       6.21012e-02       6.23593e-02       2.58075e-04         6206320     AbsErrMet
% 3     2       6.62509e-01       6.62570e-01       6.02562e-05         1193449             BothErrMet
% 3     3       2.17371e-02       2.17014e-02       3.57372e-05         160438     AbsErrMet
% 3     4       6.20919e-01       6.20903e-01       1.52024e-05         1596268             BothErrMet
% 3     5       3.82965e-01       3.83055e-01       8.95927e-05         1042167             BothErrMet
% Warning: At step 2, tried to evaluate at 27836133 samples, which is more than the remaining
% 18770218 samples. We will use all the sample left to estimate the mean. 
% > In meanMC_g>meanMC_g_err at 493
%   In meanMC_g>meanmctolfun at 243
%   In meanMC_g at 208
%   In cubMC_g at 211
%   In Test_cubMC_g at 30 
% 3     6       4.40912e-01       4.40984e-01       7.16979e-05         18790230             BothErrMet
% 3     7       2.16887e+00       2.16831e+00       5.59135e-04         28582950             BothErrMet
% 4     1       -3.52095e-01       -3.51764e-01       3.31374e-04         6517309             BothErrMet
% 4     2       5.88689e-01       5.88680e-01       9.02213e-06         1225452             BothErrMet
% 4     3       3.76752e-03       3.80556e-03       3.80333e-05         40451     AbsErrMet
% 4     4       5.43599e-01       5.43373e-01       2.25590e-04         1516358             BothErrMet
% 4     5       2.86912e-01       2.86844e-01       6.78729e-05         853151             BothErrMet
% 4     6       1.25321e-01       1.25251e-01       7.07149e-05         11047422             BothErrMet
% 4     7       2.16562e+00       2.16593e+00       3.11327e-04         93662781             BothErrMet
% 5     1       -6.49227e-01       -6.49331e-01       1.03821e-04         4553095             BothErrMet
% 5     2       5.13536e-01       5.13409e-01       1.27207e-04         1212257             BothErrMet
% 5     3       5.44026e-04       5.67130e-04       2.31041e-05         20012     AbsErrMet
% 5     4       4.64590e-01       4.64603e-01       1.29513e-05         1455926             BothErrMet
% 5     5       2.09824e-01       2.09952e-01       1.28909e-04         680797             BothErrMet
% 5     6       2.77841e-02       2.77308e-02       5.33210e-05         3955141     AbsErrMet
% Warning: At step 2, tried to evaluate at 363333179 samples, which is more than the remaining
% 156107201 samples. We will use all the sample left to estimate the mean. 
% > In meanMC_g>meanMC_g_err at 493
%   In meanMC_g>meanmctolfun at 243
%   In meanMC_g at 208
%   In cubMC_g at 211
%   In Test_cubMC_g at 30 
% 5     7       1.13488e+00       1.13532e+00       4.41697e-04         156127213             BothErrMet
% 6     1       -7.69317e-01       -7.69376e-01       5.96515e-05         2903377             BothErrMet
% 6     2       4.41269e-01       4.41474e-01       2.04762e-04         1134235             BothErrMet
% 6     3       8.33814e-05       7.34937e-05       9.88770e-06         20012     AbsErrMet
% 6     4       3.90190e-01       3.90227e-01       3.73760e-05         1342631             BothErrMet
% 6     5       1.51006e-01       1.50939e-01       6.71353e-05         508574             BothErrMet
% 6     6       5.35235e-03       5.02927e-03       3.23080e-04         998551     AbsErrMet
% Warning: At step 2, tried to evaluate at 369500174 samples, which is more than the remaining
% 138997141 samples. We will use all the sample left to estimate the mean. 
% > In meanMC_g>meanMC_g_err at 493
%   In meanMC_g>meanmctolfun at 243
%   In meanMC_g at 208
%   In cubMC_g at 211
%   In Test_cubMC_g at 30 
% 6     7       -2.32827e+00       -2.32730e+00       9.70403e-04         139017153             BothErrMet
% 7     1       -6.97541e-01       -6.97824e-01       2.82984e-04         4242250             BothErrMet
% 7     2       3.75388e-01       3.75484e-01       9.57243e-05         1067298             BothErrMet
% 7     3       8.01023e-06       8.42590e-06       4.15666e-07         20012     AbsErrMet
% 7     4       3.23245e-01       3.23235e-01       9.80560e-06         1195993             BothErrMet
% 7     5       1.06979e-01       1.06978e-01       1.26088e-06         359100             BothErrMet
% 7     6       8.80028e-04       7.72320e-04       1.07708e-04         434550     AbsErrMet
% 7     7       -1.10547e+01       -1.10568e+01       2.12771e-03         97127597     RelErrMet
% 8     1       -4.66861e-01       -4.67036e-01       1.74775e-04         7976409             BothErrMet
% 8     2       3.16567e-01       3.16602e-01       3.49663e-05         943146             BothErrMet
% 8     3       9.37483e-07       8.66209e-07       7.12746e-08         20012     AbsErrMet
% 8     4       2.64938e-01       2.64801e-01       1.37326e-04         1044264             BothErrMet
% 8     5       7.49463e-02       7.49531e-02       6.81415e-06         254752             BothErrMet
% 8     6       5.80123e-04       1.02833e-04       4.77290e-04         20012     AbsErrMet
% 8     7       -3.06063e+01       -3.06091e+01       2.80834e-03         41618674     RelErrMet