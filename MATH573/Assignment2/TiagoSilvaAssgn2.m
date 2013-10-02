% Name: Tiago Silva
% Email: tsilva@hawk.iit.edu
%
% PC OS Windows Vista 64-bit, Matlab 7.14.0.739 (R2012a)

tic;
N=10;
M = 2*eye(N)-tril(ones(N,N));
M(:,N)=1;
x=ones(N,1);
y=M*x;
minresqlp(M,y);
toc;

% Enter minresqlp.  Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||
% n      =     10   ||b||    = 1.265e+01   shift    = 0.000e+00   rtol     = 2.220e-16
% maxit  =     40   maxxnorm = 1.000e+07   Acondlim = 1.000e+15   TranCond = 1.000e+07
% precon =      0
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
% P      0   1.26e+01   2.95e+01   1.00e+00   1.00e+00   0.00e+00   1.00e+00   0.00e+00
%        1   1.26e+01   3.37e+01   8.83e-01   3.22e-01   2.33e+00   1.00e+00   6.70e-01
%        2   1.20e+01   3.82e+01   5.61e-01   2.52e-01   8.34e+00   3.58e+00   1.04e+00
%        3   1.19e+01   4.21e+01   3.92e-01   1.99e-01   1.27e+01   5.43e+00   1.41e+00
%        4   1.17e+01   4.52e+01   2.81e-01   1.88e-01   1.77e+01   7.62e+00   1.63e+00
%        5   1.16e+01   4.88e+01   2.24e-01   1.48e-01   2.06e+01   8.84e+00   1.91e+00
%        6   1.15e+01   5.11e+01   1.61e-01   1.42e-01   2.84e+01   1.22e+01   2.07e+00
%        7   1.15e+01   5.46e+01   1.36e-01   1.21e-01   3.13e+01   1.35e+01   2.29e+00
%        8   1.14e+01   5.65e+01   1.06e-01   1.18e-01   3.93e+01   1.69e+01   2.42e+00
%        9   1.13e+01   5.98e+01   9.20e-02   1.05e-01   4.23e+01   1.81e+01   2.62e+00
%       10   1.13e+01   6.14e+01   7.55e-02   1.02e-01   5.02e+01   2.16e+01   2.72e+00
%       11   1.12e+01   6.46e+01   6.74e-02   9.39e-02   5.32e+01   2.28e+01   2.89e+00
%       12   1.12e+01   6.60e+01   5.73e-02   9.20e-02   6.12e+01   2.63e+01   2.98e+00
%       13   1.12e+01   6.91e+01   5.21e-02   8.58e-02   6.42e+01   2.76e+01   3.14e+00
%       14   1.11e+01   7.04e+01   4.54e-02   8.42e-02   7.22e+01   3.10e+01   3.22e+00
%       15   1.11e+01   7.33e+01   4.18e-02   7.95e-02   7.52e+01   3.23e+01   3.36e+00
%       16   1.10e+01   7.45e+01   3.71e-02   7.82e-02   8.32e+01   3.57e+01   3.43e+00
%       17   1.10e+01   7.74e+01   3.45e-02   7.45e-02   8.62e+01   3.70e+01   3.57e+00
%       18   1.10e+01   7.84e+01   3.10e-02   7.34e-02   9.42e+01   4.04e+01   3.63e+00
%       19   1.10e+01   8.13e+01   2.91e-02   7.04e-02   9.72e+01   4.17e+01   3.76e+00
%       20   1.09e+01   8.22e+01   2.65e-02   6.95e-02   1.05e+02   4.52e+01   3.81e+00
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
%       21   1.09e+01   8.50e+01   2.49e-02   6.69e-02   1.08e+02   4.64e+01   3.93e+00
%       22   1.09e+01   8.59e+01   2.29e-02   6.61e-02   1.16e+02   4.99e+01   3.98e+00
%       23   1.09e+01   8.86e+01   2.17e-02   6.40e-02   1.19e+02   5.12e+01   4.10e+00
%       24   1.09e+01   8.94e+01   2.01e-02   6.33e-02   1.27e+02   5.46e+01   4.14e+00
%       25   1.08e+01   9.21e+01   1.91e-02   6.14e-02   1.30e+02   5.59e+01   4.25e+00
%       26   1.08e+01   9.29e+01   1.78e-02   6.08e-02   1.38e+02   5.93e+01   4.30e+00
%       27   1.08e+01   9.54e+01   1.70e-02   5.92e-02   1.41e+02   6.06e+01   4.40e+00
%       28   1.08e+01   9.62e+01   1.60e-02   5.86e-02   1.49e+02   6.41e+01   4.44e+00
%       29   1.08e+01   9.87e+01   1.53e-02   5.72e-02   1.52e+02   6.54e+01   4.54e+00
%       30   1.08e+01   9.95e+01   1.44e-02   5.67e-02   1.60e+02   6.88e+01   4.58e+00
%       31   1.07e+01   1.02e+02   1.39e-02   5.54e-02   1.63e+02   7.01e+01   4.67e+00
%       32   1.07e+01   1.03e+02   1.31e-02   5.49e-02   1.71e+02   7.35e+01   4.71e+00
%       33   1.07e+01   1.05e+02   1.26e-02   5.38e-02   1.74e+02   7.48e+01   4.80e+00
%       34   1.07e+01   1.06e+02   1.20e-02   5.34e-02   1.82e+02   7.83e+01   4.83e+00
%       35   1.07e+01   1.08e+02   1.16e-02   5.24e-02   1.85e+02   7.96e+01   4.92e+00
%       36   1.07e+01   1.09e+02   1.10e-02   5.19e-02   1.93e+02   8.30e+01   4.95e+00
%       37   1.07e+01   1.11e+02   1.06e-02   5.10e-02   1.96e+02   8.43e+01   5.04e+00
%       38   1.06e+01   1.12e+02   1.01e-02   5.07e-02   2.04e+02   8.77e+01   5.07e+00
%       39   1.06e+01   1.14e+02   9.83e-03   4.98e-02   2.07e+02   8.90e+01   5.15e+00
%       40   3.09e+01   1.15e+02D  1.65e-02   1.72e-02D  2.15e+02   9.25e+01   8.61e+00
% 
% 
% Exit minresqlp.   flag  =      8   The iteration limit was reached                        
% Exit minresqlp.   iter  =     40   (MINRES      0, MINRES-QLP     40)
% Exit minresqlp.   rnorm =  3.0883e+01     rnorm  direct =  3.0883e+01
% Exit minresqlp.                           Arnorm direct =  1.1464e+02
% Exit minresqlp.   xnorm =  8.6130e+00     xnorm  direct =  8.6130e+00
% Exit minresqlp.   Anorm =  2.1536e+02     Acond         =  9.2463e+01
% Elapsed time is 0.042093 seconds.

% If we are in a scenario were we were not able to find a solution, it
% would be interesting to inform growth factor or condition number
[L,U]=lu(M);
maxU = max(max(abs(U))); maxM = max(max(abs(M)));
growthFactor = maxU/maxM;
condNumber = cond(M);
fprintf('\nGrowth factor is %g and condition number is %g\n\n', growthFactor,condNumber);
%It might help the user understand the situation.
