% Lan Jiang ljiang14@hawk.iit.edu
% MACI64 
% matlab version 8.1.0.604 (R2013a)
n = 10;
e = ones(n,1);
A = spdiags([e -2*e e], -1:1, n, n);
b = sum(A,2);
tstart = tic;
res = minresqlp(A,b);
time = toc(tstart);

% Enter minresqlp.  Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||
% n      =     10   ||b||    = 1.414e+00   shift    = 0.000e+00   rtol     = 2.220e-16
% maxit  =     40   maxxnorm = 1.000e+07   Acondlim = 1.000e+15   TranCond = 1.000e+07
% precon =      0
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
% P      0   1.41e+00   3.16e+00   1.00e+00   1.00e+00   0.00e+00   1.00e+00   0.00e+00
%        1   6.32e-01   1.02e+00   2.36e-01   5.56e-01   2.24e+00   1.00e+00   5.66e-01
%        2   3.78e-01   5.05e-01   9.65e-02   4.61e-01   2.90e+00   2.22e+00   8.63e-01
%        3   2.58e-01   3.02e-01   5.64e-02   4.03e-01   2.90e+00   3.32e+00   1.09e+00
%        4   1.91e-01   2.57e-02   3.71e-02   4.65e-02   2.90e+00   4.58e+00   1.28e+00
%        5   1.75e-15   5.25e-15D  1.65e-16   1.04e+00D  2.90e+00   3.55e+01   3.16e+00
% 
% 
% Exit minresqlp.   flag  =      1   A solution to Ax = b found, given rtol                 
% Exit minresqlp.   iter  =      5   (MINRES      0, MINRES-QLP      5)
% Exit minresqlp.   rnorm =  1.7484e-15     rnorm  direct =  1.7484e-15
% Exit minresqlp.                           Arnorm direct =  5.2545e-15
% Exit minresqlp.   xnorm =  3.1623e+00     xnorm  direct =  3.1623e+00
% Exit minresqlp.   Anorm =  2.8983e+00     Acond         =  3.5544e+01
%
%
% Suggestions:"EXAMPLE 4: Use this matrix-vector product function in file
% \Workouts\Afun.m:", Afun could not be found.