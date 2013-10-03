% Xuan Zhou, xzhou23@hawk.iit.edu
%
% Windows 7 Ultimate Service Pack 1, Matlab 8.1.0.604 (R2013a)

n = 10;
RndMat = rand(n);
A = RndMat'*RndMat;
b = rand(n,1);
tic
x = minresqlp(A,b);
time = toc

% Enter minresqlp.  Min-length solution of symmetric (A-sI)x = b or min ||(A-sI)x - b||
% n      =     10   ||b||    = 1.719e+00   shift    = 0.000e+00   rtol     = 2.220e-16
% maxit  =     40   maxxnorm = 1.000e+07   Acondlim = 1.000e+15   TranCond = 1.000e+07
% precon =      0
% 
%     iter     rnorm     Arnorm   Compatible     LS        Anorm      Acond      xnorm
% P      0   1.72e+00   4.38e+01   1.00e+00   1.00e+00   0.00e+00   1.00e+00   0.00e+00
%        1   8.04e-01   7.57e-01   2.48e-01   3.24e-02   2.55e+01   1.00e+00   5.96e-02
%        2   4.37e-01   2.58e-01   1.71e-02   2.04e-02   2.90e+01   3.40e+01   8.19e-01
%        3   3.36e-01   2.49e-01   8.47e-03   2.56e-02   2.90e+01   6.20e+01   1.31e+00
%        4   2.83e-01   1.71e-01   5.44e-03   2.08e-02   2.90e+01   8.11e+01   1.73e+00
%        5   2.51e-01   1.28e-01   3.92e-03   1.76e-02   2.90e+01   1.12e+02   2.15e+00
%        6   2.20e-01   1.30e-02   2.70e-03   2.03e-03   2.90e+01   1.73e+02   2.76e+00
%        7   1.11e-01   1.30e-02   3.85e-04   4.05e-03   2.90e+01   1.16e+03   9.84e+00
%        8   9.23e-02   2.16e-02   2.91e-04   8.06e-03   2.90e+01   5.48e+02   1.09e+01
%        9   7.47e-02   9.45e-01   2.16e-04   4.36e-01   2.90e+01   1.16e+03   1.18e+01
%       10   6.24e-02   3.65e-03   1.72e-04   2.02e-03   2.90e+01   5.48e+02   1.24e+01
%       11   1.35e-11   2.56e-11   3.17e-14   6.54e-02   2.90e+01   1.42e+03   1.46e+01
%       12   4.60e-12   8.84e-12   1.08e-14   6.61e-02   2.90e+01   1.42e+03   1.46e+01
%       13   1.14e-12   1.37e-11   2.67e-15   4.15e-01   2.90e+01   1.42e+03   1.46e+01
%       14   9.67e-13   3.96e-12   2.27e-15   1.41e-01   2.90e+01   1.42e+03   1.46e+01
%       15   5.22e-14   4.66e-13D  1.22e-16   3.07e-01D  2.90e+01   1.42e+03   1.46e+01
% 
% 
% Exit minresqlp.   flag  =      1   A solution to Ax = b found, given rtol                 
% Exit minresqlp.   iter  =     15   (MINRES      0, MINRES-QLP     15)
% Exit minresqlp.   rnorm =  5.2162e-14     rnorm  direct =  5.2162e-14
% Exit minresqlp.                           Arnorm direct =  4.6550e-13
% Exit minresqlp.   xnorm =  1.4604e+01     xnorm  direct =  1.4604e+01
% Exit minresqlp.   Anorm =  2.9043e+01     Acond         =  1.4220e+03
% 
% time =
% 
%     1.0835



%Suggestion for improvement:  shminresqlp.m seems to be indentical to minresqlp.m.  