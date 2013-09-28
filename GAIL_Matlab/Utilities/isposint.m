function b=isposint(a)
%
% ISPOSINT To judge input is positive integer or not
b = (ceil(a)== a) && (a>=0);