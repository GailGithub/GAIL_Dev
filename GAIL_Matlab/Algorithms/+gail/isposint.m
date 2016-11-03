function b=isposint(a)
%
% ISPOSINT To judge if input is a positive integer or not
b = (ceil(a)==a) && (a>0);