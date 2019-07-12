function b=ispositive(a)
%
% ISPOSITIVE To judge if a variable is positive or not
b = isnumeric(a) && (a>0);