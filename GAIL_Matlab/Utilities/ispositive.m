function b=ispositive(a)
%
% ISPOSITIVE To judge a variable is positive or not
b = isnumeric(a) && (a>=0);