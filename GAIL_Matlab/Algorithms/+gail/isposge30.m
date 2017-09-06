function b=isposge30(a)
%
% ISPOSGE3 To judge if a variable is greater than or equal to 30
b = isnumeric(a) && (a>=30) && (ceil(a)==a);