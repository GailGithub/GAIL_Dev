function b=isfcn(h)
%
%ISFCN To judge input is a function or not
b = isa(h, 'function_handle');