function b=isfcn(h)
%
%ISFCN To judge if input is a function handle or not
b = isa(h, 'function_handle');