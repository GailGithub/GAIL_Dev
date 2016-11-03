function f_eval = producttestfun(x,center)
f_eval = ones(size(x,1),1);
for i = 1:size(center,2)
f_eval = f_eval.*(x(:,i).^2+center(i));
end
end