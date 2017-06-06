% randchoiceKeister - randomly generates the parameters in Keister test
% functions in Section 4.5.2
function [testfun,fun,param]=randchoiceKeister(fun,param,rchparam,irep)
%Choose Keister function
%param.dim=rchparam.dim(irep);
param.dim=rchparam.dim(irep);
param.interval=[zeros(1,param.dim); ones(1,param.dim)];
param.exactintegral = Keistertrue(param.dim);
testfun = @(x) Keisterfun(x,param.dim); 
end

function f=Keisterfun(x,dim) %a must be bigger than sqrt(2)
f = cos(sqrt(sum(norminv(x(:,1:dim)).^2,2)/2))*pi^(dim/2);
end


