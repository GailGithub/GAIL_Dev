function [testfun,fun,param]=randomchoiceKeister(fun,param,rchparam,irep)
%Choose Keister function
param.a=rchparam.a;
param.dim=rchparam.dim(irep);
param.interval=[zeros(1,param.dim); ones(1,param.dim)];
param.exactintegral = Keistertrue_check(param.dim);
testfun = @(x) Keisterfun(x,param.dim,param.a)./param.exactintegral; 
   %normalize to unity
param.exactintegral=1;
%keyboard
end

function f=Keisterfun(x,d,a) %a must be bigger than sqrt(2)
sumsq=sum(norminv(x(:,1:d)).^2,2);
f=((sqrt(2*pi)/a).^d).*exp(-(1./(a.^2)-1./2)*sumsq).*cos(sqrt(sumsq)/a);
%keyboard
end


