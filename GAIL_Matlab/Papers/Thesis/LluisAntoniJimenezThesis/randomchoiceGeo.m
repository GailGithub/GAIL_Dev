function [testfun,fun,out_param]=randomchoiceGeo(fun,in_param,rchparam,irep)
out_param=in_param;
out_param.dim=rchparam.dim(irep);
out_param.interval=[-Inf(1,out_param.dim); Inf(1,out_param.dim)];
%keyboard
[testfun,out_param]=geoAsianCall(fun,out_param);

