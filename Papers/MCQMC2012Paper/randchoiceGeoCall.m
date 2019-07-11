function [testfun,fun,param]=randchoiceGeoCall(fun,param,rchparam,irep)

param.dim=rchparam.dim(irep);
param.interval=[-Inf(1,param.dim); Inf(1,param.dim)];
[testfun,param]=geomMeanAsianCall(fun,param);

