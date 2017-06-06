% randchoicegaussian randomly generates the parameters in test functions
% used by choosetestfun for multivariate integration cases.
function [testfun,fun,param]=randchoicegaussian(fun,param,rchparam,irep)

shapevec=rchparam.shapeoverall(irep,:);
scalevec=rchparam.scaleoverall(irep,:);
addcvec=rchparam.addcoverall(irep,:);
overmultcvec=rchparam.overmultcoverall(irep,:);
overaddcvec=rchparam.overaddcoverall(irep,:);
centercvec=rchparam.centercoverall(irep,:);
dimvec=rchparam.dimoverall(irep,:);
fun.shape=shapevec;
fun.scale=scalevec;
fun.addc=addcvec;
fun.center=centercvec;
fun.overmultc=overmultcvec;
fun.overaddc=overaddcvec;
param.dim=dimvec;
param.interval=[zeros(1,dimvec);ones(1,dimvec)];
param.bmina=param.interval(2,:)-param.interval(1,:);
%fun.shift=rand(1,param.dim);
%param.dim=rchparam.dim(irep);
[testfun,param]=choosetestfun(fun,param);