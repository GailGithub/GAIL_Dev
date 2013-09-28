function [testfun,fun,out_param]=randchoicegaussiand1(fun,in_param,rchparam,irep)

shapevec=rchparam.shapeoverall(irep,:);
scalevec=rchparam.scaleoverall(irep,:);
addcvec=rchparam.addcoverall(irep,:);
overmultcvec=rchparam.overmultcoverall(irep,:);
overaddcvec=rchparam.overaddcoverall(irep,:);
centercvec=rchparam.centercoverall(irep,:);
fun.shape=shapevec;
fun.scale=scalevec;
fun.addc=addcvec;
fun.center=centercvec;
fun.overmultc=overmultcvec;
fun.overaddc=overaddcvec;
%fun.shift=rand(1,param.dim);
%param.dim=rchparam.dim(irep);
[testfun,out_param]=choosetestfun(fun,in_param);