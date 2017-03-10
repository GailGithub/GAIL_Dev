% randchoiceproduct - randomly generates the parameters in product test
% functions in Section 3.5.1 and 4.5.1
function [testfun,fun,param]=randchoiceproduct(fun,param,rchparam,irep)
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
[testfun,param]=choosetestfun(fun,param);