%Run TestcubMC on the Keister integrand
clear all, close all
format compact

fun.funtype='Keister';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.03;
param.measure='uniform';
param.impyes=false;
param.n0=1024;
param.errtype='comb';
param.theta=1; % Pure absolute error

test.nrep=500;
test.howoftenrep=10;
test.randch.a=sqrt(2);
test.randch.dim=floor(20.^rand(test.nrep,1)); %random dimensions to 20
test.randchoicefun=@randchoiceKeister;
%test.whichsample={'iidheavy','cubSobol'};
%test.whichsample={'iid'}; param.abstol=2e-3;
test.whichsample={'cubLattice'}; param.abstol=2e-3;
%test.whichsample={'cubSobol','cubLattice'}; param.abstol=2e-3;
TestCubatureDiffSettings(test,fun,param);
