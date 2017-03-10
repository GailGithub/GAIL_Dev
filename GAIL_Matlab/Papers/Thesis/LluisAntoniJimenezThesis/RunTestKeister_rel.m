%Run TestcubMC on the Keister integrand
clear all, close all
format compact

fun.funtype='Keister';
param.measure='uniform';
param.impyes=false;
param.mmax = 20;
param.abstol=2*1e-3;
param.reltol=2*1e-3;
param.transform = 'id';
param.method = {'cubLattice'};
param.toltype = 'comb';
param.theta = 0.0;

test.nrep=500;
test.howoftenrep=10;
test.randch.a=sqrt(2);
test.randch.dim=floor(20.^rand(test.nrep,1)); %random dimensions to 20
test.randchoicefun=@randomchoiceKeister;
test.whichsample= param.method;
TestSettings_rel(test,fun,param);