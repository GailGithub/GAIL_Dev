%Run TestcubMC on the step function
clear all, close all
format compact

fun.funtype='geomean';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.03;
param.measure='normal';
param.impyes=false;
param.mmax = 20;
param.abstol=1*1e-2;
param.reltol=0;
param.transform = 'Baker';
param.method = {'cubLattice'};

test.nrep=500;
test.howoftenrep=10;
dimchoice=[1 2 4 8 16 32 64]';
ndim=size(dimchoice,1);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));
sigmin=0.1;
sigmax=0.7;
%test.randch.sigoverall=sigmin*(sigmax/sigmin).^rand(test.nrep,1);
test.randch.sigoverall=sigmin+(sigmax-sigmin).*rand(test.nrep,1);
test.randchoicefun=@randomchoiceGeo;
test.whichsample= param.method;
TestSettings(test,fun,param);
