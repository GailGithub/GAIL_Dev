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
param.abstol=2e-2;
param.n_sigma=1e4;
param.errtype='max';
param.theta=1; % Pure absolute error

test.nrep=500;
test.howoftenrep=10;
dimchoice=[1 2 4 8 16 32 64]';
ndim=size(dimchoice,1);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));
sigmin=0.1;
sigmax=0.7;
%test.randch.sigoverall=sigmin*(sigmax/sigmin).^rand(test.nrep,1);
test.randch.sigoverall=sigmin+(sigmax-sigmin).*rand(test.nrep,1);
test.randchoicefun=@randchoiceGeoCall;
%test.whichsample={'iidheavy','cubSobol'};
test.whichsample={'cubLattice'};
%test.whichsample={'iid'};
TestCubatureDiffSettings(test,fun,param);
