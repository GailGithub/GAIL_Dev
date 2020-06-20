%Run TestcubMC on the step function
%clear all, close all
function [res,test,fun,param] = RunTestcubMConGeoAsianCall()
format compact

fun.funtype='geomean';
fun.S0=100;
fun.K=100;
fun.T=1;
fun.r=0.03;
param.measure='Gaussian';
param.impyes=false;
param.abstol=2e-2;
param.n0=1024;

%test.nrep=500;%in the paper, we use 500 repilcation numbers
test.nrep = 50;
test.howoftenrep=10;
dimchoice=[1 2 4 8 16 32 64]';
ndim=size(dimchoice,1);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));
sigmin=0.1;
sigmax=3;
test.randch.sigoverall=sigmin*(sigmax/sigmin).^rand(test.nrep,1);
test.randchoicefun=@randchoiceGeoCall;
%test.whichsample={'iid','iidheavy','Sobol'};
test.whichsample={'iid','iidheavy','Sobol','Sobolheavy'};
res = TestcubMCDiffSettings(test,fun,param);
end
