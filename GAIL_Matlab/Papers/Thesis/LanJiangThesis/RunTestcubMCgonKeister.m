% RunTestcubMCgonKeister generates the Gaussian test function in Section
% 4.5.2 and produces Figure 4.3.

clear all, close all
format compact
test.nrep=500;
fun.funtype='Keister';
param.measure='uniform';
param.impyes=false;
param.abstol=1e-2;
param.reltol = 0;
param.nSig=2^13;
param.mmin = 13;
param.nbudget=1e10;

test.howoftenrep=5;
dimchoice=1:8;
ndim=size(dimchoice,2);
test.randch.dim=dimchoice(randi(ndim,test.nrep,1));

test.randchoicefun=@randchoiceKeister;
test.whichsample={'iid','cubLattice','cubSobol'}; 
%test.whichsample={'cubSobol'}; 
%test.whichsample={'cubSobol','cubLattice'}; param.abstol=2e-3;
%res=TestcubMCgDiffSettings(test,fun,param)
res=TestcubMCgRELTOLDiffSettings(test,fun,param);
