% RunTestcubMConMVNP generates the product test function in Section 4.5.3,
% and produces Figure 4.4 
clear all;close all;clc;
tstart=tic;
param.funtype = 'MVNP';
param.nrep = 500;
tempinitial = nan(param.nrep,1);
res.iidQ = tempinitial;
res.LatticeQ = tempinitial;
res.SobolQ = tempinitial;
res.iidtime=tempinitial;
res.Soboltime = tempinitial;
res.Latticetime = tempinitial;
res.exactQcubMC = tempinitial;
res.exactQ = tempinitial;
res.iiderr= tempinitial;
res.Latticeerr= tempinitial;
res.Sobolerr= tempinitial;
res.tol = tempinitial;
res.iidexit= tempinitial;
res.Sobolexit= tempinitial;
res.Latticeexit= tempinitial;
param.measure = 'uniform';
param.abstol = 1e-4;
param.reltol = 1e-3;
param.mmin = 13;
param.nSig=2^13;
param.nbudget = 1e10;

for i = 1:param.nrep
    if round(i/5)==i/5,i, end
    param.dima=2;
    param.dimb=20;
param.dim=randi([param.dima,param.dimb],1);
hyperbox = [-Inf*ones(1,param.dim) ; sqrt(param.dim)*rand(1,param.dim)];
sig = rand;
cov = sig*ones(param.dim,param.dim);
cov(1:param.dim+1:param.dim*param.dim) = 1;
out = MultivarNorProb(hyperbox,0,cov,param);
exactQ = integral(@(t)MVNPexact(t,hyperbox(2,:),sig),...
-inf, inf,'Abstol',1e-8,'RelTol',1e-8)/sqrt(2*pi);
res.exactQ(i) = exactQ;
res.iidQ(i) = out.iidQ;
res.iiderr(i) = abs(exactQ-out.iidQ);
res.iidtime(i)=out.iidparam.time;
res.iidexit(i) = out.iidparam.exit;
res.LatticeQ(i) = out.LatticeQ;
res.Latticetime(i)=out.Latticeparam.time;
res.Latticeerr(i) = abs(exactQ-out.LatticeQ);
res.Latticeexit(i) = out.Latticeparam.exitflag(1);
res.SobolQ(i) = out.SobolQ;
res.Soboltime(i)=out.Sobolparam.time;
res.Sobolerr(i) = abs(exactQ-out.SobolQ);
res.Sobolexit(i) = out.Sobolparam.exitflag(1);

end
res.kurtmax = out.iidparam.kurtmax;
res.tol=max(param.abstol,param.reltol*abs(res.exactQ));
subdir = 'LanThesisOutput';
filename = ['TestcubMCgFixedCovon-' param.funtype '-N'...
    int2str(param.nrep) 'd' int2str(param.dim)  ...
    'abstol' num2str(param.abstol) 'rel' num2str(param.reltol)];
gail.save_mat(subdir,filename,true,res,param);
tottime = toc(tstart)
%[exact_p , e_out] = cubSobol_g(@(t) prod(gail.stdnormcdf(bsxfun(@plus, hyperbox(2,:),sqrt(sig)*t)/sqrt(1-sig)),2),[-Inf;Inf],'normal',1e-9,0);
%disp([num2str(approx_p) ' in ' num2str(a_out.time) ' seconds.'])
% disp([num2str(exact_p) ' in ' num2str(e_out.time) ' seconds.'])
% disp(['Error is ' num2str(abs(approx_p-exact_p))])

