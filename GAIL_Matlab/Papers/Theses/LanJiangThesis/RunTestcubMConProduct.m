% RunTestcubMConProduct generates the product test function in Section
% 3.5.1 and produces Figure 3.4
clear all;close all;clc;
tstart=tic;
param.funtype = 'product';
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
res.kurtvec = tempinitial;
param.measure = 'uniform';
param.abstol = 1e-3;
param.reltol = 1e-3;
param.mmin = 13;
param.nSig=2^13;
param.nbudget = 1e10;

for i = 1:param.nrep
    if round(i/5)==i/5,i, end
    param.dima=2;
    param.dimb=20;
param.dim=randi([param.dima,param.dimb],1);
%param.dim = 2;
hyperbox = [zeros(1,param.dim) ; ones(1,param.dim)];
center = rand(1,param.dim).*(4/3);

out=productfun(center,hyperbox,param);

exactQ = prod(1/3+center);
exactQ2 = prod(1/5+2/3*center+center.^2);
exactQ3 = prod(1/7+3/5*center+center.^2+center.^3);
exactQ4 = prod(1/9+4/7*center+6/5*center.^2+4/3*center.^3+center.^4);
exactvar = exactQ2-exactQ.^2;
exactkurt = (exactQ4-4*exactQ3.*exactQ+6*exactQ2.*exactQ.^2-3*exactQ.^4)./exactvar.^2;
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
res.kurtmax = out.iidparam.kurtmax;
res.kurtvec(i) = exactkurt;
end
res.tol=max(param.abstol,param.reltol*abs(res.exactQ));
subdir = 'LanThesisOutput/ProductFun';
filename = ['TestcubMCgFixedCovon-' param.funtype '-N'...
    int2str(param.nrep) 'd' int2str(param.dim)  ...
    'abstol' num2str(param.abstol) 'rel' num2str(param.reltol)];
gail.save_mat(subdir,filename,true,res,param);
tottime = toc(tstart)
%[exact_p , e_out] = cubSobol_g(@(t) prod(gail.stdnormcdf(bsxfun(@plus, hyperbox(2,:),sqrt(sig)*t)/sqrt(1-sig)),2),[-Inf;Inf],'normal',1e-9,0);
%disp([num2str(approx_p) ' in ' num2str(a_out.time) ' seconds.'])
% disp([num2str(exact_p) ' in ' num2str(e_out.time) ' seconds.'])
% disp(['Error is ' num2str(abs(approx_p-exact_p))])

