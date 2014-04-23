% This script is to estimate Bernoulli random variable p using Central Limit Theorem
function [p_clt,out_param]=meanMCBernoulli(Yrand,in_param)
out_param.n_clt = ceil(norminv(1-in_param.alpha/2)/(4*in_param.abstol^2));
out_param.npcmax = 1e6;
out_param.nmax = 1e9;
out_param.n_clt = min(out_param.nmax,out_param.n_clt);
p_clt = SplitColumnMean(Yrand,out_param.n_clt,out_param.npcmax);
end

function p =SplitColumnMean(RV,n,npcmax)
%%  Split The Param.n into columns
nopt=min(npcmax,n);
% numbers of samples per loop step
nn=floor(n/nopt); % number of loop steps
nremain=n-nn*nopt;
% number of samples in last loop step
nloop=repmat(nopt,1,nn);
%vector of numbers of samples per loop step
if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
sumY=0;
for iloop=1:nn %loops to save memory
    sumY=sumY+sum(RV(nloop(iloop)));
end
%%  Estimate p
p=sumY/n; %calculate the mean
end