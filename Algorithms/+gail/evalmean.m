function p = evalmean(RV,n,npcmax)
%%  Split n samples into columns and evaluate the mean recursively
% RV --- the function to generate the random variables
% n --- the number of samples
% npcmax --- the maximum samples per loop
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