function [exactci,cltci]=binomialCI(ntot,nsuc,alpha)
%Author: Todor Markov

if nargin==0; help binomialCI, return, end %forgot to give inputs
obsprob=nsuc/ntot; %observed probability of success
if nargin < 3 %if no alpha input
   alpha = 0.01; %this is the default
end
al2=alpha/2; %half of uncertainty
if nsuc==0; %no successes observed
   plo=0; %the lower bound must be zero
else %use a nonlinear equation solver
   plo=fzero(@(x) binocdf(nsuc-1,ntot,x)-1+al2,[0,obsprob]);
end
if nsuc==ntot; %no failures observed 
   pup=1; %the upper bound must be one
else %use a nonlinear equation solver
   pup=fzero(@(x) binocdf(nsuc,ntot,x)-al2,[obsprob,1]); %nonlinear equation
end
exactci=[plo,pup]; %confidence interval based on exact probabilities
if nargout>1 %get the CLT one also
   cltci=obsprob+ ... %CLT confidence interval
   [-1 1]*(-norminv(al2)*sqrt(obsprob*(1-obsprob))/sqrt(ntot));
end

