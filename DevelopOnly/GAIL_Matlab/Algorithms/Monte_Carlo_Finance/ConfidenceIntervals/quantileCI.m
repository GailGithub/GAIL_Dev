oofunction quantci=quantileCI(quant,Xsample,extremes,alpha)
% quantci=QUANTILECI(quant,Xsample,alpha,extremes) 
%   computes  1-alpha  confidence intervals 
%   for the  quant  quantile of a random variable  X 
%   with extreme values  extremes
%   using IID data  Xsample
if nargin==0; help quantileci, return %forgot to give inputs
elseif nargin < 4 %if no alpha input
   alpha = 0.01; %this is the default
   if nargin < 3
      extremes = [-Inf Inf];
   end
end
n=length(Xsample); %number of samples
Xorder=[extremes(1); sort(Xsample); extremes(2)]; %order statistics
al2=alpha/2; %half significance level
lo=1+binoinv_bs_ver2(al2,n,quant); %position of the lower bound
up=2+binoinv_bs_ver2(1-al2,n,quant); %position of the upper bound
   %using binoinv_bs_ver2 because binoinv is too slow
quantci=[Xorder(lo),Xorder(up)]; %confidence interval for quantile