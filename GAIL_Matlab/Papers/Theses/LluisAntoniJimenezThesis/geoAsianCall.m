function [testfun,param]=geoAsianCall(fun,param)
%   This function chooses and sets up a test function from the parameters
%      input by the user and contained in the structures
%      fun and param
%   fun.funtype         = type of test function
%   param.interval      = domain of test function
%   param.dim           = dimension of the domain
%   param.measure       = probability density function for integration
%   param.exactintegral = exact value of the integral (scalar)
%   fun.shape           = shape parameter (1 x param.dim)
%   fun.scale           = scale parameter (1 x param.dim)
%   fun.addc            = additive constant (1 x param.dim)
%   fun.overaddc        = overall additive constant (scalar)
%   fun.overmultc       = overall multiplicative constant (scalar)

if nargin < 2; %give the basic default parameters 
    param.interval=[0;1]; %default integration interval
    if nargin < 1; fun.funtype='exp'; end %exponential test function
end
if ~isfield(param,'interval'); param.interval=[0;1]; end %default interval

%[~,param]=cubMCparam([],param,'fun'); %check validity of some parameters
%keyboard

if strcmp(fun.funtype,'geomean') %geometric mean Asian call test function
    [testfun,param]=makeGeometricMeanTestFun(fun,param);
else
    error('Function type not recognized')
end
end

%% Geometric Mean Call Integrand
function [testfun,param]=makeGeometricMeanTestFun(fun,param)
%Create the exponential test function
fun=verifyparam(fun,{'S0','K','T','d','r','sigma'}, ...
    {[1 1],[1 1],[1 1],[1 1],[1 1],[1 1]}, ...
    {100,100,1,1,0.03,0.5});
param.fun=fun; %copy of the function parameters

%if false %Force bb
if strcmp(param.sample,'iid') %Time discretization calculations
   rT=fun.r*fun.T;
   rTfac=((fun.r-fun.sigma^2/2)*fun.T/param.dim)*(1:param.dim);
   Delta=fun.sigma*sqrt(fun.T/param.dim);
   %keyboard
   testfun=@(x) geomeancalltd(x,fun.S0,fun.K,rT,param.dim,rTfac,Delta);
else %Brownian bridge calculations for qMC
   log2d=log2(param.dim);
   rT=fun.r*fun.T;
   rTfac=((fun.r-fun.sigma^2/2)*fun.T/param.dim)*(1:param.dim);
   dov2l=param.dim;
   ramp=zeros(log2d+1,param.dim);
   ramp(1,:)=(fun.sigma*sqrt(fun.T)/dov2l)*(1:dov2l);
   twolmin1=1;
   Deltal=fun.T/4;
   dov2lmin1=zeros(log2d,1);
   for l=1:log2d
       dov2lmin1(l)=dov2l;
       dov2l=dov2l/2;
       ramppc=(fun.sigma*sqrt(Deltal)/dov2l)*[1:dov2l dov2l-1:-1:0];
       ramp(l+1,:)=repmat(ramppc,1,twolmin1);
       Deltal=Deltal/2;
       twolmin1=twolmin1*2;
   end
   testfun=@(x) geomeancallbb(x,fun.S0,fun.K,rT,param.dim,log2d,rTfac,ramp,dov2lmin1);
end

%% Compute exact integral of this function
rTmod=(fun.r-fun.sigma^2/2)*fun.T/2;
sigrootTmod=fun.sigma*sqrt(fun.T*(1-1/(2*param.dim)^2)/3);
x0=(log(fun.K/fun.S0)-rTmod)/sigrootTmod;
param.exactintegral=exp(-fun.r*fun.T)*(fun.S0*exp(rTmod+sigrootTmod^2/2)...
    *normcdf(sigrootTmod-x0) - fun.K*normcdf(-x0));
end

function payoff=geomeancalltd(x,S0,K,rT,d,rTfac,Delta)
    n=size(x,1);
    %Create stock paths
    %keyboard
    Smatrix=repmat(rTfac,n,1)+cumsum(Delta*x,2);
    %Compute payoff
    loggmean=sum([Smatrix(:,1:d-1) Smatrix(:,d)/2],2)/d;
    %keyboard
    payoff=max(S0*exp(loggmean)-K,zeros(n,1))*exp(-rT);
end 

function payoff=geomeancallbb(x,S0,K,rT,d,log2d,rTfac,ramp,dov2lmin1)
    n=size(x,1);
    %Create stock paths
    Smatrix=repmat(rTfac,n,1);
    Smatrix=Smatrix+x(:,1)*ramp(1,:);
    twolmin1=1;
    for l=1:log2d
        jvec=twolmin1+(1:twolmin1);
        %keyboard
        weight=reshape(repmat(x(:,jvec),dov2lmin1(l),1),n,d);
        %keyboard
        Smatrix=Smatrix+weight.*repmat(ramp(l+1,:),n,1);
        twolmin1=twolmin1*2;
    end
    %Compute payoff
    loggmean=sum([Smatrix(:,1:d-1) Smatrix(:,d)/2],2)/d;
    %keyboard
    payoff=max(S0*exp(loggmean)-K,zeros(n,1))*exp(-rT);
end


