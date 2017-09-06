function [testfun,param]=geomMeanAsianCall(fun,param)
%   This function chooses and sets up a test function from the parameters
%      input by the user and contained in the structures
%      fun and param
%   fun.funtype         = type of test function
%   param.interval      = domain of test function
%   param.dim           = dimension of the domain
%   param.measure           = probability density function for integration
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

[~,param]=cubMCparam([],param,'fun'); %check validity of some parameters

if strcmp(fun.funtype,'geomean') %exponential test function
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
%Preliminary calculations
log2d=log2(fun.d);
rT=fun.r*fun.T;
rTfac=((fun.r-fun.sigma^2/2)*fun.T/fun.d)*(1:fun.d);
dov2l=fun.d;
ramp=zeros(log2d+1,fun.d);
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

testfun=@(x) geomeancall(x,fun.S0,fun.K,rT,fun.d,log2d,rTfac,ramp,dov2lmin1);

%% Compute exact integral of this function
rTmod=(fun.r-fun.sigma^2/2)*fun.T/2;
sigrootTmod=fun.sigma*sqrt(fun.T*(1-1/(2*fun.d)^2)/3);
x0=(log(fun.K/fun.S0)-rTmod)/sigrootTmod;
param.exactintegral=exp(-fun.r*fun.T)*(fun.S0*exp(rTmod+sigrootTmod^2/2)...
    *normcdf(sigrootTmod-x0) - fun.K*normcdf(-x0));

%% Prepare text description of the exponential function
paramendstring=char(10);
param.funDescribe=['   f(x) = ' ...
    'payoff of a geometric mean Asian call option with' paramendstring ...
    '       initial stock price = ' num2str(fun.S0) paramendstring ...
    '              strike price = ' num2str(fun.K) paramendstring ...
    '            time to expiry = ' num2str(fun.T) paramendstring ...
    '   # of monitoring periods = ' int2str(fun.d) paramendstring ...
    '             interest rate = ' num2str(fun.r) paramendstring ...
    '                volatility = ' num2str(fun.sigma) paramendstring];
end

function payoff=geomeancall(x,S0,K,rT,d,log2d,rTfac,ramp,dov2lmin1)
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
    

