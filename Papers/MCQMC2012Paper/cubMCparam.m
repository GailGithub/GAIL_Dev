function [interval,param]=cubMCparam(interval,param,whfield)
%   This function checks the inputs of cubMCparam to insure that 
%   they make sense, and that all of the necessary parameters are set
%
%   interval = integration interval
%   param = parameter structure
%   whfield = the names of the fields in the param structure to be checked
%       if whfield='fun', then only the function-related parameters are checked

param.exit=0; %success! until found otherwise

%Default values
 measureDefault='uniform'; %distribution to integrate against
     tolDefault=1e-3;      %absolute error tolerance
  sampleDefault='iid';     %sampling scheme
 errmethDefault='rep';     %error estimation method
   ndmaxDefault=1e9;       %maximum number of scalar sample points
 ndpcmaxDefault=1e5;       %maximum size of a vector of samples
      n0Default=2^13;      %initial number of samples
   n0MinDefault=30;        %minimum allowed initial number of samples
   alphaDefault=0.01;      %uncertainty
   npartDefault=8;         %number of parts for quasi-standard error
npartMinDefault=4;         %minimum allowed
npartMaxDefault=32;        %maximum allowed
   fudgeDefault=1.2;       %fudge factor for error estimate
impscaleDefault=1;         %importance sampling scale default
impshiftDefault=0;         %importance sampling shift default

if nargin<3; whfield=''; end %check everything

%Check the integration interval
[two, param.dim]=size(interval); %interval should be 2 x dimension
if two==0 && isfield(param,'interval'); %if interval specified through param structure
    interval=param.interval; %then get it from there
    [two, param.dim]=size(interval); %and get the dimension
end
if any(isnan(interval(:))); %check interval for not a number
    param.exit=10; param=cubMCerr(param); return; 
end
if two~=2 %if interval is given as row vector for dimension 1, fix that
    if param.dim==2; param.dim=two; interval=interval'; 
    else param.exit=11; param=cubMCerr(param); return; %else return an error
    end
end
interval=[min(interval,[],1); max(interval,[],1)]; %ensure left and right endpoints are in order
if any(interval(1,:)==interval(2,:)); %interval is a point in one direction
    param.exit=12; param=cubMCerr(param); return;
end
param.interval=interval; %copy interval into the param structure

%Check the parameter values

%Probability density to be integrated against
if isfield(param,'measure');
    param.measure=validatestring(param.measure,{'uniform','normal','Gaussian'});
    if strcmp(param.measure,'Gaussian'); param.measure='normal'; end
else
    param.measure=measureDefault;
end
if strcmp(param.measure,'uniform')&&~all(isfinite(interval(:)))
    %cannot integrate on an infinite interval with the uniform distribution
    param.exit=13; param=cubMCerr(param); return;
end
if strcmp(param.measure,'normal')&&any(isfinite(interval(:)))
    %must integrate on an infinite interval with the normal distribution
    param.exit=14; param=cubMCerr(param); return;
end
if strcmp(whfield,'fun'); return; end

%Error tolerance
if isfield(param,'tol');
    param.tol=abs(param.tol);
else
    param.tol=tolDefault; 
end

%Sampling scheme
if isfield(param,'sample');
    param.sample=validatestring(param.sample,{'iid','Sobol','net'});
    if strcmp(param.sample,'net'); param.sample='Sobol'; end
else
    param.sample=sampleDefault; 
end

if strcmp(param.sample,'Sobol')
    if ~isfield(param,'scramble'); param.scramble=false; end
end

%Method of estimating errors
if isfield(param,'errmeth');
    param.errmeth=validatestring(param.errmeth,{'rep','qse'});
else
    param.errmeth=errmethDefault; 
end

%Maxinum number of scalar values of x, i.e. number of samples x dimension
if isfield(param,'ndmax');
    param.ndmax=abs(param.ndmax);
else
    param.ndmax=ndmaxDefault;
end

%Maxinum number of scalar values of x per vector
if isfield(param,'ndpcmax');
    param.ndpcmax=abs(param.ndpcmax);
else
    param.ndpcmax=ndpcmaxDefault;
end

%Initial sample size
if isfield(param,'n0');
    param.n0=max(ceil(min(abs(param.n0),param.ndmax/param.dim)),n0MinDefault);
else
    param.n0=n0Default; 
end

%Uncertainty
if isfield(param,'alpha');
    param.alpha=min(max(abs(param.alpha),0),1);
else
    param.alpha=alphaDefault; 
end

%Number of parts to split sample into for quasi-standard error
if isfield(param,'npart');
    param.npart=ceil(max(min(abs(param.npart),npartMaxDefault),npartMinDefault));
else
    param.npart=npartDefault;
end

%Fudge factor to multiply error estimate by
if isfield(param,'fudge');
    param.fudge=max(abs(param.fudge));
else
    param.fudge=fudgeDefault; 
end

%Importance sampling parameters
if isfield(param,'impyes');
    if isfield(param,'impscale');
        if any(param.impscale<=0); %scale must be positive 
            param.exit=20; param=cubMCerr(param); return;
        end
        [one,dim]=size(param.impscale);
        if dim==param.dim
            param.impscale=param.impscale(1,:);
        else
            if one==param.dim; param.impscale=transpose(param.dim(:,1));
            else param.impscale=param.impscale(1,1)*ones(1,param.dim);
            end
        end
    else
        param.impscale=impscaleDefault*ones(1,param.dim); %default value
    end
    if isfield(param,'impshift');
        [one,dim]=size(param.impshift);
        if dim==param.dim
            param.impshift=param.impshift(1,:);
        else
            if one==param.dim; param.impshift=transpose(param.dim(:,1));
            else param.impshift=param.impshift(1,1)*ones(1,param.dim);
            end
        end
    else
        param.impshift=impshiftDefault*ones(1,param.dim); %default value
    end
else
    param.impyes=false; 
end

end