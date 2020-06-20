function [Q,param,x]=cubMC(f,interval,param)
%   CUBMC evaluates a multidimensional integral 
%     by (quasi-)Monte Carlo integration
%     to a specified tolerance (under some assumptions)
%
%          f = the function handle of the integrand
%   interval = the 2 x d multidimensional interval (box) of integration
%      param = a structure with various parameters used in the calculation
%      param.tol      = the error tolerance
%      param.measure  = the name of a measure 
%                         uniform (default) or normal
%                         if normal, the interval is R^d
%      param.sample   = the kind of sampling scheme 
%                         iid (default), lattice, or Sobol
%      param.errmeth  = the kind of error estimation
%                         replications (default) or quasi-standard error
%      param.ndmax    = maximum number of function evaluations
%      param.ndpcmax  = number of elements in an array of optimal size
%      param.n0       = the initial sample size
%      param.fudge    = the fudge factor used to estimate the error
%      param.impyes   = perform importance sampling 
%                         (only for normal measure)
%      param.impscale = standard deviation for importance sampling
%      param.impshift = mean shift for importance sampling
%      param.exit     = the state of program when exiting
%      param.fun      = parameters that define the test function (optional)
%
%   The integral to be evaluated takes the form
%       ?   f(x) rho(x) dx
%      R^d
%   where f(x) is the integrand
%         rho(x) is the Radon-Nikodym derivative of the measure
%                = 1 on a finite interval for 'uniform' or
%                = (2*pi)^(-d/2) exp(-x'*x/2) for 'normal'

tstart=tic; %start clock
if nargin<2; interval=[0,1]; end %default interval
if nargin<3; param.tol=1e-2; end %default tolerance
[interval,param]=cubMCparam(interval,param); %check validity of inputs
function newf=transformIntegrand(oldf,interval,param) 
        %function to transform the integrand
 %    so that the sample points don't have to be transformed

if strcmp(param.measure,'uniform') %uniform measure
    a=interval(1,:); %left endpoint
    b=interval(2,:); %right endpoint
    if all(a==0) && all(b==1) %no change needed
        newf=oldf; 
    else %transform points and integrand
        bmina=b-a; %interval width
        volbox=prod(bmina); %volume of the interval
        newf=@(x) oldf(x.*repmat(bmina,size(x,1),1)+repmat(a,size(x,1),1))...
            .*volbox;
        
       %stretch and shift, then multiply by volume
    end
elseif strcmp(param.measure,'normal')
    if strcmp(param.sample,'iid') %iid sampling
        if param.impyes
            if all(param.impscale==1) %no expansion or contraction 
                if all(param.impshift==0) %no change needed
                    newf=oldf; 
                else %just a shift
                    newf=@(x) oldf(x+repmat(param.impshift,size(x,1),1)) ...
                        .* exp(-x.*param.impshift'- ...
                        sum(param.impshift.*param.impshift)/2);
                end
            else %expansion or contraction needed
                %keyboard
                if all(param.impshift==0) %no shift
                    newf=@(x) cvfunscaleonly(x,oldf,param.impscale);
                else %need a shift
                    newf=@(x) cvfunscaleshift(x,oldf,...
                        param.impscale,param.impshift);
                end
            end
        else
            newf=oldf;
        end
    else %quasi-Monte Carlo sampling
        if param.impyes
            if all(param.impscale==1) %no expansion or contraction 
                if all(param.impshift==0) %no change needed
                    newf=@(x) oldf(norminv(x)); 
                else %just a shift
                    newf=@(x) cvfunqmcshiftonly(x,oldf,param.impshift);
                end
            else %expansion or contraction needed
                if all(param.impshift==0) %no shift
                    newf=@(x) cvfunscaleonly(norminv(x),oldf,param.impscale);
                else %need a shift
                    newf=@(x) cvfunscaleshift(norminv(x),oldf,...
                        param.impscale,param.impshift);
                end
            end
        else
            newf=@(x) oldf(norminv(x));
        end
    end
end
end
f=transformIntegrand(f,interval,param); 
%make transformations of the integrand so the sample points don't have to
%   be transformed

if strcmp(param.sample,'iid') %iid sampling
    tpcstart=tic; %start clock for sample generation
    if strcmp(param.measure,'uniform')
        x=rand(param.n0,param.dim);%uniform random numbers
        %save('xmat.mat','x')
    else 
        x=randn(param.n0,param.dim); %normal random numbers
    end
    tsample=toc(tpcstart); %time for sample generation
    tpcstart=tic; %start clock for integrand evaluation
    fx=f(x);%evaluate integrand
    %keyboard
    tintegrand=toc(tpcstart); %time for integrand evaluation
    sig0=std(fx); %sample standard deviation
    
    param.estvari=var(fx)./(prod(param.interval(2,:)-param.interval(1,:))).^2;
    sig0up=param.fudge*sig0; %upper bound on standard deviation
    %estimate sample size needed
    alpha1=1-sqrt(1-param.alpha);
    param.kurtmax=(param.n0-3)/(param.n0-1) ...
        + ((alpha1*param.n0)/(1-alpha1))*(1-1/param.fudge^2)^2;
    %sample size for estimating mean, no smaller than n0
    if sig0up==0;
        param.n=param.n0;
    else
       toloversig=param.tol/sig0up;
       ncheb=ceil(1/(alpha1*toloversig.^2));
       
       A=18.1139;
       A1=0.3322;
       A2=0.429;
       A3= 0.3031;
       A4= 0.646;
       A5= 0.469; % Six constants in Berry-Esseen inequality
       M3upper=param.kurtmax^(3/4);%using Jensen inequality to
       % bound the third moment
       BEfun=@(logsqrtn)normcdf(-exp(logsqrtn).*toloversig)...
           +exp(-logsqrtn).*min([A1*(M3upper+A2),A3*(M3upper+A4),A5*M3upper, ...
           A*M3upper./(1+(exp(logsqrtn).*toloversig).^3)])- alpha1/2;
       
%        A=18.1139;
%        A1=0.3328;
%        A2=0.429; %constant in Berry-Esseen inequality
%        M3upper=param.kurtmax^(3/4);
%         BEfun=@(logsqrtn) ...
%                 normcdf(-exp(logsqrtn).*toloversig)...
%                 +exp(-logsqrtn).*min(A1*(M3upper+A2), ...
%                 A*M3upper./(1+(exp(logsqrtn).*toloversig).^3))...
%                 - alpha1/2;
            
        if BEfun(log(sqrt(param.n0))) <= 0 || ncheb <=param.n0;
            param.n=param.n0;
        else
            logsqrtnCLT=log(norminv(1-alpha1/2)/toloversig);
            param.n=min(ncheb,ceil(exp(2*fzero(BEfun,logsqrtnCLT))));
        end
    end
    nd=param.n*param.dim; %total number of scalars in sample
    if nd>param.ndmax %too many samples for the time allowed
        param.exit=1;
        cubMCerr(param,tstart); %print warning message
        param.n=floor(param.ndmax/param.dim);
    end
    %# of samples per loop step
    nopt=min(floor(param.ndpcmax/param.dim),param.n); 
    %Put all samples in one or more vectors
    nn=floor(param.n/nopt); %number of loop steps
    nremain=param.n-nn*nopt; %# of samples in last loop step
    nloop=repmat(nopt,1,nn); %vector of # of samples per loop step
    if nremain>0; nloop=[nloop nremain]; nn=nn+1; end
    sumfx=0;
    for iloop=1:nn %loops to save memory  
        tpcstart=tic; %start clock for sample generation
        if strcmp(param.measure,'uniform')
            x=rand(nloop(iloop),param.dim); %uniform random numbers
        else
            x=randn(nloop(iloop),param.dim); %normal random numbers
        end
        tsample=tsample+toc(tpcstart); %time for sample generation
        tpcstart=tic; %start clock for integrand evaluation
        sumfx=sumfx+sum(f(x)); %evaluate integrand and sum values
        tintegrand=tintegrand+toc(tpcstart); %time for integrand evaluation
    end
    param.Q=sumfx/param.n; %approximate integral is sample mean
    param.n=param.n+param.n0; %total number of samples used

elseif strcmp(param.sample,'Sobol') %Sobol' sampling
    stream=sobolset(param.dim); %create Sobol set
    if param.scramble %scramble it
        stream=scramble(stream,'MatousekAffineOwen'); %scrambled Sobol set
    end
    stream=qrandstream(stream); %prepare stream for the rand command
    log2npart=ceil(log2(param.npart)); %log # of pieces of sample
    param.npart=2^log2npart; %number of pieces of sample
    newnpart=param.npart; %number of new parts of the sample   
    param.n=max(2^ceil(log2(param.n0)),param.npart);
        %initial sample size, a power of 2
    newn=param.n; %number of samples to take next
    partLength=newn/param.npart; %number of samples per part
    notDone=true; %whether to keep iterating    
    nopt=2^floor(log2(param.ndpcmax/param.dim));
        %optimal number of samples per vector
    meanPart=zeros(1,param.npart); %vector of means of parts of sample
    nparti=0; %index to add new parts
    
    tsample=0;
    tintegrand=0;
    while notDone 
%        nd=newn*param.dim; %total number of scalars in sample
        if partLength<=nopt %fit at least one part in one vector
            npts=min(newn,nopt); %number of points per step
            nn=newn/npts; %number of steps in loop
            npartloop=npts/partLength; %number of parts in step
            for iloop=1:nn
                tpcstart=tic; %start clock for sample generation
                x=rand(stream,npts,param.dim); %sample next Sobol' points
                tsample=tsample+toc(tpcstart); %time for sample generation
                tpcstart=tic; %start clock for integrand evaluation
                fx=f(x); %evaluate integrand
                tintegrand=tintegrand+toc(tpcstart); 
                    %time for integrand evaluation
                fxPart=reshape(fx,partLength,npartloop); 
                    %sample means of parts
                meanPart(nparti+(1:npartloop))=mean(fxPart,1);
                nparti=nparti+npartloop; %increment part index
            end            
        else %fit one part in several vectors
%        elseif nd<=param.ndmax %fit one part in several vectors
%            keyboard
            nn=partLength/nopt; %number of steps per part
            for ipartloop=1:newnpart
                sumfx=0; %initialize sum
                nparti=nparti+1; %increment part part index
                for iloop=1:nn
                    tpcstart=tic; %start clock for sample generation
                    x=rand(stream,nopt,param.dim); %next Sobol' points
                    tsample=tsample+toc(tpcstart); %time for sample generation
                    tpcstart=tic; %start clock for integrand evaluation
                    fx=f(x); %evaluate integrand
                    tintegrand=tintegrand+toc(tpcstart); 
                        %time for integrand evaluation
                    sumfx=sumfx+sum(fx); %sum of integrand values
                end
                meanPart(nparti)=sumfx/partLength; %mean of part
            end
%            keyboard
        end
        qse=std(meanPart)/sqrt(param.npart); %quasi-standard error for mean
        errEst=param.fudge*qse; %corrected by a fudge factor
        if errEst>param.tol; %tolerance not yet satisfied
            newn=param.n; %number of samples to take next
            param.n=2*param.n; %double total sample size for next round
            nd=param.n*param.dim;
            if nd>param.ndmax %too many samples for the time allowed
                notDone=false; %so have to stop
                param.exit=1; %set warning flag
                cubMCerr(param,tstart); %print warning message
            else %prepare for next iteration
                partLength=partLength*2; %double the length of each part
                meanPart(1:param.npart/2)=(meanPart(1:2:param.npart-1)...
                    +meanPart(2:2:param.npart))/2; %update meanPart vector
                newnpart=param.npart/2; %need only half parts next
                nparti=newnpart; %initial index starts halfway in future
            end
        else %terminate
            notDone=false;
        end
    end
    param.Q=mean(meanPart); %compute approximate integral
end
Q=param.Q; %assign answer
param.time=toc(tstart); %elapsed time
param.tsample=tsample;
param.tintegrand=tintegrand;
end



function y=cvfunscaleonly(x,oldf,scale)
%function to perform expansion or contraction only for control variates
    n=size(x,1); %number of points
    %keyboard
    y=oldf(x.*repmat(scale,n,1)) ...
        .* exp(-sum(x.*x.*repmat(scale.*scale-1,n,1),2)/2) ...
        .* prod(scale);
end
        
function y=cvfunscaleshift(x,oldf,scale,shift)
%function to perform expansion or contraction + mean shift
%   for control variates
    n=size(x,1); %number of points
    y=oldf(x.*repmat(scale,n,1)+repmat(shift,n,1)) ...
        .* exp(-((sum(x.*x.*repmat(scale.*scale-1,n,1),2) ...
        + shift*shift'))/2 + x*transpose(shift.*scale)) ...
        .* prod(scale);
end

function y=cvfunqmcshiftonly(x,oldf,shift)
%function to perform mean shift only
%   for control variates with quasi-Monte Carlo sampling
    n=size(x,1); %number of points
    z=norminv(x);  %inverse normal cumulative distribution function
    y=oldf(z+repmat(shift,n,1)) ...
        .* exp(-(shift*shift'/2+z*shift')).* prod(scale).^2;
end

        



    
