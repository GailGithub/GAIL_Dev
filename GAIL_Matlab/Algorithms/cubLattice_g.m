function [q,out_param]=cubLattice_g(varargin)
% (f,d,abstol,density,shift,mmin, mmax, fudge,diff)
% cubLattice is a Quasi-Monte Carlo method to evaluate a multidimensional integral
% given a specified absolute error tolerance with guaranteed absolute
% error, using rank-1 lattices.
% mmin<mmax
% Example:
%   f=@(x) x(:,1).*x(:,2);
%   d=2;
%   abstol=10^-5;
%   trueint=1/4;
%   error=abs(trueint-cubLattice(f,d,abstol));
% 
% [q,err,time,n]=cubLattice(f,d,abstol,density,shift) estimates the integral of f over the
% d-dimensional unit cube using rank-1 Lattice rules...
% f --- the integrand...
% 
% Guarantee: see papers...
% 
% Examples
% 
% Example 2: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2), where x is a vector x = [x1 x2].
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2); d=2;
% >> Q = cubLattice(f,d,1e-3)
% Q = 0.55***
% 
% See also cubSobol, cubMC_g, meanMC_g, integral_g
% 
% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
% conservative fixed width confidence intervals via Monte Carlo sampling,
% Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
% Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
% arXiv:1208.4318 [math.ST]
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
% Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
% 1.3.0)" [MATLAB Software], 2014. Available from
% http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.

%%
% Check and initialize parameters
[f,out_param] = cubLattice_g_param(varargin{:});

tic
%% Initialize parameters   
if strcmp(out_param.density,'normal')
   f=@(x) f(norminv(x));
end
if strcmp(out_param.diff,'Baker')
    f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(out_param.diff,'C0')
    f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(out_param.diff,'C1')
    f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(out_param.diff,'C1sin')
    f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
end

%% Main algorithm
mlag=4;
Stilde=zeros(out_param.mmax-out_param.mmin+1,1);
errest=zeros(out_param.mmax-out_param.mmin+1,1);
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1);

%% Initial points and FWT
out_param.n=2^out_param.mmin;
xpts=mod(lattice_gen(1,out_param.n,out_param.d)+out_param.shift,1); n0=out_param.n;
y=f(xpts);
yval=y;

%% Compute initial FFT
for l=0:out_param.mmin-1
   nl=2^l;
   nmminlm1=2^(out_param.mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;
end % 'y' now contains the FWT coefficients

%% Approximate integral
q=mean(yval);
appxinteg(1)=q;

%% Create kappanumap
kappanumap=(1:out_param.n)'; %initialize map
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %
   temp=kappanumap(nl+1+flip);
   kappanumap(nl+1+flip)=kappanumap(1+flip);
   kappanumap(1+flip)=temp;
end

%% Compute Stilde
nllstart=int64(2^(out_param.mmin-mlag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
out_param.err=out_param.fudge*2^(-out_param.mmin)*Stilde(1);
errest(1)=out_param.err;
if out_param.err <= out_param.abstol; out_param.time=toc; return, end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax 
   out_param.n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=mod(lattice_gen(nnext+1,2*nnext,out_param.d)+out_param.shift,1);
   n0=n0+nnext;
   ynext=f(xnext);
   yval=[yval; ynext];

   %% Compute initial FFT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
      coefv=repmat(coef,nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+coefv.*oddval)/2;
      ynext(~ptind)=(evenval-coefv.*oddval)/2;
   end

   %% Compute FFT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)/(2*nl))';
   coefv=repmat(coef,nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+coefv.*oddval)/2;
   y(~ptind)=(evenval-coefv.*oddval)/2;

   %% Update kappanumap
   kappanumap=[kappanumap; (nnext+1:out_param.n)']; %initialize map
   for l=m-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone);
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end

   %% Compute Stilde
   nllstart=int64(2^(m-mlag-1));
   meff=m-out_param.mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   out_param.err=out_param.fudge*2^(-m)*Stilde(meff);
   errest(meff)=out_param.err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if out_param.err <= out_param.abstol; out_param.time=toc; break
    elseif m==out_param.mmax; warning('The maximum budget is attained without reaching the tolerance.');
   end

end
out_param.time=toc;
end





%% Parsing for the input of cubLattice_g
function [f, out_param] = cubLattice_g_param(varargin)

% Default parameter values
default.d  = 1;
default.abstol  = 1e-4;
default.density  = 'uniform';
default.shift  = rand;
default.mmin  = 10;
default.mmax  = 24;
default.fudge = 3;
default.diff = 'C1sin';

if isempty(varargin)
    help cubLattice_g
    warning('At least, function f must be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
else
    f = varargin{1};
    out_param.f=f;
end;

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin{3};
    validvarargin=(isnumeric(in3) || isstruct(in3) ...
        || ischar(in3));
end

if ~validvarargin
    % If only one input f or in2 is nothing above, use all the default parameters
    % warning(['MATLAB:cubLattice_g:mininputarg',' At least function and dimension must be specified and all parameters must be numeric or strings.'])
    out_param.d = default.d;
    out_param.abstol = default.abstol;
    out_param.density = default.density;
    out_param.shift = default.shift;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
    out_param.diff = default.diff;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in3)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'d',default.d,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'shift',default.shift,@isnumeric);
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'diff',default.diff,...
            @(x) any(validatestring(x, {'Baker','C0','C1','C1sin'})));
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'d',default.d,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addParamValue(p,'shift',default.shift,@isnumeric);
        addParamValue(p,'mmin',default.mmin,@isnumeric);
        addParamValue(p,'mmax',default.mmax,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'diff',default.diff,...
            @(x) any(validatestring(x, {'Baker','C0','C1','C1sin'})));
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% For dimension to be positive integer
if (~isposint(out_param.d)) % Dimension should be a postitive integer
    warning('MATLAB:cubLattice_g:dnotposint',...
        ['The dimension should be a positive integer,'...
        'We take the ceil of the the absolute value.'])
    out_param.d = ceil(abs(out_param.d));
end

% Force error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['MATLAB:cubLattice_g:abstolnonpos','Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

%%%% Check density

% Force mmin to be integer greater than 0
if (out_param.mmin < 1 || ~isposint(out_param.mmin))
    warning(['MATLAB:cubLattice_g:lowmmin',' The minimum starting exponent should be an integer greater or equal than 1.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater or equal than
% mmin an smaller than 28
if ~(isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=27)
    warning(['MATLAB:cubLattice_g:wrongmmax',' The maximum exponent for the budget should be an integer smaller or equal to 27.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if (out_param.fudge <= 0 )
    warning(['MATLAB:cubLattice_g:fudgenonpos','The fudge factor should be greater than 0.' ...
            ' Using default fudge factor ' num2str(default.fudge)])
    out_param.fudge = default.fudge;
end


%%%%% Check diff

end

