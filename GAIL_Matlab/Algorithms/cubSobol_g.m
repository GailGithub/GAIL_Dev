function [q,out_param]=cubSobol_g(varargin)
% cubSobol_g is a Quasi-Monte Carlo method using Sobol' cubature over the
% d-multidimensional region to integrate within a specified absolute error 
% tolerance with guarantees under Fourier coefficients cone decay assumptions.
%
% [q,out_param] = cubSobol_g(f,d) estimates the integral of f over the
% d-dimensional region to within a predefined absolute error tolerance
% 1e-4 and with guaranteed error. Input f is a function handle. f should
% accept an n x d matrix input, where d is the dimension of the hypercube,
% and n is the number of points being evaluated simultaneously. The input d
% is the dimension in which the function f is defined. Given the
% construction of Sobol', d must be a positive integer with 1<=d<=1111.
%
% q = cubSobol_g(f,d,abstol,density,mmin,mmax,fudge)
% estimates the integral of f over the d-dimensional region. The answer
% is given within the absolute error tolerance abstol. All parameters
% should be input in the order specified above. If an input is not specified,
% the default value is used. Note that if an input is not specified,
% the remaining tail can not be specified either.
%
% q = cubSobol_g(f,d,'abstol',abstol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge)
% estimates the integral of f over the d-dimensional region. The answer
% is given within the absolute error tolerance abstol. All the field-value
% pairs are optional and can be supplied with any order. If an input is not
% specified, the default value is used.
%
% q = cubSobol_g(f,d,in_param) estimates the integral of f over the
% d-dimensional region. The answer is given within the absolute error 
% tolerance in_param.abstol.
%
% f --- the integrand.
%
% d --- dimension where f is defined. d must be a positive integer 1<=d<=1111.
%
% in_param.abstol --- the absolute error tolerance, abstol>0. By default is 1e-4. 
%
% in_param.density --- for f(x), we can define x uniform in [0,1)^d or
% normally distributed with covariance matrix Id^d. By default is
% 'uniform'. The only possible values are 'uniform' or 'normal'.
%
% in_param.mmin --- the minimum number of points to start is 2^mmin. The
% cone condition on the Fourier coefficients decay requires a minimum
% number of points to start. The advice is to consider at least mmin=10.
% mmin needs to be a positive integer with mmin<=mmax. By default is 10.
%
% in_param.mmax --- the maximum budget is 2^mmax. By construction of the
% Sobol' generator, mmax is a positive integer such that mmin<=mmax<=53.
% The default value is 24.
%
% in_param.fudge --- the constant multiplying the cone of functions. For more
% information about this parameter, refer to the references. It should be a
% real positve number. By default is 3.
%
% q --- the estimated value of the integral.
%
% out_param.overbudget --- string stating whether the max budget is
% attained without reaching the guaranteed error tolerance.
%
% out_param.n --- number of points used when calling cubSobol_g for f.
%
% out_param.pred_err --- predicted bound on the error based on the cone
% condition. If the function lies in the cone, the real error should be
% smaller than this predicted error.
%
% out_param.time --- time elapsed when calling cubSobol_g for f.
%
%
% GUARANTEE
% ---------
% The guarantee is based on the decay of the Fourier coefficients of the
% given inegrand. For more details on how the cone is defined, please refer
% to the references below.
%
%
% EXAMPLES
% --------
%
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the
% interval [0,1)^2:
%
% >> f=@(x) x(:,1).*x(:,2); d=2;
% >> q=cubSobol_g(f,d,1e-5,'uniform')
% q = 0.25***
% 
% 
% Example 2:
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:
%
% >> f=@(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d=3;
% >> q=cubSobol_g(f,d,1e-3,'normal')
% q = 1.1***
% 
%
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2); d=2;
% >> q=cubSobol_g(f,d,1e-3,'uniform')
% q = 0.55***
%
%
% Example 4: 
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
% 
% >> f=@(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d=1;
% >> q=cubSobol_g(f,d,1e-4,'normal')
% q = 2.05***
%
%
%
%
% See also cubLattice_g, cubMC_g, meanMC_g, integral_g
% 
% References
%
% [1]  F. J. Hickernell, Lluis Antoni Jimenez Rugama
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluís Antoni Jiménez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% "GAIL: Guaranteed Automatic Integration Library (Version 2.0.0)"
% [MATLAB Software], 2014. Available from http://code.google.com/p/gail/
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.


tic
%% Check and initialize parameters
[f,out_param] = cubSobol_g_param(varargin{:});

if strcmp(out_param.density,'normal')
   f=@(x) f(norminv(x));
end

%% Main algorithm
mlag=4; %distance between coefficients summed and those computed
sobstr=sobolset(out_param.d); %generate a Sobol' sequence
sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde=zeros(out_param.mmax-out_param.mmin+1,1);
errest=zeros(out_param.mmax-out_param.mmin+1,1);
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1);
out_param.overbudget='Max budget reached with no guarantees.';

%% Initial points and FWT
out_param.n=2^out_param.mmin;
xpts=sobstr(1:out_param.n,1:out_param.d); n0=out_param.n;
y=f(xpts);
yval=y;

%% Compute initial FWT
for l=0:out_param.mmin-1
   nl=2^l;
   nmminlm1=2^(out_param.mmin-l-1);
   ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
end
%y now contains the FWT coefficients

%% Approximate integral
q=mean(yval);
appxinteg(1)=q;

%% Create kappanumap
kappanumap=(1:out_param.n)'; %initialize map
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone);
   temp=kappanumap(nl+1+flip);
   kappanumap(nl+1+flip)=kappanumap(1+flip);
   kappanumap(1+flip)=temp;
end

%% Compute Stilde
nllstart=2^(out_param.mmin-mlag-1);
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
out_param.pred_err=out_param.fudge*2^(-out_param.mmin)*Stilde(1);
errest(1)=out_param.pred_err;
if out_param.pred_err <= out_param.abstol; out_param.overbudget='Max budget not reached.'; out_param.time=toc; return, end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax 
   out_param.n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=sobstr(n0+(1:nnext),1:out_param.d); 
   n0=n0+nnext;
   ynext=f(xnext);
   yval=[yval; ynext];

   %% Compute initial FWT on next points
   for l=0:mnext-1
      nl=2^l;
      nmminlm1=2^(mnext-l-1);
      ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
      evenval=ynext(ptind);
      oddval=ynext(~ptind);
      ynext(ptind)=(evenval+oddval)/2;
      ynext(~ptind)=(evenval-oddval)/2;
   end

   %% Compute FWT on all points
   y=[y;ynext];
   nl=2^mnext;
   ptind=[true(nl,1); false(nl,1)];
   evenval=y(ptind);
   oddval=y(~ptind);
   y(ptind)=(evenval+oddval)/2;
   y(~ptind)=(evenval-oddval)/2;
   %disp('line 100')

   %% Update kappanumap
   kappanumap=[kappanumap; (nnext+1:out_param.n)']; %initialize map
   for l=m-1:-1:1
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone); %
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end

   %% Compute Stilde
   nllstart=2^(m-mlag-1);
   meff=m-out_param.mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   out_param.pred_err=out_param.fudge*2^(-m)*Stilde(meff);
   errest(meff)=out_param.pred_err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if out_param.pred_err <= out_param.abstol; out_param.overbudget='Max budget not reached.'; out_param.time=toc; return, end

end
out_param.time=toc;
end


%% Parsing for the input of cubSobol_g
function [f, out_param] = cubSobol_g_param(varargin)

% Default parameter values
default.abstol  = 1e-4;
default.density  = 'uniform';
default.mmin  = 10;
default.mmax  = 24;
default.fudge = 3;

if numel(varargin)<2
    help cubSobol_g
    warning('MATLAB:cubSobol_g:fdnotgiven',...
        'At least, function f and dimension d of f must be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    out_param.d=1;
else
    f = varargin{1};
    out_param.f=f;
    d = varargin{2};
    if ~isnumeric(d) || ~isposint(d)
        warning('MATLAB:cubSobol_g:dnotposint',...
            'The dimension d of f must be a positive integer. Example for f(x)=x^2:')
        f = @(x) x.^2;
        out_param.f=f;
        out_param.d=1;
    else
    out_param.d=d;
    end
end;

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin{3};
    validvarargin=(isnumeric(in3) || isstruct(in3) ...
        || ischar(in3));
end

if ~validvarargin
    % If only one input f or in2 is nothing above, use all the default parameters
    %warning(['MATLAB:cubSobol_g:mininputarg',' Optional parameters must be numeric or strings.'])
    out_param.abstol = default.abstol;
    out_param.density = default.density;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    addRequired(p,'d',@isnumeric);
    if isnumeric(in3) || ischar(in3) %if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addParamValue(p,'mmin',default.mmin,@isnumeric);
        addParamValue(p,'mmax',default.mmax,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
    end
    parse(p,f,d,varargin{3:end})
    out_param = p.Results;
end;

% Force error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['MATLAB:cubSobol_g:abstolnonpos','Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

%%%% Check density (to code)

% Force mmin to be integer greater than 0
if (out_param.mmin < 1 || ~isposint(out_param.mmin))
    warning(['MATLAB:cubSobol_g:lowmmin',' The minimum starting exponent should be an integer greater or equal than 1.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater or equal than
% mmin an smaller than 54
if ~(isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning(['MATLAB:cubSobol_g:wrongmmax',' The maximum exponent for the budget should be an integer smaller or equal to 53.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if (out_param.fudge <= 0 )
    warning(['MATLAB:cubSobol_g:fudgenonpos','The fudge factor should be greater than 0.' ...
            ' Using default fudge factor ' num2str(default.fudge)])
    out_param.fudge = default.fudge;
end

end