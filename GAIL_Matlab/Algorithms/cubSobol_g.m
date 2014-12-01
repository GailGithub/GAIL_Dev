function [q,out_param] = cubSobol_g(varargin)
%CUBSOBOL_G is a Quasi-Monte Carlo method using Sobol' cubature over the
%d-dimensional region to integrate within a specified absolute error 
%tolerance with guarantees under Walsh-Fourier coefficients cone decay assumptions.
%
%   [q,out_param] = CUBSOBOL_G(f,d) estimates the integral of f over the
%   d-dimensional region with an error guaranteed not to be greater than the
%   predefined error tolerance 1e-4. Input f is a function handle. f should
%   accept an n x d matrix input, where d is the dimension of the hypercube,
%   and n is the number of points being evaluated simultaneously. The input d
%   is the dimension in which the function f is defined. Given the
%   construction of Sobol', d must be a positive integer with 1<=d<=1111.
%
%   q = CUBSOBOL_G(f,d,abstol,density,mmin,mmax,fudge)
%   estimates the integral of f over a d-dimensional region. The answer
%   is given within the absolute error tolerance abstol. All parameters
%   should be input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail can not be specified either.
%
%   q = CUBSOBOL_G(f,d,'abstol',abstol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge)
%   estimates the integral of f over a d-dimensional region. The answer
%   is given within the absolute error tolerance abstol. All the field-value
%   pairs are optional and can be supplied with any order. If an input is not
%   specified, the default value is used.
%
%   q = CUBSOBOL_G(f,d,in_param) estimates the integral of f over the
%   d-dimensional region. The answer is given within the absolute error 
%   tolerance in_param.abstol.
% 
%   Input Arguments
%
%     f --- the integrand whose input should be a matrix nxd where n is the
%     number of data points and d the dimension. By default it is the
%     quadratic function.
%
%     d --- dimension of domain on which f is defined. d must be a positive
%     integer 1<=d<=1111. By default it is 1.
%
%     in_param.abstol --- the absolute error tolerance, abstol>0. By 
%     default it is 1e-4.
%
%     in_param.density --- for f(x)*mu(dx), we can define mu(dx) to be the
%     density function of a uniformly distributed random variable in [0,1)^d
%     or normally distributed with covariance matrix I_d. By default it 
%     is 'uniform'. The only possible values are 'uniform' or 'normal'.
%
%     in_param.mmin --- the minimum number of points to start is 2^mmin. The
%     cone condition on the Fourier coefficients decay requires a minimum
%     number of points to start. The advice is to consider at least mmin=10.
%     mmin needs to be a positive integer with mmin<=mmax. By default it is 10.
%
%     in_param.mmax --- the maximum budget is 2^mmax. By construction of the
%     Sobol' generator, mmax is a positive integer such that mmin<=mmax<=53.
%     The default value is 24.
%
%     in_param.fudge --- the positive function multiplying the finite 
%     sum of Fast Walsh coefficients specified in the cone of functions.
%     For more information about this parameter, refer to the references.
%     By default it is @(x) 5*2^-x.
%
%   Output Arguments
%
%     q --- the estimated value of the integral.
%
%     out_param.overbudget --- boolean stating whether the max budget is
%     attained without reaching the guaranteed error tolerance. Output 1
%     means we have overrun our budget.
%
%     out_param.n --- number of points used when calling cubSobol_g for f.
%
%     out_param.pred_err --- predicted bound on the error based on the cone
%     condition. If the function lies in the cone, the real error should be
%     smaller than this predicted error.
%
%     out_param.time --- time elapsed in seconds when calling cubSobol_g for f.
% 
%  Guarantee
% This algorithm computes the integral of real valued functions in [0,1)^d 
% with a prescribed absolute error tolerance. The Walsh-Fourier 
% coefficients of the integrand are assumed to be absolutely convergent.
% If the algorithm terminates without warning messages, the output is 
% given with guarantees under the assumption that the integrand lies inside
% a cone of functions. The guarantee is based on the decay rate of the 
% Walsh-Fourier coefficients. For more details on how the cone is defined, 
% please refer to the references below.
% 
%  Examples
% 
% Example 1:
% Estimate the integral with integrand f(x) = x1.*x2 in the interval [0,1)^2:
%
% >> f=@(x) x(:,1).*x(:,2); d=2; q = cubSobol_g(f,d,1e-5,'uniform')
% q = 0.25***
% 
% 
% Example 2:
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:
%
% >> f=@(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d=3; q = cubSobol_g(f,d,1e-3,'normal')
% q = 1.1***
% 
%
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2); d=2; q = cubSobol_g(f,d,1e-3,'uniform')
% q = 0.55***
%
%
% Example 4: 
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
% 
% >> f=@(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d=1; q = cubSobol_g(f,d,1e-4,'normal','fudge',@(x) 2^-(2*x))
% q = 2.05***
%
%
%   See also CUBLATTICE_G, CUBMC_G, MEANMC_G, INTEGRAL_G
% 
%  References
%
%   [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama: Reliable adaptive
%   cubature using digital sequences (2014). Submitted for publication:
%   arXiv:1410.8615.
%
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   "GAIL: Guaranteed Automatic Integration Library (Version 2.0.0)"
%   [MATLAB Software], 2014. Available from http://code.google.com/p/gail/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.

tic
%% Check and initialize parameters
[f,out_param] = cubSobol_g_param(varargin{:});

if strcmp(out_param.density,'normal')
   f=@(x) f(gail.stdnorminv(x));
end

%% Main algorithm
mlag=4; %distance between coefficients summed and those computed
sobstr=sobolset(out_param.d); %generate a Sobol' sequence
sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde=zeros(out_param.mmax-out_param.mmin+1,1); %initialize sum of DFWT terms
errest=zeros(out_param.mmax-out_param.mmin+1,1); %initialize error estimates
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
out_param.overbudget=true; %we have overrun our budget, until indicated otherwise

%% Initial points and FWT
out_param.n=2^out_param.mmin; %total number of points to start with
n0=out_param.n; %initial number of points
xpts=sobstr(1:n0,1:out_param.d); %grab Sobol' points
y=f(xpts); %evaluate integrand
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

%% Create kappanumap implicitly from the data
kappanumap=(1:out_param.n)'; %initialize map
for l=out_param.mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %which in the pair are the larger ones
   temp=kappanumap(nl+1+flip); %then flip 
   kappanumap(nl+1+flip)=kappanumap(1+flip); %them
   kappanumap(1+flip)=temp; %around
end

%% Compute Stilde
nllstart=2^(out_param.mmin-mlag-1);
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
out_param.pred_err=out_param.fudge(out_param.mmin)*Stilde(1);
errest(1)=out_param.pred_err;
if out_param.pred_err <= out_param.abstol
   out_param.overbudget=false; 
   out_param.time=toc; 
   return
end

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
%   for l=m-1:-1:1
   for l=m-1:-1:m-mlag-1 %update just some, not exactly sure about this
      nl=2^l;
      oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
      newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
      flip=find(newone>oldone);
      temp=kappanumap(nl+1+flip);
      kappanumap(nl+1+flip)=kappanumap(1+flip);
      kappanumap(1+flip)=temp;
   end

   %% Compute Stilde
   nllstart=2^(m-mlag-1);
   meff=m-out_param.mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   out_param.pred_err=out_param.fudge(m)*Stilde(meff);
   errest(meff)=out_param.pred_err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if out_param.pred_err <= out_param.abstol
      out_param.overbudget=false; 
      out_param.time=toc; 
      return 
   end

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
default.fudge = @(x) 5*2^-x;

if numel(varargin)<2
    help cubSobol_g
    warning('MATLAB:cubSobol_g:fdnotgiven',...
        'At least, function f and dimension d need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    out_param.d=1;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('MATLAB:cubSobol_g:fnotfcn',...
            'The given input f was not a function. Example for f(x)=x^2:')
        f = @(x) x.^2;
        out_param.f=f;
        out_param.d=1;
    else
        out_param.f=f;
        d = varargin{2};
        if ~isnumeric(d) || ~gail.isposint(d) || ~(d<101)
            warning('MATLAB:cubSobol_g:dnotposint',...
                'The dimension d must be a positive integer less than 101. Example for f(x)=x^2:')
            f = @(x) x.^2;
            out_param.f=f;
            out_param.d=1;
        else
        out_param.d=d;
        end
    end
end;

validvarargin=numel(varargin)>2;
if validvarargin
    in3=varargin(3:end);
    for j=1:numel(varargin)-2
    validvarargin=validvarargin && (isnumeric(in3{j}) ...
        || ischar(in3{j}) || isstruct(in3{j}) || gail.isfcn(in3{j}));
    end
    if ~validvarargin
        warning('MATLAB:cubSobol_g:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{3};
end

if ~validvarargin   
    out_param.abstol = default.abstol;
    out_param.density = default.density;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'d',@isnumeric);
    if isnumeric(in3) || ischar(in3) %if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
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
        addParamValue(p,'fudge',default.fudge,@gail.isfcn);
    end
    parse(p,f,d,varargin{3:end})
    out_param = p.Results;
end;

% Force error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('MATLAB:cubSobol_g:abstolnonpos',['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force density to be uniform or normal only
if ~(strcmp(out_param.density,'uniform') || strcmp(out_param.density,'normal') )
    warning('MATLAB:cubSobol_g:notdensity',['The density can only be uniform or normal.' ...
            ' Using default density ' num2str(default.density)])
    out_param.density = default.density;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('MATLAB:cubSobol_g:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('MATLAB:cubSobol_g:wrongmmax',['The maximum exponent for the budget should be an integer smaller or equal to 53.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('MATLAB:cubSobol_g:fudgenofcn',['The fudge factor should be a positve function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end
end
