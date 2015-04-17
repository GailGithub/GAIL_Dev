function [q,out_param] = cubSobol_g_CSC(varargin)
%cubSobol_g_CSC Quasi-Monte Carlo method using Sobol' cubature over the
%d-dimensional region to integrate within a specified generalized error
%tolerance with guarantees under Walsh-Fourier coefficients cone decay assumptions
%
%   [q,out_param] = cubSobol_g_CSC(f,d) estimates the integral of f over the
%   d-dimensional region with an error guaranteed not to be greater than 
%   a specific generalized error tolerance, 
%   tolfun:=max(abstol,reltol*|integral(f)|). The generalized tolerance function can
%   aslo be cosen as tolfun:=theta*abstol+(1-theta)*reltol*|integral(f)| 
%   where theta is another input parameter. Input f is a function handle. f should
%   accept an n x d matrix input, where d is the dimension of the hypercube,
%   and n is the number of points being evaluated simultaneously. The input d
%   is the dimension in which the function f is defined. Given the
%   construction of Sobol', d must be a positive integer with 1<=d<=1111.
%
%   q = cubSobol_g_CSC(f,d,abstol,reltol,density,mmin,mmax,fudge,errtype,theta)
%   estimates the integral of f over a d-dimensional region. The answer
%   is given within the generalized error tolerance tolfun. All parameters
%   should be input in the order specified above. If an input is not specified,
%   the default value is used. Note that if an input is not specified,
%   the remaining tail cannot be specified either.
%
%   q = cubSobol_g_CSC(f,d,'abstol',abstol,'reltol',reltol,'density',density,'mmin',mmin,'mmax',mmax,'fudge',fudge,'errtype',errtype,'theta',theta)
%   estimates the integral of f over a d-dimensional region. The answer
%   is given within the generalized error tolerance tolfun. All the field-value
%   pairs are optional and can be supplied with any order. If an input is not
%   specified, the default value is used.
%
%   q = cubSobol_g_CSC(f,d,in_param) estimates the integral of f over the
%   d-dimensional region. The answer is given within the generalized error 
%   tolerance tolfun.
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
%     in_param.reltol --- the relative error tolerance, which should be
%     in (0,1]. Default value is 1e-1.
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
%     By default it is @(x) 5*2.^-x.
%
%     in_param.errtype --- this is the tolerance function. There are two
%     choices, 'max' (chosen by default) which takes
%     max(abstol,reltol*|integral(f)|) and 'comb' which is a linear combination
%     theta*abstol+(1-theta)*reltol*|integral(f)|. Theta is another 
%     parameter that can be specified (see below).
% 
%     in_param.theta --- this input is parametrizing the errtype 
%     'comb'. Thus, it is only afecting when the errtype
%     chosen is 'comb'. It stablishes the linear combination weight
%     between the absolute and relative tolerances
%     theta*abstol+(1-theta)*reltol*|integral(f)|. Note that for theta=1, 
%     we have pure absolute tolerance while for theta=0, we have pure 
%     relative tolerance. By default, theta=1.
%
%   Output Arguments
%
%     q --- the estimated value of the integral.
%
%     out_param.n --- number of points used when calling cubSobol_g_CSC for f.
%
%     out_param.pred_err --- predicted bound on the error based on the cone
%     condition. If the function lies in the cone, the real error should be
%     smaller than this predicted error.
%
%     out_param.time --- time elapsed in seconds when calling cubSobol_g_CSC for f.
%
%     out_param.exit --- this is a number defining the conditions of
%     success or failure satisfied when finishing the algorithm. The 
%     algorithm is considered successful (with out_param.exit == 1) if no 
%     other flags arise warning that the results are certainly not 
%     guaranteed. The initial value is 1 and the final value of this
%     parameter is encoded as follows:
%     
%                       1    If reaching overbudget. It states whether
%                       the max budget is attained without reaching the
%                       guaranteed error tolerance.
%      
%                       2   If the function lies outside the cone. In
%                       this case, results are not guaranteed. Note that
%                       this parameter is computed on the transformed
%                       function, not the input function. For more
%                       information on the transforms, check the input
%                       parameter in_param.transfrom.
% 
%  Guarantee
% This algorithm computes the integral of real valued functions in [0,1)^d 
% with a prescribed generalized error tolerance. The Walsh-Fourier 
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
% >> f = @(x) x(:,1).*x(:,2); d = 2; q = cubSobol_g_CSC(f,d,1e-5,1e-1,'uniform')
% q = 0.25***
% 
% 
% Example 2:
% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2+0.11
% in the interval R^3 where x1, x2 and x3 are normally distributed:
%
% >> f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2+0.11; d = 3; q = cubSobol_g_CSC(f,d,1e-3,1e-3,'normal')
% q = 1.1***
% 
%
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [0,1)^2:
% 
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); d = 2; q = cubSobol_g_CSC(f,d,1e-3,1e-1,'uniform')
% q = 0.55***
%
%
% Example 4: 
% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.
% 
% >> f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); d = 1; q = cubSobol_g_CSC(f,d,1e-4,1e-1,'normal','fudge',@(m) 2.^-(2*m))
% q = 2.05***
%
%
% Example 5:
% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the interval
% [0,1)^5 with pure absolute error 1e-5.
% 
% >> f = @(x) 8*prod(x,2); d = 5; q = cubSobol_g_CSC(f,d,1e-5,'errtype','comb','theta',1)
% q = 0.25***
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
%   "GAIL: Guaranteed Automatic Integration Library (Version 2.1)"
%   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software", Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://code.google.com/p/gail/ 
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice And Experiences
%   (WSSSPE1)", Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

tic
%% Check and initialize parameters
[f,out_param] = cubSobol_g_CSC_param(varargin{:});

if strcmp(out_param.density,'normal')
   f=@(x) f(gail.stdnorminv(x));
end

%% Main algorithm
mlag=4; %distance between coefficients summed and those computed
sobstr=sobolset(out_param.d); %generate a Sobol' sequence
sobstr=scramble(sobstr,'MatousekAffineOwen'); %scramble it
Stilde=zeros(out_param.mmax-out_param.mmin+1,1); %initialize sum of DFWT terms
StildeNC=zeros(out_param.mmax-out_param.mmin+1,mlag); %initialize various sums of DFWT terms for necessary conditions
cond1=(1+out_param.fudge(mlag))*(1+2*out_param.fudge(mlag-(1:mlag)))./(1+out_param.fudge(mlag-(1:mlag))); % Factors for the necessary conditions
cond2=(1+out_param.fudge(mlag-(1:mlag)))*(1+2*out_param.fudge(mlag))/(1+out_param.fudge(mlag)); % Factors for the necessary conditions
errest=zeros(out_param.mmax-out_param.mmin+1,1); %initialize error estimates
appxinteg=zeros(out_param.mmax-out_param.mmin+1,1); %initialize approximations to integral
exit_len = 2;
out_param.exit=logical(zeros(exit_len,1)); %we start the algorithm with all warning flags down
outside_cone=false; %internal flag that becomes true if the function lies outside the cone


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
for i = 1:mlag % Storing the information for the necessary conditions
    nllstart=2*nllstart;
    StildeNC(i,i)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
end
out_param.pred_err=out_param.fudge(out_param.mmin)*Stilde(1);
errest(1)=out_param.pred_err;

deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.errtype)+gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.errtype));
deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
    out_param.reltol,out_param.theta,abs(q-errest(1)),...
    out_param.errtype)-gail.tolfun(out_param.abstol,out_param.reltol,...
    out_param.theta,abs(q+errest(1)),out_param.errtype));

is_done = false;
if out_param.pred_err <= deltaplus
   q=q+deltaminus;
   appxinteg(1)=q;
   out_param.time=toc;
   is_done = true;
elseif out_param.mmin == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
   out_param.exit(1)=1;
end

%% Loop over m
for m=out_param.mmin+1:out_param.mmax 
   if is_done,
       break;
   end    
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

   %% Update kappanumap
   kappanumap=[kappanumap; (nnext+1:out_param.n)']; %initialize map
   for l=m-1:-1:m-mlag-1
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
   for i = 1:mlag % Storing the information for the necessary conditions
       nllstart=2*nllstart;
       StildeNC(i+meff-1,i)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   end
   % disp((Stilde(meff)*cond1(1:min(meff-1,mlag)))>=(StildeNC(meff-1,1:min(meff-1,mlag)))) % Displaying necessary condition 1 results (1 if satisfied)
   % disp((StildeNC(meff-1,1:min(meff-1,mlag)).*cond2(1:min(meff-1,mlag)))>=(Stilde(meff)*ones(1,min(meff-1,mlag)))) % Displaying necessary condition 2 results (1 if satisfied)
   if ~(all((Stilde(meff)*cond1(1:min(meff-1,mlag)))>=(StildeNC(meff-1,1:min(meff-1,mlag))))*   ...
    all((StildeNC(meff-1,1:min(meff-1,mlag)).*cond2(1:min(meff-1,mlag)))>=(Stilde(meff)*ones(1,min(meff-1,mlag)))))  ...
    && outside_cone == false,
        outside_cone = true; % We are outside the cone
        out_param.exit(2) = 1;
   end
   out_param.pred_err=out_param.fudge(m)*Stilde(meff);
   errest(meff)=out_param.pred_err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   
    deltaplus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.errtype)+gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.errtype));
    deltaminus = 0.5*(gail.tolfun(out_param.abstol,...
        out_param.reltol,out_param.theta,abs(q-errest(meff)),...
        out_param.errtype)-gail.tolfun(out_param.abstol,out_param.reltol,...
        out_param.theta,abs(q+errest(meff)),out_param.errtype));
   
   if out_param.pred_err <= deltaplus
      q=q+deltaminus;
      appxinteg(meff)=q;
      out_param.time=toc;
      is_done = true;
   elseif m == out_param.mmax % We are on our max budget and did not meet the error condition => overbudget
      out_param.exit(1) = 1;
   end

end

exit_str='';
if sum(out_param.exit) == 0
  exit_str = '0';
else
  for i=1:exit_len
    if out_param.exit(i)==1,
      if i<exit_len 
        exit_str = strcat(exit_str,{num2str(i)}, {' '});
      else
        exit_str = strcat(exit_str,{num2str(i)});
      end
    end
  end
end
out_param.exitflag = exit_str;
out_param = rmfield(out_param,'exit');
out_param.time=toc;
end


%% Parsing for the input of cubSobol_g_CSC
function [f, out_param] = cubSobol_g_CSC_param(varargin)

% Default parameter values
default.abstol  = 1e-4;
default.reltol  = 1e-1;
default.density  = 'uniform';
default.mmin  = 10;
default.mmax  = 24;
default.fudge = @(m) 5*2.^-m;
default.errtype  = 'max';
default.theta  = 1;

if numel(varargin)<2
    help cubSobol_g_CSC
    warning('MATLAB:cubSobol_g_CSC:fdnotgiven',...
        'At least, function f and dimension d need to be specified. Example for f(x)=x^2:')
    f = @(x) x.^2;
    out_param.f=f;
    out_param.d=1;
else
    f = varargin{1};
    if ~gail.isfcn(f)
        warning('MATLAB:cubSobol_g_CSC:fnotfcn',...
            'The given input f was not a function. Example for f(x)=x^2:')
        f = @(x) x.^2;
        out_param.f=f;
        out_param.d=1;
    else
        out_param.f=f;
        d = varargin{2};
        if ~isnumeric(d) || ~gail.isposint(d) || ~(d<101)
            warning('MATLAB:cubSobol_g_CSC:dnotposint',...
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
        warning('MATLAB:cubSobol_g_CSC:validvarargin','Optional parameters must be numeric or strings. We will use the default optional parameters.')
    end
    in3=varargin{3};
end

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
  f_addParamVal = @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin   
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.density = default.density;
    out_param.mmin = default.mmin;
    out_param.mmax = default.mmax;  
    out_param.fudge = default.fudge;
    out_param.errtype = default.errtype;
    out_param.theta = default.theta;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'d',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        addOptional(p,'mmin',default.mmin,@isnumeric);
        addOptional(p,'mmax',default.mmax,@isnumeric);
        addOptional(p,'fudge',default.fudge,@gail.isfcn);
        addOptional(p,'errtype',default.errtype,...
            @(x) any(validatestring(x, {'max','comb'})));
        addOptional(p,'theta',default.theta,@isnumeric);
    else
        if isstruct(in3) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'density',default.density,...
            @(x) any(validatestring(x, {'uniform','normal'})));
        f_addParamVal(p,'mmin',default.mmin,@isnumeric);
        f_addParamVal(p,'mmax',default.mmax,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@gail.isfcn);
        f_addParamVal(p,'errtype',default.errtype,...
            @(x) any(validatestring(x, {'max','comb'})));
        f_addParamVal(p,'theta',default.theta,@isnumeric);
    end
    parse(p,f,d,varargin{3:end})
    out_param = p.Results;
end;

% Force absolute tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('MATLAB:cubSobol_g_CSC:abstolnonpos',['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Force relative tolerance greater than 0 and smaller than 1
if (out_param.reltol <= 0) || (out_param.reltol > 1)
    warning('MATLAB:cubSobol_g_CSC:reltolnonunit',['Relative tolerance should be chosen in (0,1].' ...
            ' Using default relative tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end

% Force density to be uniform or normal only
if ~(strcmp(out_param.density,'uniform') || strcmp(out_param.density,'normal') )
    warning('MATLAB:cubSobol_g_CSC:notdensity',['The density can only be uniform or normal.' ...
            ' Using default density ' num2str(default.density)])
    out_param.density = default.density;
end

% Force mmin to be integer greater than 0
if (~gail.isposint(out_param.mmin) || ~(out_param.mmin < out_param.mmax+1))
    warning('MATLAB:cubSobol_g_CSC:lowmmin',['The minimum starting exponent ' ...
            'should be an integer greater than 0 and smaller or equal than the maxium.' ...
            ' Using default mmin ' num2str(default.mmin)])
    out_param.mmin = default.mmin;
end

% Force exponent budget number of points be a positive integer greater than
% or equal to mmin an smaller than 54
if ~(gail.isposint(out_param.mmax) && out_param.mmax>=out_param.mmin && out_param.mmax<=53)
    warning('MATLAB:cubSobol_g_CSC:wrongmmax',['The maximum exponent for the budget should be an integer smaller or equal to 53.' ...
            ' Using default mmax ' num2str(default.mmax)])
    out_param.mmax = default.mmax;
end

% Force fudge factor to be greater than 0
if ~((gail.isfcn(out_param.fudge) && (out_param.fudge(1)>0)))
    warning('MATLAB:cubSobol_g_CSC:fudgenofcn',['The fudge factor should be a positve function.' ...
            ' Using default fudge factor ' func2str(default.fudge)])
    out_param.fudge = default.fudge;
end

% Force errtype to be max or comb
if ~(strcmp(out_param.errtype,'max') || strcmp(out_param.errtype,'comb') )
    warning('MATLAB:cubSobol_g_CSC:noterrtype',['The error type can only be max or comb.' ...
            ' Using default errtype ' num2str(default.errtype)])
    out_param.errtype = default.errtype;
end

% Force theta to be in [0,1]
if (out_param.theta < 0) || (out_param.theta > 1)
    warning('MATLAB:cubSobol_g_CSC:thetanonunit',['Theta should be chosen in [0,1].' ...
            ' Using default theta ' num2str(default.theta)])
    out_param.theta = default.theta;
end

end