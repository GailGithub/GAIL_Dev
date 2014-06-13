function [q,err,time,n]=cubLattice(f,d,abstol,density,shift,mmax, fudge,diff)
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
[f,out_param] = cubLattice_param(varargin{:});

tic
%% Initialize parameters
if nargin < 8
   diff='C1sin';
    if nargin < 7
        fudge=3;
        if nargin < 6
            mmax=24; % Maximum budget of points is 2^mmax. To check with lattice_gen.
            if nargin < 5
                shift=rand; % Random shift to our lattice rule.
                if nargin < 4
                   density='uniform';
                   if nargin < 3
                      abstol=1e-4;
                      if nargin < 2
                         d=1;
                      end
                   end
                end
            end
        end
    end
end
if strcmp(density,'normal')
   f=@(x) f(norminv(x));
end
if strcmp(diff,'Baker')
    f=@(x) f(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(diff,'C0')
    f=@(x) f(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(diff,'C1')
    f=@(x) f(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(diff,'C1sin')
    f=@(x) f(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
end

mmin=10; %initial number of points is 2^mmin
mlag=4;
Stilde=zeros(mmax-mmin+1,1);
errest=zeros(mmax-mmin+1,1);
appxinteg=zeros(mmax-mmin+1,1);

%% Initial points and FWT
n=2^mmin;
xpts=mod(lattice_gen(1,n,d)+shift,1); n0=n;
y=f(xpts);
yval=y;

%% Compute initial FFT
for l=0:mmin-1
   nl=2^l;
   nmminlm1=2^(mmin-l-1);
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
kappanumap=(1:n)'; %initialize map
for l=mmin-1:-1:1
   nl=2^l;
   oldone=abs(y(kappanumap(2:nl))); %earlier values of kappa, don't touch first one
   newone=abs(y(kappanumap(nl+2:2*nl))); %later values of kappa, 
   flip=find(newone>oldone); %
   temp=kappanumap(nl+1+flip);
   kappanumap(nl+1+flip)=kappanumap(1+flip);
   kappanumap(1+flip)=temp;
end

%% Compute Stilde
nllstart=int64(2^(mmin-mlag-1));
Stilde(1)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
err=fudge*2^(-mmin)*Stilde(1);
errest(1)=err;
if err <= abstol; time=toc; return, end

%% Loop over m
for m=mmin+1:mmax 
   n=2^m;
   mnext=m-1;
   nnext=2^mnext;
   xnext=mod(lattice_gen(nnext+1,2*nnext,d)+shift,1);
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
   kappanumap=[kappanumap; (nnext+1:n)']; %initialize map
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
   meff=m-mmin+1;
   Stilde(meff)=sum(abs(y(kappanumap(nllstart+1:2*nllstart))));
   err=fudge*2^(-m)*Stilde(meff);
   errest(meff)=err;

   %% Approximate integral
   q=mean(yval);
   appxinteg(meff)=q;
   if err <= abstol; time=toc; break
    elseif m==mmax; warning('The maximum budget is attained without reaching the tolerance.');
   end

end
time=toc;







function [f, out_param] = integral_g_param(varargin)
% Parsing for the input of cubLattice

% Default parameter values
default.d  = 1;
default.abstol  = 1e-4;
default.density  = 'uniform';
default.shift  = rand;
default.mmax  = 24;
default.fudge = 3;
default.diff = 'C1sin';

if isempty(varargin)
    help integral_g
    warning('Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
    out_param.f=f;
else
    f = varargin{1};
    out_param.f=f;
end;

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
    %if only one input f, use all the default parameters
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.absabstol = default.abstol;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;    
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'nlo',default.nlo,@isnumeric);
        addParamValue(p,'nhi',default.nhi,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('MATLAB:integral_g:anoinfinity',['a can not be infinity nor NaN. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('MATLAB:integral_g:bnoinfinity',['b can not be infinity not Nan. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;
if (out_param.b < out_param.a)
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
    flip=1;
end

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['MATLAB:integral_g:abstolnonpos','Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let initial number of points be a positive integer
if (~isposint(out_param.nlo))
    if isposge3(out_param.nlo)
        warning('MATLAB:integral_g:lowinitnotint',['Lowest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo))])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:integral_g:lowinitlt3',['Lowest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~isposint(out_param.nhi))
    if isposge3(out_param.nhi)
        warning('MATLAB:integral_g:highinitnotint',['Highest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:integral_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    if isposge3(out_param.nhi)
        warning('MATLAB:integral_g:nlobtnhi',['Highest initial number of points should be at least equal to to lowest initial number of points.' ...
            ' Using ', num2str(ceil(out_param.nhi)), ' as nlo'])
        out_param.nlo = ceil(out_param.nhi);
    else
        warning('MATLAB:integral_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

out_param.ninit = max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+(out_param.b-out_param.a)))),3);
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('MATLAB:integral_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:integral_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end
