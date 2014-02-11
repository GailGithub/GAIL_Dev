function [q,out_param] = integralab_g_CSC(varargin)
%  integralab_g_CSC 1-D guaranteed function integration using trapezoidal rule
% 
%  Description
%
%   q = integralab_g_CSC(f) computes q, the definite integral of function f
%   on the interval [0,1] by trapezoidal rule with 
%   in a guaranteed absolute error of 1e-6. Default starting number of
%   sample points taken is 100 and default cost budget is 1e7. Input f is a 
%   function handle. The function y = f(x) should accept a vector argument 
%   x and return a vector result y, the integrand evaluated at each element
%   of x.
%
%   q = integralab_g_CSC(f,in_param) computes q, the definite integral of 
%   function f by trapezoidal rule within a guaranteed absolute error
%   in_param.abstol, starting number of points in_param.ninit, and cost 
%   budget in_param.nmax. If a field is not specified, the default value is
%   used.
%   
%   in_param.abstol --- absolute error tolerance required by user.
% 
%   in_param.nlo --- lowest initial number of function values used
%
%   in_param.nhi --- highest initial number of function values used
%
%   in_param.nmax --- cost budget (maximum number of function values)
% 
%   in_param.a --- low end of the integral
%
%   in_param.b --- high end of the integral
%
%   q = INTEGRAL_G(f,'abstol',abstol,'nlo',nlo,'nhi',nhi,,'a',a,'b',b,'nmax',nmax) computes
%   q, the definite integral of function f by trapezoidal rule within a 
%   guaranteed absolute error tolerance abstol, starting number of points 
%   ninit, and cost budget nmax. All three field-value pairs are optional 
%   and can be supplied.
%
%   [q, out_param] = INTEGRAL_G(f,...) returns the approximated 
%   integration q and output structure out_param, which includes the 
%   fileds in_param plus the following fields:
%
%   out_param.exceedbudget --- it is true if the algorithm tries to use 
%   more points than cost budget, false otherwise.
% 
%   out_param.tauchange --- it is true if the cone constant has been
%   changed, false otherwise. See [1] for details. If true, you may wish to
%   change the input in_param.ninit to a larger number.
% 
%   out_param.npoints --- number of points we need to 
%   reach the guaranteed absolute error tolerance abstol.
%
%   out_param.errest --- approximation error defined as the differences
%   between the true value and the approximated value of the integral.
%
%   out_param.nlo --- lowest initial number of function values
%
%   out_param.nhi --- highest initial number of function values
%
%   out_param.nmax --- cost budget (maximum number of function values)
% 
%   out_param.a --- low end of the integral
%
%   out_param.b --- high end of the integral
% 
%
%   Examples
%
%   Example 1: 
%   >> q = integralab_g_CSC(@(x) x.^2)
%   q = 0.3333
%
%
%   Example 2:
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.1353 
%
%
%   Example 3:
%   >> q = integralab_g_CSC()
%   Warning: Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.
%   >  In ***
%   q = 0.3333
%
%
%
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'abstol',1e-5,'nhi',52,'nmax',1e7)
%   q = 0.7468
%
%
%
%   >> f = @(x) exp(-(x-1).^2); q = integralab_g_CSC(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.7468
%
%  
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',0,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.8821
%
%
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',0,'b',3,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.8862 
%
%
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',-1,'b',3,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 1.6330 
%
%
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',-4.5,'b',1.5,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 1.7424
%
%
%  >>  [~,out]=integralab_g_CSC(@(x) x.^2) 
% 
% out = 
%     abstol: 1.0000e-06
%                    f: @(x)x.^2
%                ninit: 100
%                 nmax: 10000000
%                  tau: 197
%         exceedbudget: 0
%            tauchange: 0
%                    q: 0.3333
%              npoints: 3565
%               errest: ***e-07
%
%
%   >> f = @(x) exp(-x.^2); [~, out_param] = integralab_g_CSC(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%         out_param = 
%        abstol: 1.0000e-05
%                      f: @(x)exp(-x.^2)
%                  ninit: 1000
%                   nmax: 10000000
%                    tau: 1997
%           exceedbudget: 0
%              tauchange: 0
%                      q: 0.1353
%                npoints: 2998
%                 errest: 7.3718e-06    
%
%
%   >> f = @(x) exp(-x.^2); q = integralab_g_CSC(f,1,2,1e-5,100,10000)
%   q = 0.1353 
%
%
%  >>  inparam.a=0; inparam.b=3; inparam.abstol=1e-13; q=integralab_g_CSC(@(x) exp(2*x),inparam)
%  q =  201.2144
% 
%
% See also funappxab_g, cubMC_g
%
% Reference:
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%      The complexity of guaranteed automatic algorithms: Cones, not
%      balls, Journal of Complexity 2013, to appear, DOI
%      10.1016/j.jco.2013.09.002.

%%
% check parameter satisfy conditions or not
[f,out_param] = integralab_g_CSC_param(varargin{:});
intervallen=out_param.b-out_param.a;
if (intervallen==0)
  q=out_param.q;
  return;
end

%% main alg
out_param.tau=ceil((out_param.ninit-1)*2-1); % computes the minimum requirement of number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of q is less than cost budget
out_param.tauchange=false;  % if the cone constant has been changed
newf=transformIntegrand(f,out_param);
[q,out_param]=integral_g(newf,out_param);


function [f, out_param] = integralab_g_CSC_param(varargin)
% parse the input to the integral_g function

% Default parameter values
default.a = 0;
default.b = 1;
default.abstol  = 1e-6;
default.nlo = 10;
default.nhi = 1000;
default.nmax  = 1e7;


if isempty(varargin)    
    warning('MATLAB:integralab_g_CSC:nofunction','Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    help integralab_g_CSC
    f = @(x) x.^2;
    out_param.f = f;
else
    f = varargin{1};
     out_param.f = f;
end;

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) || ischar(in2));
end

if ~validvarargin
    %if only one input f, use all the default parameters
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
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
    warning('MATLAB:integralab_g_CSC:anoinfinity',['a cannot be infinity or NaN. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('MATLAB:integralab_g_CSC:anoinfinity',['b cannot be infinity or NaN. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('MATLAB:integralab_g_CSC:blea','b cannot be smaller than a; exchange these two. The result should be negative of q.')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('MATLAB:integralab_g_CSC:beqa','b equals a.')
    out_param.q=0;
    out_param.npoints=0;
    out_param.errest=0;
    return;
end;

% let cost budget be a positive integer
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('MATLAB:integralab_g_CSC:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:integralab_g_CSC:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% let cost budget be a positive integer
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('MATLAB:integralab_g_CSC:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:integralab_g_CSC:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

% let initial number of points be a positive integer
if (out_param.nlo > out_param.nhi)
    warning('MATLAB:integralab_g_CSC:logrhi', 'Lower bound of initial number of points is larger than upper bound of initial number of points; exchange these two')
    temp = out_param.nlo;
    out_param.nlo = out_param.nhi;
    out_param.nhi = temp;
end;

if (~isposint(out_param.nlo))
    if isposge3(out_param.nlo)
        warning('MATLAB:integralab_g_CSC:lowinitnotint',['Lower initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo))])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:integralab_g_CSC:lowinitlt3',['Lower initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end

if (~isposint(out_param.nhi))
    if isposge3(out_param.nhi)
        warning('MATLAB:integralab_g_CSC:highinitnotint',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:integralab_g_CSC:highinitlt3',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

if (out_param.nlo > out_param.nhi)
    if isposge3(out_param.nhi)
        warning('MATLAB:integralab_g_CSC:nlobtnhi',['Highest initial number of points should be at least equal to to lowest initial number of points.' ...
            ' Using ', num2str(ceil(out_param.nhi)), ' as nlo'])
        out_param.nlo = ceil(out_param.nhi);
    else
        warning('MATLAB:integralab_g_CSC:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+(out_param.b-out_param.a))));
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('MATLAB:integralab_g_CSC:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:integralab_g_CSC:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

if 0,
if (~isfinite(f(out_param.a)))
  warning('MATLAB:integralab_g_CSC:unboundedfa', ['f(a) seems to be unbounded or undefined. Use a = ' num2str(out_param.a+eps)]);
  out_param.a = out_param.a+eps;
end

if (~isfinite(f(out_param.b)))
  warning('MATLAB:integralab_g_CSC:unboundedfb', ['f(b) seems to be unbounded or undefined. Use b = ' num2str(out_param.b-eps)]);
  out_param.b = out_param.b-eps;
end
end

function newf=transformIntegrand(oldf,out_param)
  % Transform integrand linearly  
  b=out_param.b;
  a=out_param.a;
  len = b-a;
  newf = @(x) (len*oldf(x.*len+a));
 


%% doctest results%%
% doctest integralab_g_CSC
% TAP version 13
% 1..3
% ok 1 - "q = integralab_g_CSC(@(x) x.^2)"
% ok 2 - "f = @(x) exp(-x.^2); q = integralab_g_CSC(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)"
% ok 3 - "q = integralab_g_CSC()"
