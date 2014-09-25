function [pp,out_param]=funappx_g(varargin)
%FUNAPPX_G 1-D guaranteed function recovery on a closed interval [a,b]
%
%   pp = FUNAPPX_G(f) approximates function f on the default interval [0,1]
%   by a piecewise polynomial structure pp within the guaranteed absolute
%   error tolerance of 1e-6. Default initial number of points is 100 and
%   default cost budget is 1e7.  Input f is a function handle. The
%   statement y = f(x) should accept a vector argument x and return a
%   vector y of function values that is of the same size as x. Output pp
%   may be evaluated via PPVAL.
%
%   pp = FUNAPPX_G(f,a,b,abstol,nlo,nhi,nmax) for a given function f and
%   the ordered input parameters that define the finite interval [a,b], a
%   guaranteed absolute error tolerance abstol, a lower bound of initial
%   number of points nlo, an upper bound of initial number of points nhi,
%   and a cost budget nmax.
%
%   pp = FUNAPPX_G(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax)
%   recovers function f on the finite interval [a,b], given a guaranteed
%   absolute error tolerance abstol, a lower bound of initial number of
%   points nlo, an upper bound of initial number of points nhi, and a cost
%   budget nmax. All six field-value pairs are optional and can be supplied
%   in different order.
%
%   pp = FUNAPPX_G(f,in_param) recovers function f on the finite interval
%   [in_param.a,in_param.b], given a guaranteed absolute error tolerance
%   in_param.abstol, a lower bound of initial number of points
%   in_param.nlo, an upper bound of initial number of points in_param.nhi,
%   and a cost budget in_param.nmax. If a field is not specified, the
%   default value is used.
%
%   [pp, out_param] = FUNAPPX_G(f,...) returns a piecewise polynomial
%   structure pp and an output structure out_param.
%
%   Input Arguments
%
%     in_param.a --- left end point of interval, default value is 0
%
%     in_param.b --- right end point of interval, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6
%
%     in_param.nlo --- lower bound of initial number of points we used,
%     default value is 10
%
%     in_param.nhi --- upper bound of initial number of points we used,
%     default value is 1000
%
%     in_param.nmax --- cost budget, default value is 1e7
%
%   Output Arguments
%
%     pp.form --- pp means piecewise polynomials
%
%     pp.breaks --- show the location of interpolation points
%
%     pp.coefs --- coefficients for piecewise linear polynomials
%
%     pp.pieces --- number of piecewise linear polynomials
%
%     pp.order --- be 2 as we use piecewise linear polynomials
%
%     pp.dim --- be 1 as we do univariate approximation
%
%     pp.orient --- always be 'first'
%
%     out_param.exceedbudget --- it is 0 if the number of points used in
%     the construction of pp is less than cost budget, 1 otherwise.
%
%     out_param.ninit --- initial number of points we use
%
%     out_param.npoints --- number of points we need to reach the
%     guaranteed absolute error tolerance
%
%     out_param.errorbound --- an upper bound of the absolute error
%
%     out_param.nstar --- final value of the parameter defining the cone of
%     functions for which this algorithm is guaranteed; nstar = ninit-2
%     initially and is increased as necessary
%
%     out_param.a --- left end point of interval
%
%     out_param.b --- right end point of interval
%
%     out_param.abstol --- guaranteed absolute error tolerance
%
%     out_param.nlo --- a lower bound of initial number of points we use
%
%     out_param.nhi --- an upper bound of initial number of points we use
%
%     out_param.nmax --- cost budget
%
%  Guarantee
%
%  If the function to be approximated, f, satisfies the cone condition
%                          2 nstar   ||     f(b)-f(a)  ||
%      ||f''||        <=  ---------  ||f'- ----------- ||
%             \infty        b - a    ||       b - a    ||\infty,
%  then the pp output by this algorithm is guaranteed to satisfy
%      ||f-ppval(pp, )||\infty <= abstol,
%  and the upper bound of the cost is
%          ____________________________
%         / nstar*(b-a)^2 ||f''||\infty
%        / ---------------------------- + 2 nstar + 4
%      \/          2 abstol
%
%  provided the flag exceedbudget = 0.
%
%
%  Examples
%
%   Example 1:
%
%
%   >> f = @(x) x.^2; [pp, out_param] = funappx_g(f)
%
%  pp =
%
%       form: 'pp'
%     breaks: [1x9901 double]
%      coefs: [9900x2 double]
%     pieces: 9900
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param =
%
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%             nmax: 10000000
%            ninit: 100
%            nstar: 98
%     exceedbudget: 0
%          npoints: 9901
%       errorbound: 2.5508e-09
%
%
%   Example 2:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappx_g(f,-2,2,1e-7,10,10,1000000)
%
% pp =
%
%       form: 'pp'
%     breaks: [1x33733 double]
%      coefs: [33732x2 double]
%     pieces: 33732
%      order: 2
%        dim: 1
%     orient: 'first'
%
% out_param =
%
%                a: -2
%           abstol: 1.0000e-07
%                b: 2
%                f: @(x)x.^2
%              nhi: 10
%              nlo: 10
%             nmax: 1000000
%            ninit: 10
%            nstar: 8
%     exceedbudget: 0
%          npoints: 33733
%       errorbound: 3.5154e-09
%
%
%   Example 3:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappx_g(f,'a',-2,'b',2,'nhi',100,'nlo',10)
%
% pp =
%
%       form: 'pp'
%     breaks: [1x31249 double]
%      coefs: [31248x2 double]
%     pieces: 31248
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param =
%
%                a: -2
%           abstol: 1.0000e-06
%                b: 2
%                f: @(x)x.^2
%              nhi: 100
%              nlo: 10
%             nmax: 10000000
%            ninit: 64
%            nstar: 62
%     exceedbudget: 0
%          npoints: 31249
%       errorbound: 4.0965e-09
%
%
%   Example 4:
%
%   >> in_param.a = -10; in_param.b = 10;
%   >> in_param.abstol = 10^(-7); in_param.nlo = 10; in_param.nhi = 100;
%   >> in_param.nmax = 10^6; f = @(x) x.^2;
%   >> [pp, out_param] = funappx_g(f,in_param)
%
% pp =
%
%       form: 'pp'
%     breaks: [1x590071 double]
%      coefs: [590070x2 double]
%     pieces: 590070
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param =
%
%                a: -10
%           abstol: 1.0000e-07
%                b: 10
%                f: @(x)x.^2
%              nhi: 100
%              nlo: 10
%             nmax: 1000000
%            ninit: 90
%            nstar: 88
%     exceedbudget: 0
%          npoints: 590071
%       errorbound: 2.8721e-10
%
%
%   See also INTEGRAL_G, MEANMC_G, CUBMC_G
%
%  References
%
%   [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%   Yizhi Zhang, The Cost of Deterministic, Adaptive, Automatic
%   Algorithms: Cones, Not Balls, Journal of Complexity 30 (2014),
%   pp. 21-45.
%
%
%   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   and Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library
%   (Version 1.3)" [MATLAB Software], 2014. Available from
%   http://code.google.com/p/gail/
%
%   If you find GAIL helpful in your work, please support us by citing
%   the above paper and software.
%

% check parameter satisfy conditions or not
[f, out_param] = funappx_g_param(varargin{:});

MATLABVERSION= gail.matlab_version;
if MATLABVERSION >= 8.3
    warning('off', 'MATLAB:interp1:ppGriddedInterpolant');
end;

%% main algorithm

% initialize number of points
n = out_param.ninit;
% initialize nstar
out_param.nstar = n - 2;
% cost budget flag
out_param.exceedbudget = 1;
% tau change flag
tauchange = 0;
% length of interval
len = out_param.b-out_param.a;
% add flag
flag = 0;

while n < out_param.nmax;
    
    if(flag==0)
        x = out_param.a:len/(n-1):out_param.b;
        y = f(x);
    else
        xnew = repmat(x(1:end-1),m-1,1)+repmat((1:m-1)'*len/(n-1),1,(n-1)/m);
        ynew = f(xnew);
        xnew = [x(1:end-1); xnew];
        x = [xnew(:); x(end)]';
        ynew = [y(1:end-1); ynew];
        y = [ynew(:); y(end)]';
    end;
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = (n-1)/len*max(abs(diff_y-(y(n)-y(1))/(n-1)));
    %approximate the stronger norm of input function
    fn = (n-1)^2/len^2*max(abs(diff(diff_y)));
    
    % Stage 2: satisfy necessary condition
    if out_param.nstar*(2*gn+fn*len/(n-1)) >= fn*len;
        % Stage 3: check for convergence
        errbound = 4*out_param.abstol*(n-1)*(n-1-out_param.nstar)...
            /out_param.nstar/len;
        % satisfy convergence
        if errbound >= gn;
            out_param.exceedbudget = 0; break;
        end;
        % otherwise increase number of points
        m = max(ceil(1/(n-1)*sqrt(gn*out_param.nstar*...
            len/4/out_param.abstol)),2);
        n = m*(n-1)+1;
        flag = 1;
        % Stage2: do not satisfy necessary condition
    else
        % increase tau
        out_param.nstar = fn/(2*gn/len+fn/(n-1));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= out_param.nstar+2;
            % true, go to Stage 3
            errbound = 4*out_param.abstol*(n-1)*(n-1-out_param.nstar)...
                /out_param.nstar/len;
            if errbound >= gn;
                out_param.exceedbudget = 0; break;
            end;
            m = max(ceil(1/(n-1)*sqrt(gn*out_param.nstar*...
                len/4/out_param.abstol)),2);
            n = m*(n-1)+1;
            flag = 1;
        else
            % otherwise increase number of points, go to Stage 1
            n = 2 + ceil(out_param.nstar);
        end;
    end;
end;

if tauchange == 1;
    warning('MATLAB:funappx_g:peaky','This function is peaky relative to nlo and nhi. You may wish to increase nlo and nhi for similar functions.')
end;

% Check cost budget flag
if out_param.exceedbudget == 1;
    n = 1 + (n-1)/m*floor((out_param.nmax-1)*m/(n-1));
    warning('MATLAB:funappx_g:exceedbudget','funappx_g attempted to exceed the cost budget. The answer may be unreliable.')
    out_param.npoints = n;
    out_param.errorbound = fn*len^2/(8*(n-1)^2);
    %out_param.errbound = fn/(8*(n-1)^2);
    %out_param.ballradius = 2*out_param.abstol*(out_param.nmax-2)*(out_param.nmax...
    %    -2-out_param.tau)/out_param.tau;
    x1 = out_param.a:len/(out_param.npoints-1):out_param.b;
    y1 = f(x1);
    pp = interp1(x1,y1,'linear','pp');
else
    out_param.npoints = n;
    out_param.errorbound = fn*len^2/(8*(n-1)^2);
    pp = interp1(x,y,'linear','pp');
    
    
end;
if MATLABVERSION >= 8.3
    warning('on', 'MATLAB:interp1:ppGriddedInterpolant');
end;


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.nlo = 10;
default.nhi = 1000;
default.nmax  = 1e7;


if isempty(varargin)
    warning('MATLAB:funappx_g:nofunction','Function f must be specified. Now GAIL is using f(x)=x^2 and unit interval [0,1].')
    help funappx_g
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
    addRequired(p,'f',@gail.isfcn);
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

% let end point of interval not be infinity

% flagv = 0;
% if (length(out_param.a)==2)
%     flagv = 1;
%     if (length(varargin)==5)
%         out_param.nhi = out_param.nlo;
%         out_param.nlo = out_param.abstol;
%         out_param.abstol = out_param.b;
%         out_param.b = out_param.a(2);
%         out_param.a = out_param.a(1);
%     elseif(length(varargin)==6)
%         out_param.nmax = out_param.nhi;
%         out_param.nhi = out_param.nlo;
%         out_param.nlo = out_param.abstol;
%         out_param.abstol = out_param.b;
%         out_param.b = out_param.a(2);
%         out_param.a = out_param.a(1);
%     elseif(length(varargin)==4)
%         out_param.nlo = out_param.abstol;
%         out_param.abstol = out_param.b;
%         out_param.b = out_param.a(2);
%         out_param.a = out_param.a(1);
%     elseif(length(varargin)==3)
%         out_param.abstol = out_param.b;
%         out_param.b = out_param.a(2);
%         out_param.a = out_param.a(1);
%     elseif(length(varargin)==2)
%         out_param.b = out_param.a(2);
%         out_param.a = out_param.a(1);
%     end;
% end;


if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('MATLAB:funappx_g:anoinfinity',['a can not be infinity. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('MATLAB:funappx_g:bnoinfinity',['b can not be infinity. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('MATLAB:funappx_g:blea','b can not be smaller than a; exchange these two. ')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('MATLAB:funappx_g:beqa',['b can not equal a. Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['MATLAB:funappx_g:abstolnonpos ','Error tolerance should be greater than 0.' ...
        ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.isposintive(out_param.nmax)
        warning('MATLAB:funappx_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:funappx_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

% let initial number of points be a positive integer

% if (length(out_param.nlo) == 2 && flagv == 1)
%     if(length(varargin) ==4)
%         out_param.nhi = out_param.nlo(2);
%         out_param.nlo = out_param.nlo(1);
%     else
%         out_param.nmax = out_param.nhi;
%         out_param.nhi = out_param.nlo(2);
%         out_param.nlo = out_param.nlo(1);
%     end;
% elseif (length(out_param.nlo) == 2 && flagv == 0)
%     if(length(varargin) == 5)
%         out_param.nhi = out_param.nlo(2);
%         out_param.nlo = out_param.nlo(1);
%     else
%         out_param.nmax = out_param.nhi;
%         out_param.nhi = out_param.nlo(2);
%         out_param.nlo = out_param.nlo(1);
%     end;
% end;

if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('MATLAB:funappx_g:lowinitnotint',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:funappx_g:lowinitlt3',[' Lower bound of initial number of points should be a positive integer greater than 3.' ...
            ' Using 3 as nlo'])
        out_param.nlo = 3;
    end
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('MATLAB:funappx_g:hiinitnotint',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:funappx_g:hiinitlt3',[' Upper bound of points should be a positive integer greater than 3.' ...
            ' Using default number of points ' int2str(default.nhi) ' as nhi' ])
        out_param.nhi = default.nhi;
    end
end

if (out_param.nlo > out_param.nhi)
    warning('MATLAB:funappx_g:logrhi', 'Lower bound of initial number of points is larger than upper bound of initial number of points; Use nhi as nlo')
    %     temp = out_param.nlo;
    %     out_param.nlo = out_param.nhi;
    %     out_param.nhi = temp;
    out_param.nhi = out_param.nlo;
end;
if (out_param.nlo > out_param.nmax)
    warning('MATLAB:funappx_g:logecost',['Lower bound of initial number of points should be smaller than cost budget.' ...
        ' Using ', num2str(ceil(out_param.nmax/2))])
    out_param.nlo = out_param.nmax/2;
end;
if (out_param.nhi > out_param.nmax)
    warning('MATLAB:funappx_g:higecost',['Upper bound of initial number of points should be smaller than cost budget.' ...
        ' Using ', num2str(out_param.nlo)])
    out_param.nhi = out_param.nlo;
end;

h = out_param.b - out_param.a;
out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+h)));

