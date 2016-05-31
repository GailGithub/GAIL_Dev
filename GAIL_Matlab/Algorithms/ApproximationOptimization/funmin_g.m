function [fmin,out_param]=funmin_g(varargin)
%funmin_g 1-D guaranteed locally adaptive function optimization
%   on [a,b]
%
%   fmin = FUNMIN_G(f) finds minimum value of function f on the  
%   default interval [0,1] within the guaranteed absolute error tolerance 
%   of 1e-6. Input f is a function handle.
%
%   fmin = FUNMIN_G(f,a,b,abstol) finds minimum value of
%   function f with ordered input parameters that define the finite
%   interval [a,b], and a guaranteed absolute error tolerance abstol.
%
%   fmin = FUNMIN_G(f,'a',a,'b',b,'abstol',abstol) finds minimum
%   value of function f on the interval [a,b] with a guaranteed absolute
%   error tolerance. All four field-value pairs are optional and can be
%   supplied in different order.
%
%   fmin = FUNMIN_G(f,in_param) finds minimum value of function f  
%   on the interval [in_param.a,in_param.b] with a guaranteed absolute
%   error tolerance in_param.abstol. If a field is not specified, the
%   default value is used.
%
%   [fmin, out_param] = FUNMIN_G(f,...) returns minimum value fmin
%   of function f and an output structure out_param.
%
%   Input Arguments
%
%     f --- input function
%
%     in_param.a --- left end point of interval, default value is 0
%
%     in_param.b --- right end point of interval, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6.
%
%
%   Optional Input Arguments
%
%     in_param.nlo --- lower bound of initial number of points we used,
%     default value is 10
%
%     in_param.nhi --- upper bound of initial number of points we used,
%     default value is 1000
%
%     in_param.nmax --- cost budget, default value is 1e7.
%
%   Output Arguments
%
%     out_param.f --- input function
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
%     out_param.ninit --- initial number of points we use
%
%     out_param.npoints --- number of points needed to reach the guaranteed
%     absolute error tolerance or the guaranteed X tolerance
%
%     out_param.exitflag --- the state of program when exiting
%              0  Success
%              1  Number of points used is greater than out_param.nmax
%
%     out_param.errest --- estimation of the absolute error bound
%
%     out_param.intervals --- the intervals containing point(s) where the
%     minimum occurs. Each column indicates one interval where the first
%     row is the left point and the second row is the right point.
%
%  Guarantee
%
%  For [a,b] there exists a partition, P={[t_0,t_1], [t_1,t_2], ...,
%  [t_{L-1},t_L]}, where a=t_0 < t_1 < ... < t_L=b. If the function to be
%  minimized, f, satisfies the cone condition
%  for each sub interval [t_{l-1},t_l], where 1 <= l <= L, then the output
%  fappx by this algorithm is guaranteed to satisfy
%      |
%
%
%  Examples
%
%  Example 1:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2); [fmin,out_param] = funmin_g(f)
%
%  fmin =
%
%      1
% 
%  out_param = 
% 
%             f: @(x)exp(0.01*(x-0.5).^2)
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%           nlo: 10
%           nhi: 1000
%         ninit: 100
%          nmax: 10000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 2
%       npoints: 123
%        errest: 3.7879e-07
%     intervals: [2x1 double]
%
%
%  Example 2:
%
%  >> f = @(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funmin_g(f,-2,2,1e-7,10,10,1000000)
%
%  fmin =
% 
%      1
% 
%  out_param = 
% 
%             f: @(x)exp(0.01*(x-0.5).^2)
%             a: -2
%             b: 2
%        abstol: 1.0000e-07
%           nlo: 10
%           nhi: 10
%         ninit: 10
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 9
%       npoints: 57
%        errest: 9.1055e-08
%     intervals: [2x1 double]
%
%
%  Example 3:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> in_param.a = -13; in_param.b = 8;
%  >> in_param.abstol = 10^(-7);
%  >> in_param.nlo = 10; in_param.nhi = 100;
%  >> in_param.nmax = 10^6;
%  >> [fmin,out_param] = funmin_g(f,in_param)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)exp(0.01*(x-0.5).^2)
%             a: -13
%             b: 8
%        abstol: 1.0000e-07
%           nlo: 10
%           nhi: 100
%         ninit: 91
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 9
%       npoints: 149
%        errest: 2.5117e-08
%     intervals: [2x1 double]
%
%
%  Example 4:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funmin_g(f,'a',-2,'b',2,'nhi',100,'nlo',10,'nmax',1e6,'abstol',1e-5)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)exp(0.01*(x-0.5).^2)
%             a: -2
%             b: 2
%        abstol: 1.0000e-05
%           nlo: 10
%           nhi: 100
%         ninit: 64
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 4
%       npoints: 88
%        errest: 2.5064e-06
%     intervals: [2x1 double]
%
%
%  See also FMINBND, FUNAPPX_G, INTEGRAL_G
%
%  References
%   [1]  Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
%   Univariate Function Minimization," MS thesis, Illinois Institute of
%   Technology, 2014.
%
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
%   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
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
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%


% check parameter satisfy conditions or not
[f, in_param] = funmin_g_param(varargin{:});
MATLABVERSION = gail.matlab_version;
out_param = in_param;


%% main algorithm
a = out_param.a;
b = out_param.b;
abstol = out_param.abstol;
n = out_param.ninit;
x = a:(b-a)/(n-1):b;
y = f(x);
iSing = find(isinf(y));
if ~isempty(iSing)
    error('GAIL:funmin_g:yInf',['Function f(x) = Inf at x = ', num2str(x(iSing))]);
end
if length(y) == 1
    % probably f is a constant function and Matlab would
    % reutrn only a value fmin
    fmin = y;
    max_errest = 0;
end
iter = 0;
exit_len = 2;
fh = 3*(b-a)/(n-2);
C0 = 3;
C = @(h) (C0*fh)./(fh-h);
max_errest = 1;

% we start the algorithm with all warning flags down
out_param.exitflag = false(1,exit_len);


while n < out_param.nmax
    %% Stage 1: compute length of each subinterval and approximate |f''(t)|
    len = diff(x(1:n));
    deltaf = 2 * diff(diff(y(1:n))./len) ./ (len(1:end-1) + len(2:end));
    h = x(4:n)-x(1:n-3);
    Br = [abs(deltaf(2:end)).*C(h) 0 0];
    Bl = [0 0 abs(deltaf(1:end-1)).*C(h)];
    Un=min(y);
    diff_y=diff(y);
    min_int = (y(1:n-1)+y(2:n)-abs(diff_y))./2;
    errest = Un+len.^2/8.*max(Br,Bl)-min_int;
    max_errest = max(errest);
    lowerbound = Un-errest;
    
    
    
    %% Stage 2: compute bound of |f''(t)| and estimate error
    
    % update iterations
    iter = iter + 1;
    if max_errest <= abstol,
        break
    end
    
    
    %% Stage 3: find I and update x,y
    badinterval = (errest > abstol);
    badlinterval= (Un-min_int+len.^2/8.*Bl>abstol);
    badrinterval= (Un-min_int+len.^2/8.*Br>abstol);
    maybecut=(badinterval|[0 badlinterval(3:end) 0]|[badlinterval(3:end)...
        0 0]|[0 badrinterval(1:end-2) 0]|[0 0 badrinterval(1:end-2)]);
    maxlength = (len>max(len(maybecut))-eps);
    whichcut = maybecut & maxlength;
    if (out_param.nmax<(n+length(find(whichcut))))
        out_param.exitflag(1) = true;
        warning('GAIL:funmin_g:exceedbudget',['funmin_g '...
            'attempted to exceed the cost budget. The answer may be '...
            'unreliable.'])
        break;
    end;
    if(iter==out_param.maxiter)
        out_param.exitflag(2) = true;
        warning('GAIL:funmin_g:exceediter',['Number of iterations has '...
            'reached maximum number of iterations.'])
        break;
    end;
    newx = x(whichcut) + 0.5 * len(whichcut);
    if n + length(newx) > length(x)
        xx(1:n) = x(1:n);
        yy(1:n) = y(1:n);
        x = xx;
        y = yy;
    end
    tt = cumsum(whichcut);
    x([1 (2:n)+tt]) = x(1:n);
    y([1 (2:n)+tt]) = y(1:n);
    tem = 2 * tt + cumsum(whichcut==0);
    x(tem(whichcut)) = newx;
    y(tem(whichcut)) = f(newx);
    n = n + length(newx);
end;

%% find intervals
index = find(lowerbound < Un);
m = size(index,2);
if m > 0
    ints = zeros(2,m);
    ints = [x(index); x(index+1)];
    leftint = find([true diff(index)~=1]);
    rightint = find([diff(index)~=1 true]);
    q = size(leftint,2);
    ints1 = zeros(2,q);
    ints1(1,:) = ints(1,leftint);
    ints1(2,:) = ints(2,rightint);
else
    ints1 = zeros(2,0);
end
interval=ints1;


%% postprocessing
fmin = Un;
out_param.iter = iter;
out_param.npoints = n;
out_param.errest = max_errest;
out_param.intervals = interval;

% control the order of out_param
out_param = orderfields(out_param, ...
{'f', 'a', 'b','abstol','nlo','nhi','ninit','nmax','maxiter',...
'exitflag','iter','npoints','errest','intervals','output_x'});

if (in_param.output_x)
  out_param.x = x;
  out_param.y = y;
end

function [f, out_param] = funmin_g_param(varargin)
% parse the input to the funmin_g function

%% Default parameter values
default.a = 0;
default.b = 1;
default.abstol = 1e-6;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;
default.maxiter = 1000;
default.output_x = 0;

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;


if isempty(varargin)
    warning('GAIL:funmin_g:nofunction',['Function f must be specified. '...
        'Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval '...
        '[0,1].'])
    help funmin_g
    f = @(x) exp(0.01*(x-0.5).^2);
    out_param.f = f;
else
    if gail.isfcn(varargin{1})
        f = varargin{1};
        out_param.f = f;
    else
        warning('GAIL:funmin_g:notfunction',['Function f must be a '...
            'function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2).'])
        f = @(x) exp(0.01*(x-0.5).^2);
        out_param.f = f;
    end
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
    out_param.abstol = default.abstol;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.nmax = default.nmax ;
    out_param.maxiter = default.maxiter;
    out_param.output_x = default.output_x;
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
        addOptional(p,'nmax',default.nmax,@isnumeric)
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
        addOptional(p,'output_x',default.output_x,@isnumeric)
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'a',default.a,@isnumeric);
        f_addParamVal(p,'b',default.b,@isnumeric);
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'nlo',default.nlo,@isnumeric);
        f_addParamVal(p,'nhi',default.nhi,@isnumeric);
        f_addParamVal(p,'nmax',default.nmax,@isnumeric);
        f_addParamVal(p,'maxiter',default.maxiter,@isnumeric);
        f_addParamVal(p,'output_x',default.output_x,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let end point of interval not be infinity
if (out_param.a == inf || out_param.a == -inf)
    warning('GAIL:funmin_g:aisinf',['a cannot be infinity. '...
        'Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf || out_param.b == -inf)
    warning(['GAIL:funmin_g:bisinf','b cannot be infinity. '...
        'Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('GAIL:funmin_g:blea',['b cannot be smaller than a;'...
        ' exchange these two. '])
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('GAIL:funmin_g:beqa',['b cannot equal a. '...
        'Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% let error tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:funmin_g:tolneg', ['Error tolerance should be greater'...
        ' than 0. Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.isposintive(out_param.nmax)
        warning('GAIL:funmin_g:budgetnotint',['Cost budget should be '...
            'a positive integer. Using cost budget '...
            , num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:funmin_g:budgetisneg',['Cost budget should be '...
            'a positive integer. Using default cost budget '...
            int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('GAIL:funminglobal_g:lowinitnotint',['Lower bound of '...
            'initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('GAIL:funminglobal_g:lowinitlt3',[' Lower bound of '...
            'initial number of points should be a positive integer greater'...
            ' than 3. Using 3 as nlo'])
        out_param.nlo = 3;
    end
    warning('GAIL:funminglobal_g:lowinitnotint',['Lower bound of '...
        'initial nstar should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
    out_param.nlo = ceil(out_param.nlo);
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('GAIL:funminglobal_g:hiinitnotint',['Upper bound of '...
            'initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('GAIL:funminglobal_g:hiinitlt3',[' Upper bound of '...
            'points should be a positive integer greater than 3. Using '...
            'default number of points ' int2str(default.nhi) ' as nhi' ])
        out_param.nhi = default.nhi;
    end
    warning('GAIL:funminglobal_g:hiinitnotint',['Upper bound of '...
        'initial nstar should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
    out_param.nhi = ceil(out_param.nhi);
end

if (out_param.nlo > out_param.nhi)
    warning('GAIL:funmin_g:logrhi', ['Lower bound of initial number of'...
        ' points is larger than upper bound of initial number of '...
        'points; Use nhi as nlo'])
    out_param.nhi = out_param.nlo;
end;

h = out_param.b - out_param.a;
out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)...
    ^(1/(1+h)));

if (~gail.isposint(out_param.maxiter))
    if gail.ispositive(out_param.maxiter)
        warning('GAIL:funmin_g:maxiternotint',['Max number of '...
            'iterations should be a positive integer. Using max number '...
            'of iterations as  ', num2str(ceil(out_param.maxiter))])
        out_param.maxiter = ceil(out_param.maxiter);
    else
        warning('GAIL:funmin_g:budgetisneg',['Max number of iterations'...
            ' should be a positive integer. Using max number of '...
            'iterations as ' int2str(default.maxiter)])
        out_param.maxiter = default.maxiter;
    end;
end

