function [fmin,out_param]=funmin_g(varargin)
%funmin_g 1-D guaranteed locally adaptive function optimization on [a,b]
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
%   error tolerance. All three field-value pairs are optional and can be
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
%     in_param.ninit --- initial number of subintervals. Default to 20.
%
%     in_param.nmax --- cost budget, default value is 1e7.
%
%     in_param.maxiter --- max number of iterations, default value is 1000.
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
%     out_param.nmax --- cost budget
%
%     out_param.ninit --- initial number of subintervals
%
%     out_param.npoints --- number of points needed to reach the guaranteed
%     absolute error tolerance
%
%     out_param.exit --- this is a vector with two elements, for
%     tracking important warnings in the algorithm. The algorithm is considered successful (with
%     out_param.exit == [0 0]) if no flags arise warning that the
%     results are not guaranteed. The initial value is [0 0] and
%     the final value of this parameter is encoded as follows:
%
%                      [1 0]:   If reaching overbudget. It states whether
%                      the max budget is attained without reaching the
%                      guaranteed error tolerance.
%
%                      [0 1]:   If reaching overiteration. It states whether
%                      the max iterations is attained without reaching the
%                      guaranteed error tolerance.
%
%     out_param.errest --- estimation of the absolute error bound
%
%     out_param.iter --- number of iterations
%
%     out_param.intervals --- the intervals containing point(s) where the
%     minimum occurs. Each column indicates one interval where the first
%     row is the left point and the second row is the right point.
%
%
%  Guarantee
%
%  Please check the details of the guarantee in [1].
%
%
%  Examples
%
%  Example 1:
%
%  >> f = @(x) exp(0.01*(x-0.5).^2); [fmin,out_param] = funmin_g(f)
%
%  fmin =
%
%      1
%
%  out_param =***
%
%             f: @(x)exp(0.01*(x-0.5).^2)
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%         ninit: 20
%          nmax: 10000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 5
%       npoints: 69
%        errest: 2.5955e-07
%     intervals: [2***1 double]
%
%
%  Example 2:
%
%  >> f = @(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funmin_g(f,-2,2,1e-7,10,1000000)
%
%  fmin =
%
%      1
%
%  out_param =***
%
%            f: @(x)exp(0.01*(x-0.5).^2)
%            a: -2
%            b: 2
%       abstol: 1.0000e-07
%        ninit: 10
%         nmax: 1000000
%      maxiter: 1000
%     exitflag: [0 0]
%         iter: 9
%      npoints: 79
%       errest: 6.1251e-08
%    intervals: [2***1 double]
%
%
%  Example 3:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> in_param.a = -13; in_param.b = 8;
%  >> in_param.abstol = 10^(-7);
%  >> in_param.ninit = 100;
%  >> in_param.nmax = 10^6;
%  >> [fmin,out_param] = funmin_g(f,in_param)
%
%  fmin =
%
%     1.0000
%
%  out_param =***
%
%            f: @(x)exp(0.01*(x-0.5).^2)
%            a: -13
%            b: 8
%       abstol: 1.0000e-07
%        ninit: 100
%         nmax: 1000000
%      maxiter: 1000
%     exitflag: [0 0]
%         iter: 8
%      npoints: 203
%       errest: 6.7816e-08
%    intervals: [2***1 double]
%
%
%  Example 4:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funmin_g(f,'a',-2,'b',2,'ninit',64,'nmax',1e6,'abstol',1e-5)
%
%  fmin =
%
%     1
%
%  out_param =***
%
%            f: @(x)exp(0.01*(x-0.5).^2)
%            a: -2
%            b: 2
%       abstol: 1.0000e-05
%        ninit: 64
%         nmax: 1000000
%      maxiter: 1000
%     exitflag: [0 0]
%         iter: 3
%      npoints: 107
%       errest: 8.0997e-06
%    intervals: [2***1 double]
%
% >> out_param(:).intervals
% ans =
%
%    0.3594
%    0.6406
%
%
%  See also FMINBND, FUNAPPX_G, INTEGRAL_G
%
%
%  References
%
%   [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Xin Tong, "Local
%   Adaption for Approximation and Minimization of Univariate Functions,"
%   Journal of Complexity 40, pp. 17-33, 2017.
%
%   [2] Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
%   Univariate Function Minimization," MS thesis, Illinois Institute of
%   Technology, 2014.
%
%   [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   [4] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%


% check parameter satisfy conditions or not
[f, in_param] = funmin_g_param(varargin{:});
out_param = in_param;
out_param = rmfield(out_param,'output_x');

%% main algorithm
a = out_param.a;
b = out_param.b;
abstol = out_param.abstol;
n = out_param.ninit+1;
x = zeros(1, max(100,ceil(out_param.nmax/100))); % preallocation
y = x;
x(1:n) = a:(b-a)/(n-1):b;
y(1:n) = f(x(1:n));
iSing = find(isinf(y));
if ~isempty(iSing)
    out_param.exitflag(2) = true;
    error('GAIL:funmin_g:yInf',['Function f(x) = Inf at x = ', num2str(x(iSing))]);
end
if length(y) == 1
    % probably f is a constant function and Matlab would
    % reutrn only a value fmin
    M_hat = y;
    %max_errest = 0;
end
iter = 0;
exit_len = 2;
% we start the algorithm with all warning flags down
out_param.exitflag = false(1,exit_len);

fh = 3*(b-a)/(n-2);
C0 = 10;
C = @(h) (C0*fh)./(fh-h);
indexI_p = ([0 0 ones(1,n-3) 0]>0);
indexI_m = ([0 ones(1,n-3) 0 0]>0);
max_errest = 1;

while n < out_param.nmax
    %% Stage 1: compute length of each subinterval and approximate |f''(t)|
    len = diff(x(1:n));
    deltaf = diff(diff(y(1:n)));
    h = x(2:n-1) - x(1:n-2);
    diff_y=diff(y);
    min_int = (y(1:n-1)+y(2:n)-abs(diff_y(1:n-1)))./2;
    M_hat=min(y(1:n));

    err_p=zeros(1,n); err_hat_p=zeros(1,n);
    err_m=zeros(1,n); err_hat_m=zeros(1,n);
    err_p(2:end-1) = abs(1/8 * C(3*h)...
                           .* deltaf.*indexI_p(2:end-1));
    indexI_tilde_p = (err_p > abstol);
    err_m(2:end-1) = abs(1/8 * C(3*h)...
                           .* deltaf.*indexI_m(2:end-1));
    indexI_tilde_m = (err_m > abstol);
    err_hat_p(3:end-1) = indexI_p(3:end-1).*...
        (err_p(3:end-1)+M_hat-min_int(1:end-2));
    err_hat_m(2:end-2) = indexI_m(2:end-2).*...
        (err_m(2:end-2)+M_hat-min_int(3:end));

    indexI_hat_p = (err_hat_p > abstol)| ...
        ((indexI_tilde_p.*[0 0 0 indexI_tilde_m(1:end-3)]) & ...
        (([indexI_tilde_p(4:end) 0 0 0].*err_hat_m)>abstol));
    indexI_hat_m = (err_hat_m > abstol)| ...
        ((indexI_tilde_m.*[indexI_tilde_p(4:end) 0 0 0]) & ...
        (([0 0 0 indexI_tilde_m(1:end-3)].*err_hat_p)>abstol));

    % update iterations
    iter = iter + 1;
    if max(indexI_hat_p|indexI_hat_m) == 0,
        max_errest = max(max(err_p),max(err_m));
        break;
    end

    %% Stage 2: Split the subintervals as needed
    %find the index of the subinterval which is needed to be cut
    midpoint_p = [indexI_hat_p(3:end) 0 0 ] | [indexI_hat_p(2:end) 0];
    midpoint_m = indexI_hat_m | [0 indexI_hat_m(1:end-1)];
    whichcut = midpoint_p(1:end-1) | midpoint_m(1:end-1);

    %check to see if exceed the cost budget
    if (out_param.nmax<(n+length(find(whichcut))))
        out_param.exitflag(1) = true;
        warning('GAIL:funmin_g:exceedbudget',['funmin_g '...
            'attempted to exceed the cost budget. The answer may be '...
            'unreliable.'])
        break;
    end

    %check to see if exceed the maximumber number of iterations
    if(iter==out_param.maxiter)
        out_param.exitflag(2) = true;
        warning('GAIL:funmin_g:exceediter',['Number of iterations has '...
            'reached maximum number of iterations.'])
        break;
    end

    %generate split points for x
    newx=x(whichcut)+0.5*len(whichcut);

    %relocate the space for new x
    if n + length(newx) > length(x)
      xx = zeros(1, out_param.nmax);
      yy = xx;
      xx(1:n) = x(1:n);
      yy(1:n) = y(1:n);
      x = xx;
      y = yy;
    end

    %update x and y
    tt = cumsum(whichcut);
    x([1 (2:n)+tt]) = x(1:n);
    y([1 (2:n)+tt]) = y(1:n);
    tem = 2 * tt + cumsum(whichcut==0);
    x(tem(whichcut)) = newx;
    y(tem(whichcut)) = f(newx);

    %update the set I to consist of the new indices
    newindex_p([1 (2:n)+tt]) = [indexI_hat_p(2:end) 0];
    newindex_p(tem) = indexI_hat_p(2:end);
    indexI_p = ([0 0 newindex_p(3:end-1) 0]>0);
    newindex_m([1 (2:n)+tt]) = [0 indexI_hat_m(1:end-1)];
    newindex_m(tem) = indexI_hat_m(1:end-1);
    indexI_m = ([0 newindex_m(2:end-2) 0 0]>0);

    %update # of points
    n = n + length(newx);
end

%% postprocessing
fmin = M_hat;
out_param.iter = iter;
out_param.npoints = n;
out_param.errest = max_errest;


% control the order of out_param
out_param = orderfields(out_param, ...
{'f', 'a', 'b','abstol','ninit','nmax','maxiter',...
'exitflag','iter','npoints','errest'});

ints = find(ismember( x, newx ));
leftint = ints(logical([true diff(ints)>2]));
rightint = ints(logical([diff(ints)>2 true]));
q = size(leftint,2);
ints1 = zeros(2,q);
ints1(1,:) = x(leftint);
ints1(2,:) = x(rightint);
out_param.intervals = ints1;

if (in_param.output_x)
   out_param.x = x(1:n);
   out_param.y = y(1:n);
end

function [f, out_param] = funmin_g_param(varargin)
% parse the input to the funmin_g function

%% Default parameter values
default.a = 0;
default.b = 1;
default.abstol = 1e-6;
default.ninit = 20;
default.nmax = 1e7;
default.maxiter = 1000;
default.output_x = false;

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
end

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
    out_param.ninit = default.ninit;
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
        addOptional(p,'ninit',default.ninit,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
        addOptional(p,'maxiter',default.maxiter,@isnumeric);
        addOptional(p,'output_x',default.output_x,@logical);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'a',default.a,@isnumeric);
        f_addParamVal(p,'b',default.b,@isnumeric);
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'ninit',default.ninit,@isnumeric);
        f_addParamVal(p,'nmax',default.nmax,@isnumeric);
        f_addParamVal(p,'maxiter',default.maxiter,@isnumeric);
        f_addParamVal(p,'output_x',default.output_x,@logical);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end

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
    end
end

if (~gail.isposint(out_param.ninit))
    if gail.isposge3(out_param.ninit)
        warning('GAIL:funming_g:lowinitnotint',['Lower bound of '...
            'initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.ninit)) ' as ninit '])
        out_param.ninit = ceil(out_param.ninit);
    else
        warning('GAIL:funmin_g:lowinitlt3',[' Lower bound of '...
            'initial number of points should be a positive integer greater'...
            ' than 3. Using 3 as ninit'])
        out_param.ninit = 3;
    end
    warning('GAIL:funmin_g:lowinitnotint',['Lower bound of '...
        'initial nstar should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.ninit)) ' as ninit '])
    out_param.ninit = ceil(out_param.ninit);
end

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
    end
end
if (out_param.output_x~=true&&out_param.output_x~=false)
    warning('GAIL:funmin_g:output_x', ['Input of output_x'...
        ' can only be true or false; use default value false'])
    out_param.output_x = false;
end
