function [fmin,out_param]=funminNoPenalty_g(varargin)
%funminNoPenalty_g 1-D guaranteed locally adaptive function optimization 
%   on [a,b]
%
%   fmin = FUNMINNOPENALTY_G(f) finds minimum value of function f on the default
%   interval [0,1] within the guaranteed absolute error tolerance of 1e-6
%   and the X tolerance of 1e-3. Default initial number of points is 100
%   and default cost budget is 1e7. Input f is a function handle.
%
%   fmin = FUNMINNOPENALTY_G(f,a,b,abstol,TolX) finds minimum value of
%   function f with ordered input parameters that define the finite
%   interval [a,b], a guaranteed absolute error tolerance abstol and a
%   guaranteed X tolerance TolX.
%
%   fmin = FUNMINNOPENALTY_G(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%   finds minimum value of function f on the interval [a,b] with a 
%   guaranteed absolute error tolerance abstol and a guaranteed X tolerance 
%   TolX. All five
%   field-value pairs are optional and can be supplied in different order.
%
%   fmin = FUNMINNOPENALTY_G(f,in_param) finds minimum value of function f on the
%   interval [in_param.a,in_param.b] with a guaranteed absolute error
%   tolerance in_param.abstol and a guaranteed X tolerance in_param.TolX.
%   If a field is not specified, the default value is used.
%
%   [fmin, out_param] = FUNMINNOPENALTY_G(f,...) returns minimum value fmin of
%   function f and an output structure out_param.
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
%     in_param.TolX --- guaranteed X tolerance, default value is 1e-3.
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
%     out_param.TolX --- guaranteed X tolerance
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
%     out_param.volumeX --- the volume of intervals containing the point(s)
%     where the minimum occurs
%
%     out_param.intervals --- the intervals containing point(s) where the
%     minimum occurs. Each column indicates one interval where the first
%     row is the left point and the second row is the right point.
%
%  Guarantee
%
%
%
%  Examples
%
%  Example 1:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2); [fmin,out_param] = funminNoPenalty_g(f)
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
%          TolX: 1.0000e-03
%           nlo: 10
%           nhi: 1000
%         ninit: 100
%          nmax: 10000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 1
%       npoints: 199
%        errest: 1.1306e-07
%       volumeX: 0.0101
%     intervals: [2x1 double]
%
%
%  Example 2:
%
%  >> f = @(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funminNoPenalty_g(f,-2,2,1e-7,1e-4,10,10,1000000)
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
%          TolX: 1.0000e-04
%           nlo: 10
%           nhi: 10
%         ninit: 10
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 7
%       npoints: 1153
%        errest: 3.9920e-08
%       volumeX: 0.0069
%     intervals: [2x1 double]
%
%
%  Example 3:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> in_param.a = -13; in_param.b = 8;
%  >> in_param.abstol = 10^(-7); in_param.TolX = 1e-4;
%  >> in_param.nlo = 10; in_param.nhi = 100;
%  >> in_param.nmax = 10^6;
%  >> [fmin,out_param] = funminNoPenalty_g(f,in_param)
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
%          TolX: 1.0000e-04
%           nlo: 10
%           nhi: 100
%         ninit: 91
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 8
%       npoints: 9508
%        errest: 9.9894e-08
%       volumeX: 0.0027
%     intervals: [2x1 double]
%
%
%  Example 4:
%
%  >> f=@(x) exp(0.01*(x-0.5).^2);
%  >> [fmin,out_param] = funminNoPenalty_g(f,'a',-2,'b',2,'nhi',100,'nlo',10,'nmax',1e6,'abstol',0,'TolX',1e-4)
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
%        abstol: 0
%          TolX: 1.0000e-04
%           nlo: 10
%           nhi: 100
%         ninit: 64
%          nmax: 1000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 11
%       npoints: 129025
%        errest: 3.1664e-12
%       volumeX: 6.2004e-05
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
[f, in_param] = funminNoPenalty_g_param(varargin{:});
MATLABVERSION = gail.matlab_version;
out_param = in_param;


%% main algorithm
a = out_param.a;
b = out_param.b;
abstol = out_param.abstol;
TolX = out_param.TolX;
n = out_param.ninit;
x = a:(b-a)/(n-1):b;
y = f(x);
iSing = find(isinf(y));
if ~isempty(iSing)
    error('GAIL:funminNoPenalty_g:yInf',['Function f(x) = Inf at x = ', num2str(x(iSing))]);
end
if length(y) == 1  
    % probably f is a constant function and Matlab would  
    % reutrn only a value fmin 
    fmin = y;
    max_errest = 0;
end

fh = 4*(b-a)/(n-1);
C0 = 2.2;
C = @(h) C0*fh./(fh-h);
max_errest = 1;

iter = 1;
exit_len = 2;
% we start the algorithm with all warning flags down
out_param.exitflag = false(1,exit_len); 


while n < out_param.nmax
    %% Stage 1: compute length of each subinterval and approximate |f''(t)|
    len = diff(x);
    deltaf = diff(diff(y)./len)./(len(2:end)+len(1:end-1));
    deltaf = [0 0 abs(deltaf) 0 0];
    
    %% Stage 2: compute bound of |f''(t)| and estimate error
    h = [x(2)-a  x(3)-a  x(4:end)-x(1:end-3)  b-x(end-2)  b-x(end-1)];
    normbd = C(max(h(1:n-1),h(3:n+1))) .* max(deltaf(1:n-1),deltaf(4:n+2));
    errest = len.^2/8.*normbd;
    % find intervals contain minimum
    Un=min(y);
    diff_y=diff(y);
    ln=diff_y/2+y(1:n-1)-abs(diff_y)/2-errest;
    goodindex = find(ln < Un);
    m = size(goodindex,2);
    if m > 0
        ints = zeros(2,m);
        ints = [x(goodindex); x(goodindex+1)];
        leftint = find([true diff(goodindex)~=1]);
        rightint = find([diff(goodindex)~=1 true]);
        q = size(leftint,2);
        ints1 = zeros(2,q);
        ints1(1,:) = ints(1,leftint);
        ints1(2,:) = ints(2,rightint);
    else
        ints1 = zeros(2,0);
    end
    interval=ints1;
    volumeX = sum(interval(2,:)-interval(1,:));
    % update iterations
    max_errest = max(errest);
    if max_errest <= abstol || volumeX < TolX
        break
    end 
    iter = iter + 1;
 
    %% Stage 3: find I and update x,y
    badinterval = (errest > abstol | ln < Un);
    whichcut = badinterval | [badinterval(2:end) 0] | [0 badinterval(1:end-1)];
    if (out_param.nmax<(n+length(find(whichcut==1))))
        out_param.exitflag(1) = true;
        warning('GAIL:funminNoPenalty_g:exceedbudget',['funminNoPenalty_g'...
            'attempted to exceed the cost budget. The answer may be '...
            'unreliable.'])
        break;
    end; 
    if(iter==out_param.maxiter)
        out_param.exitflag(2) = true;
        warning('GAIL:funminNoPenalty_g:exceediter',['Number of iterations has '...
            'reached maximum number of iterations.'])
        break;
    end;
    newx = x(whichcut) + 0.5 * len(whichcut);
    tt = cumsum(whichcut); 
    x([1 (2:n)+tt]) = x;
    y([1 (2:n)+tt]) = y;
    tem = 2 * tt + cumsum(whichcut==0);
    x(tem(whichcut)) = newx;
    y(tem(whichcut)) = f(newx);
    n = length(x);
end;

%% postprocessing
fmin = Un; 
out_param.iter = iter;
out_param.npoints = n;
out_param.errest = max_errest;
out_param.volumeX = volumeX;
out_param.intervals = interval;

% control the order of out_param
out_param = orderfields(out_param, ...
            {'f', 'a', 'b','abstol','TolX','nlo','nhi','ninit','nmax','maxiter',...
             'exitflag','iter','npoints','errest','volumeX', 'intervals'});




function [f, out_param] = funminNoPenalty_g_param(varargin)
% parse the input to the funminNoPenalty_g function

%% Default parameter values
default.a = 0;
default.b = 1;
default.abstol = 1e-6;
default.TolX = 1e-3;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;
default.maxiter = 1000;



MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;

 
if isempty(varargin)
  warning('GAIL:funminNoPenalty_g:nofunction',['Function f must be specified. '...
      'Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval '...
      '[0,1].'])
  help funminNoPenalty_g
  f = @(x) exp(0.01*(x-0.5).^2);
  out_param.f = f;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:funminNoPenalty_g:notfunction',['Function f must be a '...
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
    out_param.TolX = default.TolX;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.nmax = default.nmax ;
    out_param.maxiter = default.maxiter;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'TolX',default.TolX,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric)
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'a',default.a,@isnumeric);
        f_addParamVal(p,'b',default.b,@isnumeric);
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'TolX',default.TolX,@isnumeric);
        f_addParamVal(p,'nlo',default.nlo,@isnumeric);
        f_addParamVal(p,'nhi',default.nhi,@isnumeric);
        f_addParamVal(p,'nmax',default.nmax,@isnumeric);
        f_addParamVal(p,'maxiter',default.maxiter,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let end point of interval not be infinity
if (out_param.a == inf||out_param.a == -inf)
    warning('GAIL:funminNoPenalty_g:aisinf',['a cannot be infinity. '...
        'Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf)
    warning(['GAIL:funminNoPenalty_g:bisinf','b cannot be infinity. '...
        'Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('GAIL:funminNoPenalty_g:blea',['b cannot be smaller than a;'...
        ' exchange these two. '])
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('GAIL:funminNoPenalty_g:beqa',['b cannot equal a. '...
        'Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% let error tolerance greater than 0
if (out_param.abstol < 0 )
    warning('GAIL:funminNoPenalty_g:tolneg', ['Error tolerance should be greater'...
        ' than 0. Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

%let length tolerance greater than 0
if (out_param.TolX < 0)
    warning('GAIL:funmin_g:Xtolnonpos', ['X tolerance should be greater than or equal to 0.' ...
        ' Using default X tolerance ' num2str(default.TolX)]);
    out_param.TolX = default.TolX;
end

% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.isposintive(out_param.nmax)
        warning('GAIL:funminNoPenalty_g:budgetnotint',['Cost budget should be '...
            'a positive integer. Using cost budget '...
            , num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:funminNoPenalty_g:budgetisneg',['Cost budget should be '...
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
    warning('GAIL:funminNoPenalty_g:logrhi', ['Lower bound of initial number of'...
        ' points is larger than upper bound of initial number of '...
        'points; Use nhi as nlo'])
    out_param.nhi = out_param.nlo;
end;

h = out_param.b - out_param.a;
out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)...
    ^(1/(1+h)));

if (~gail.isposint(out_param.maxiter))
    if gail.ispositive(out_param.maxiter)
        warning('GAIL:funminNoPenalty_g:maxiternotint',['Max number of '...
            'iterations should be a positive integer. Using max number '...
            'of iterations as  ', num2str(ceil(out_param.maxiter))])
        out_param.maxiter = ceil(out_param.maxiter);
    else
        warning('GAIL:funminNoPenalty_g:budgetisneg',['Max number of iterations'...
            ' should be a positive integer. Using max number of '...
            'iterations as ' int2str(default.maxiter)])
        out_param.maxiter = default.maxiter;
    end;
end

