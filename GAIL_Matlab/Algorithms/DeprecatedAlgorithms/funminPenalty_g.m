function [fmin,out_param]=funminPenalty_g(varargin)
%FUNMINPENALTY_G 1-D guaranteed global minimum value on [a,b] and the subset
%containing optimal solutions
%
%   fmin = FUNMINPENALTY_G(f) finds minimum value of function f on the default
%   interval [0,1] within the guaranteed absolute error tolerance of 1e-6
%   and the X tolerance of 1e-3. Default initial number of points is 100
%   and default cost budget is 1e7. Input f is a function handle.
%
%   fmin = funmin_G(f,a,b,abstol,TolX) finds minimum value of
%   function f with ordered input parameters that define the finite
%   interval [a,b], a guaranteed absolute error tolerance abstol and a
%   guaranteed X tolerance TolX.
%
%   fmin = FUNMINPENALTY_G(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%   finds minimum value of function f on the interval [a,b] with a 
%   guaranteed absolute error tolerance abstol and a guaranteed X tolerance 
%   TolX. All five
%   field-value pairs are optional and can be supplied in different order.
%
%   fmin = FUNMINPENALTY_G(f,in_param) finds minimum value of function f on the
%   interval [in_param.a,in_param.b] with a guaranteed absolute error
%   tolerance in_param.abstol and a guaranteed X tolerance in_param.TolX.
%   If a field is not specified, the default value is used.
%
%   [fmin, out_param] = FUNMINPENALTY_G(f,...) returns minimum value fmin of
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
%     out_param.tau --- latest value of tau
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
%     out_param.tauchange --- it is 1 if out_param.tau changes, otherwise
%     it is 0
%
%     out_param.intervals --- the intervals containing point(s) where the
%     minimum occurs. Each column indicates one interval where the first
%     row is the left point and the second row is the right point.
%
%  Guarantee
%
%   If the function to be minimized, f satisfies the cone condition
%
%   ||f''||_\infty <=  tau/(b-a)||f'-(f(b)-f(a))/(b-a)||__\infty,
%      
%   then the fmin output by this algorithm is guaranteed to satisfy
%
%       | min(f)-fmin| <= abstol,
%   or
%       volumeX <= TolX,
%
%   provided the flag exitflag = 0. 
%
%
%  Examples
%
%  Example 1:
%
%  >> f=@(x) (x-0.3).^2+1; [fmin,out_param] = funminPenalty_g(f)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)(x-0.3).^2+1
%             a: 0
%             b: 1
%        abstol: 1.0000e-06
%          TolX: 1.0000e-03
%           nlo: 10
%           nhi: 1000
%          nmax: 10000000
%         ninit: 100
%           tau: 197
%      exitflag: 0
%       npoints: 6337
%        errest: 6.1554e-07
%       volumeX: 0.0015
%     tauchange: 0
%     intervals: [2x1 double]
%
%
%  Example 2:
%
%  >> f=@(x) (x-0.3).^2+1;
%  >> [fmin,out_param] = funminPenalty_g(f,-2,2,1e-7,1e-4,10,10,1000000)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)(x-0.3).^2+1
%             a: -2
%             b: 2
%        abstol: 1.0000e-07
%          TolX: 1.0000e-04
%           nlo: 10
%           nhi: 10
%          nmax: 1000000
%         ninit: 10
%           tau: 17
%      exitflag: 0
%       npoints: 18433
%        errest: 9.5464e-08
%       volumeX: 5.4175e-04
%     tauchange: 0
%     intervals: [2x1 double]
%
%
%  Example 3:
%
%  >> f=@(x) (x+1.3).^2+1;
%  >> in_param.a = -13; in_param.b = 8;
%  >> in_param.abstol = 10^(-7); in_param.TolX = 1e-4;
%  >> in_param.nlo = 10; in_param.nhi = 100;
%  >> in_param.nmax = 10^6;
%  >> [fmin,out_param] = funminPenalty_g(f,in_param)
% 
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)(x+1.3).^2+1
%             a: -13
%             b: 8
%        abstol: 1.0000e-07
%          TolX: 1.0000e-04
%           nlo: 10
%           nhi: 100
%          nmax: 1000000
%         ninit: 91
%           tau: 179
%      exitflag: 0
%       npoints: 368641
%        errest: 7.1473e-08
%       volumeX: 5.2354e-04
%     tauchange: 0
%     intervals: [2x1 double]
%
%
%  Example 4:
%
%  >> f=@(x) (x-0.3).^2+1;
%  >> [fmin,out_param] = funminPenalty_g(f,'a',-2,'b',2,'nhi',100,'nlo',10,'nmax',1e6,'abstol',1e-4,'TolX',1e-2)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%             f: @(x)(x-0.3).^2+1
%             a: -2
%             b: 2
%        abstol: 1.0000e-04
%          TolX: 0.0100
%           nlo: 10
%           nhi: 100
%          nmax: 1000000
%         ninit: 64
%           tau: 125
%      exitflag: 0
%       npoints: 2017
%        errest: 6.2273e-05
%       volumeX: 0.0146
%     tauchange: 0
%     intervals: [2x1 double]
%
%
%   See also FMINBND, FUNAPPX_G, INTEGRAL_G
%
%  References
%   [1]  Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
%   Univariate Function Minimization," MS thesis, Illinois Institute of 
%   Technology, 2014.
% 
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
%   [MATLAB Software], 2015. Available from http://gailgithub.github.io/GAIL_Dev/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/ 
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice and Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

% Parse and check parameter satisfy conditions or not
[f,out_param] = funminPenalty_g_param(varargin{:});
%in_param = gail.funmin_g_in_param(varargin{:});
%out_param = in_param.toStruct();
%f = in_param.f;
MATLABVERSION = gail.matlab_version;


%% Main algorithm

% initialize number of points
n = out_param.ninit;
% initialize tau
out_param.tau = 2*(n-1)-1;
% cost budget flag
out_param.exitflag = 1;
% tau change flag
tauchange = 0;
% length of interval
len = out_param.b-out_param.a;


while n < out_param.nmax;
    % Stage 1: estimate weaker and stronger norm
    x = out_param.a:len/(n-1):out_param.b;
    y = f(x);
    iSing = find(isinf(y));
    if ~isempty(iSing)
         error('GAIL:funminPenalty_g:yInf',['Function f(x) = Inf at x = ', num2str(x(iSing))]);
    end
    if length(y) == 1  
%     probably f is a constant function and Matlab would  
%     reutrn only a value fmin 
    Un = y;
    interval=[out_param.a out_param.b];
    out_param.exitflag = 0;
    errest = 0;
    volumeX=len;
    
    break;
    end        
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = (n-1)*max(abs(diff_y-(y(n)-y(1))/(n-1)))/len;
    %approximate the stronger norm of input function
    fn = (n-1)^2*max(abs(diff(diff_y)))/len^2;
    
    % Stage 2: satisfy necessary condition of cone
    if out_param.tau*(gn/len+fn/(2*n-2))< fn;
        % increase tau
        out_param.tau = 2*fn/(gn/len+fn/(2*n-2));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n < ((out_param.tau+1)/2);
        % not large enough, increase points number, and go to Stage 1
            l = n;
            n = 1 + (n-1)*ceil(out_param.tau+1/(2*n-2));
        end
    end
        
    % Stage 3: check for convergence
    bn = 2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
    min_index = (2*(n-1)^2.*abs(diff_y)) < (bn*len);
    min_in = min_index.*(diff_y/2+y(1:n-1)-0.25*(2*(n-1)^2*diff_y.^2/bn/len+bn*len/2/(n-1)^2));
    min_end = (~min_index).*(diff_y/2+y(1:n-1)-abs(diff_y)/2);
    ln = min_in+min_end;
    % minimum values of each interval
    Ln = min(ln); % lower bound
    Un = min(y); % upper bound
    errest = Un-Ln;
    % find the intervals containing minimum points
    index = find(min_index ==1 & ln < Un);
    m = size(index,2);        
    if m > 0
        delta = (n-1)^2/len^2*diff_y(index).^2-2*bn/len*(diff_y(index)./2 ...
            +y(index)-bn*len/8/(n-1)^2-Un);
        ints = zeros(2,m);
        ints(1,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
            -len*sqrt(delta)./bn;
        ints(2,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
            +len*sqrt(delta)./bn;
        leftint = find([true diff(index)~=1]);
        rightint = find([diff(index)~=1 true]);
        q = size(leftint,2);
        ints1 = zeros(2,q);
        ints1(1,:) = ints(1,leftint);
        ints1(2,:) = ints(2,rightint);
    else
        ints1 = zeros(2,0);
    end
    k=size(ints1,2);
    % check if [a,a-h] contains minimum
    if abs(y(1)-Un)<out_param.abstol, 
        k=k+1;
        ints2=zeros(2,k);
        ints2(:,1)=[x(1),x(2)];
        if k>1
          ints2(:,2:end)=ints1;
        end
    else
        ints2=ints1;
    end
    % check if [b-h,b] contains minimum
    if abs(y(end)-Un)<out_param.abstol,
        k=k+1;
        ints3=zeros(2,k);
        ints3(:,end)=[x(end-1),x(end)];
        if k>1
          ints3(:,1:end-1)=ints2;
        end
    else
        ints3=ints2;
    end
    interval=ints3;    
    
    volumeX = sum(interval(2,:)-interval(1,:));
    % satisfy convergence
    if errest < out_param.abstol || volumeX < out_param.TolX
        out_param.exitflag = 0; break;
    end
    % otherwise increase points number
    l = n;
    n = 2*(n-1)+1;
end;

% check tau change flag
if tauchange == 1
    warning('GAIL:funminPenalty_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similar functions.')
end;

% check cost budget flag
if out_param.exitflag == 1
    n = l;
    warning('GAIL:funminPenalty_g:exceedbudget','funminPenalty_g attempted to exceed the cost budget. The answer may be unreliable.')
end

fmin = Un;
out_param.npoints = n;
out_param.errest = errest;
out_param.volumeX = volumeX;
out_param.tauchange = tauchange;
out_param.intervals = interval;
out_param = orderfields(out_param, ...
             { 'f','a','b','abstol','TolX','nlo','nhi','nmax','ninit','tau',...
              'exitflag','npoints','errest','volumeX','tauchange','intervals'});

function [f, out_param] = funminPenalty_g_param(varargin)
% Parse the input to the funminPenalty_g.m function

%% Default parameter values

default.a = 0;
default.b = 1;
default.abstol = 1e-6;
default.TolX = 1e-3;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;
    
if isempty(varargin)
    warning('GAIL:funminPenalty_g:nofunction','Function f must be specified. Now funminPenalty_g will use f(x)=(x-0.3)^2+1.')
    help funminPenalty_g
    f = @(x) (x-0.3).^2+1;
    out_param.f = f;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:funminPenalty_g:notfunction','Function f must be a function handle. Now funminPenalty_g will use f(x)=(x-0.3)^2+1.')
    f = @(x) (x-0.3).^2+1;
    out_param.f = f;
  end
end
     
validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) || ischar(in2));
end
        
if ~validvarargin
% There is only one input f or the second input is not satisfied our type.
% The default parameters are used.
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
    out_param.TolX = default.TolX;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2) % more inputs of numerical type. Put them in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'TolX',default.TolX,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) % second input is a structure
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
    end
    parse(p,f,varargin{2:end});
    out_param = p.Results;
end

% a and b can't be infinity
if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('GAIL:funminPenalty_g:anoinfinity',['a cannot be infinity. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('GAIL:funminPenalty_g:bnoinfinity',['b cannot be infinity. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

% b is greater than a
if (out_param.b < out_param.a)
    warning('GAIL:funminPenalty_g:blea','b cannot be smaller than a; exchange these two. ')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('GAIL:funminPenalty_g:beqa',['b cannot equal to a. Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% Check whether the error tolerance is nonnegative
if out_param.abstol < 0
    warning('GAIL:funminPenalty_g:abstolnonpos', ['Error tolerance should be greater than or equal to 0.' ...
        ' Using default error tolerance ', num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Check whether the length tolerance is nonnegative
if out_param.TolX < 0
    warning('GAIL:funminPenalty_g:Xtolnonpos', ['X tolerance should be greater than or equal to 0.' ...
        ' Using default X tolerance ' num2str(default.TolX)]);
    out_param.TolX = default.TolX;
end

% Check whether the cost budget is a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.ispositive(out_param.nmax)
        warning('GAIL:funminPenalty_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))]);
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:funminPenalty_g:budgetnotposint', ['Cost budget should be a positive integer.'...
        ' funminPenalty_g will use the default budget' num2str(out_param.nmax)]);
        out_param.nmax = default.nmax;
    end
end

% Check nlo and nhi
if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('GAIL:funminPenalty_g:lowinitnotint',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('GAIL:funminPenalty_g:lowinitlt3',[' Lower bound of initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('GAIL:funminPenalty_g:hiinitnotint',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('GAIL:funminPenalty_g:hiinitlt3',[' Upper bound of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    warning('GAIL:funminPenalty_g:logrhi', 'Lower bound of initial number of points is larger than upper bound of initial number of points; Use nhi as nlo')
    out_param.nhi = out_param.nlo;
end;
if (out_param.nlo > out_param.nmax)
    warning('GAIL:funminPenalty_g:logecost',['Lower bound of initial number of points should be smaller than cost budget.' ...
            ' Using ', num2str(ceil(out_param.nmax/2))])
    out_param.nlo = ceil(out_param.nmax/2);
    out_param.nhi = out_param.nlo;
end;

out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+out_param.b - out_param.a)));

end

end
