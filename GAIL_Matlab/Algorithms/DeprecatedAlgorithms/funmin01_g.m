function [fmin,out_param]=funmin01_g(varargin)
%FUNMIN01_G 1-D guaranteed global minimum value(s) on [0,1] and the subset 
% containing optimal solutions
%
%   fmin = FUNMIN01_G(f) finds minimum value of function f on the interval
%   [0,1] within a guatanteed absolute error tolerance of 1e-6 and X
%   tolerance of 1e-3. The default initial number of points is 52 and
%   default cost budget is 1e7. Input f is a function handle.
%
%   fmin = FUNMIN01_G(f,abstol,TolX,ninit,nmax) finds minimum value of
%   function f on the interval [0,1] with ordered input parameters:
%   guaranteed absolute error tolerance abstol, guaranteed absolute X
%   tolerance TolX, initial number of points ninit and cost budget nmax.
%
%   fmin = FUNMIN01_G(f,'abstol',abstol,'TolX',TolX,'ninit',ninit,'nmax',
%   nmax) finds minimum value of function f on the interval [0,1] with a
%   guaranteed absolute error tolerance abstol, guaranteed absolute X
%   tolerance TolX, initial number of points ninit and cost budget nmax. All
%   the three field-value pairs are optional and can be supplied in
%   different order.
%
%   fmin = FUNMIN01_G(f,in_param) finds minimum value of function f on the
%   interval [0,1] with a structure input parameters in_param. If a field is
%   not specified, the default value is used.
%
%   [fmin, out_param] = FUNMIN01_G(f,...) returns minimum value fmin of
%   function f and an output structure out_param.
%
%   Input Arguments
%
%     in_param.abstol --- guaranteed absolute error tolerance, default value
%     is 1e-6.
%
%     in_param.abstol --- guaranteed absolute error tolerance, default 
%     value is 1e-6.
%
%     in_param.TolX --- guaranteed X tolerance, default value is 1e-3.
%
%     in_param.ninit --- initial number of points, default value is 52.
%
%     in_param.nmax --- cost budget, default value is 1e7.
%
%   Output Arguments
%
%     out_param.abstol --- guaranteed absolute error tolerance
%
%     out_param.ninit --- initial number of points
%
%     out_param.nmax --- cost budget
%
%     out_param.TolX --- guaranteed X tolerance
%
%     out_param.tau --- latest value of tau
%
%     out_param.exceedbudget --- 0 if the number of points used to find the
%     minimux value is less than the cost budget; 1, otherwise.
%
%     out_param.npoints --- number of points needed to reach the guaranteed
%     absolute error tolerance or the guaranteed X tolerance
%
%     out_param.error --- estimation of the absolute error bound
%
%     out_param.intervals --- the intervals containing point(s) where the
%     minimum occurs
%
%     out_param.volumeX --- the volume of intervals containing the point(s)
%     where the minimum occurs
%
%     out_param.tauchange --- it is 1 if tau is too small, and the algorithm
%     has used a larger tau.
% 
%     out_param.tauchange --- it is 1 if tau is too small, and the 
%     algorithm has used a larger tau.
%
%  Guarantee
%    
%  If the function to be minimized, f, satisfies the cone condition
%      ||f''||_\infty <= \tau ||f'-f(1)+f(0)||_\infty,
%  then the fmin output by this algorithm is guaranteed to satisfy
%      ||min(f)-fmin||_\infty <= abstol
%  or 
%      volumeX <= TolX,
%  provided the flag exceedbudget = 0.
%
%
%  Examples
%
%  Example 1:
%
%  >> f=@(x) (x-0.3).^2+1; [fmin,out_param] = funmin01_g(f)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%           abstol: 1.0000e-06
%             TolX: 1.0000e-03
%            ninit: 52
%             nmax: 10000000
%              tau: 101
%     exceedbudget: 0
%          npoints: 6529
%            error: 2.9618e-07
%          volumeX: 9.8451e-04
%        tauchange: 0
%        intervals: [2x1 double]
%
%
%  Example 2:
%
%  >> f=@(x) (x-0.3).^2+1; in_param.abstol = 1e-8;
%  >> in_param.ninit = 10; in_param.nmax = 1e6;
%  >> in_param.TolX = 1e-4;
%  >> [fmin,out_param] = funmin01_g(f,in_param)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%           abstol: 1.0000e-08
%                f: @(x)(x-0.3).^2+1
%            ninit: 10
%             nmax: 1000000
%             TolX: 1.0000e-04
%              tau: 17
%     exceedbudget: 0
%          npoints: 18433
%            error: 5.9665e-09
%          volumeX: 1.3544e-04
%        tauchange: 0
%        intervals: [2x1 double]
%
%
%  Example 3:
%
%  >> f=@(x) (x-0.3).^2+1;
%  >> [fmin,out_param] = funmin01_g(f,'ninit',10,'nmax',1e6,'abstol',1e-4,'TolX',1e-2)
%
%  fmin =
% 
%     1.0000
% 
%  out_param = 
% 
%           abstol: 1.0000e-04
%                f: @(x)(x-0.3).^2+1
%            ninit: 10
%             nmax: 1000000
%             TolX: 0.0100
%              tau: 17
%     exceedbudget: 0
%          npoints: 145
%            error: 9.4167e-05
%          volumeX: 0.0170
%        tauchange: 0
%        intervals: [2x1 double]
% 
%
%   See also FUNAPPX_G, INTEGRAL_G
%
%  References
%   [1]  Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, 
%   and Yizhi Zhang. The Cost of Deterministic, Adaptive, Automatic 
%   Algorithms: Cones, Not Balls. Journal of Complexity, 30:21-45, 2014
%
%   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   and Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library
%   (Version 1.3.0)" [MATLAB Software], 2014. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   [3]  Xin Tong. A Guaranteed, Adaptive, Automatic Algorithm for
%   Univariate Function Minimization. 2014
%
%   If you find GAIL helpful in your work, please support us by citing
%   the above paper and software.
%


% Parse and check the validity of input parameters
[f,out_param] = funmin01_g_param(varargin{:});


%% Main algorithm

% initialize number of points
n = out_param.ninit;
% initialize tau
out_param.tau = ceil((n-1)*2-1);
% cost budget flag
out_param.exceedbudget = 1;
% tau change flag
tauchange = 0;

while n < out_param.nmax;
    % Stage 1: estimate weaker and stronger norm
    x = (0:n-1)/(n-1);
    y = f(x);
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = (n-1)*max(abs(diff_y-(y(n)-y(1))/(n-1)));
    %approximate the stronger norm of input function
    fn = (n-1)^2*max(abs(diff(diff_y)));
    
    % Stage 2: satisfy necessary condition of conev
    if out_param.tau*(gn+fn/(2*n-2))>= fn;
        % Stage 3: check for convergence
        bn = 2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
        min_index = 2*(n-1)^2.*abs(diff_y) < bn;
        min_in = min_index.*(diff_y/2+y(1:n-1)-0.25*(2*(n-1)^2*diff_y.^2/bn+bn/2/(n-1)^2));
        min_end = (~min_index).*(diff_y/2+y(1:n-1)-abs(diff_y)/2);
        ln = min_in+min_end;
        % minimum values of each interval
        Ln = min(ln); % lower bound
        Un = min(y); % upper bound
        error = Un-Ln;
        % find the intervals containing minimum points
        index = find(min_index ==1 & ln < Un);
        m = size(index,2);
        if m > 0
            delta = (n-1)^2*diff_y(index).^2-2*bn*(diff_y(index)./2 ...
                +y(index)-bn/8/(n-1)^2-Un);
            ints = zeros(2,m);
            ints(1,:)=x(index)+1/2/(n-1)-(n-1)*diff_y(index)./bn ...
                -sqrt(delta)./bn;
            ints(2,:)=x(index)+1/2/(n-1)-(n-1)*diff_y(index)./bn ...
                +sqrt(delta)./bn;
            leftint = find([1 diff(index)~=1]);
            rightint = find([diff(index)~=1 1]);
            q = size(leftint,2);
            interval = zeros(2,q);
            interval(1,:) = ints(1,leftint);
            interval(2,:) = ints(2,rightint);
        else
            interval = zeros(2,0);
        end
        volumeX = sum(interval(2,:)-interval(1,:));
        % satisfy convergence
        if error < out_param.abstol || volumeX < out_param.TolX
            out_param.exceedbudget = 0; break;
        end
        % otherwise increase points number
        l = n;
        n = 2*(n-1)+1;
        
        % Stage 2: do not satisfy necessary condition
    else
        % increase tau
        out_param.tau = 2*fn/(gn+fn/(2*n-2));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= ((out_param.tau+1)/2);
            % large enough, go to Stage 3
            bn = 2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
            min_index = 2*(n-1)^2.*abs(diff_y) < bn;
            min_in = min_index.*(diff_y/2+y(1:n-1)-0.25*(2*(n-1)^2*diff_y.^2/bn+bn/2/(n-1)^2));
            min_end = (~min_index).*(diff_y/2+y(1:n-1)-abs(diff_y)/2);
            ln = min_in+min_end;
            % minimum values of each interval
            Ln = min(ln); % lower bound
            Un = min(y); % upper bound
            error = Un-Ln;
            % find the intervals containing minimum points
            index = find(min_index ==1 & ln < Un);
            m = size(index,2);
            if m > 0
                delta = (n-1)^2*diff_y(index).^2-2*bn*(diff_y(index)./2 ...
                    +y(index)-bn/8/(n-1)^2-Un);
                ints = zeros(2,m);
                ints(1,:)=x(index)+1/2/(n-1)-(n-1)*diff_y(index)./bn ...
                    -sqrt(delta)./bn;
                ints(2,:)=x(index)+1/2/(n-1)-(n-1)*diff_y(index)./bn ...
                    +sqrt(delta)./bn;
                leftint = find([1 diff(index)~=1]);
                rightint = find([diff(index)~=1 1]);
                q = size(leftint,2);
                interval = zeros(2,q);
                interval(1,:) = ints(1,leftint);
                interval(2,:) = ints(2,rightint);
            else
                interval = zeros(2,0);
            end
            volumeX = sum(interval(2,:)-interval(1,:));
            % satisfy convergence
            if error < out_param.abstol || volumeX < out_param.TolX
                out_param.exceedbudget = 0; break;
            end
            % otherwise increase points number
            l = n;
            n = 2*(n-1)+1;
        else
            % not large enough, increase points number, and go to Stage 1
            l = n;
            n = 1 + (n-1)*ceil(out_param.tau+1/(2*n-2));
        end;
    end;
end;


% check tau change flag
if tauchange == 1
    warning('GAIL:funmin01_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similar functions.')
end;

% check cost budget flag
if out_param.exceedbudget == 1
    n = l;
    warning('GAIL:funmin01_g:exceedbudget','funmin01_g attempted to exceed the cost budget. The answer may be unreliable.')
end

fmin = Un;
out_param.npoints = n;
out_param.error = error;
out_param.volumeX = volumeX;
out_param.tauchange = tauchange;
out_param.intervals = interval;


function [f, out_param] = funmin01_g_param(varargin)
% Parse the input to the funmin01_g.m function
default.abstol = 1e-6;
default.TolX = 1e-3;
default.ninit = 52;
default.nmax = 1e7;
    
if isempty(varargin)
    help funmin01_g
    warning('Function f must be specified. Now funmin01_g will use f(x)=(x-0.3)^2+1.')
    f = @(x) (x-0.3).^2+1;
else
    f = varargin{1};
end
     
validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) || ischar(in2));
end
        
if ~validvarargin
% There is only one input f or the second input is not satisfied our type.
% The default parameters are used.
    out_param.abstol = default.abstol;
    out_param.TolX = default.TolX;
    out_param.ninit = default.ninit;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2) % more inputs of numerical type. Put them in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'TolX',default.TolX,@isnumeric);
        addOptional(p,'ninit',default.ninit,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) % second input is a structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'TolX',default.TolX,@isnumeric);
        addParamValue(p,'ninit',default.ninit,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end});
    out_param = p.Results;
end
        
% Check whether the error tolerance is positive
if out_param.abstol < 0
    warning(['Error tolerance should be greater than or equal to 0.'...
        ' funmin01_g will use the default error tolerance' ... 
        num2str(default.abstol)]);
    out_param.abstol = default.abstol;
end

% Check whether the length tolerance is positive
if out_param.TolX < 0
    warning(['Tolerance on X should be greater than or equal to 0.'...
    ' funmin01_g will use the default X tolerance ' num2str(default.TolX)]);
    out_param.abstol = default.TolX;
end
        
% Check whether the initial number of points is a positive integer
if (~gail.isposint(out_param.ninit))
    if gail.ispositive(out_param.ninit)
        warning(['Initial number of points should be a integer.' ...
            ' funmin01_g will use ' num2str(ceil(out_param.ninit))]);
        out_param.ninit = ceil(out_param.ninit);
    else
        warning(['Initial number of points should be a positive integer.' ...
            ' funmin01_g will use the default initial number of points ' ...
            num2str(out_param.ninit)]);
        out_param.ninit = default.ninit;
    end
end
        
% Check whether the cost budget is a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.ispositive(out_param.nmax)
        warning(['Cost budget should be a integer.'' funmin01_g will use ' ...
            num2str(ceil(out_param.nmax))]);
        out_param.nmax = ceil(out_param.nmax);
    else
        warning(['Cost budget should be a positive integer.'...
        ' funmin01_g will use the default budget' num2str(out_param.nmax)]);
        out_param.nmax = default.nmax;
    end
end

end

end
