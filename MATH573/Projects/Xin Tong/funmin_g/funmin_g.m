function [fmin,out_param]=funmin_g(varargin)
%FUNMIN_G Guaranteed minimum value of one-dimensional function
%on interval [0,1] and the subset containing point(s) where the minimum occurs 
%
%  fmin = FUNMIN_G(f) finds minimum value of function f on the interval
%  [0,1] within a guaranteed absolute error of 1e-6. The default initial
%  number of points is 52 and defoult cost budget is 1e7. 
%
%  fmin = FUNMIN_G(f,abstol,ninit,nmax) finds minimum value of function f
%  with ordered input parameters: guaranteed absolute error abstol, initial
%  number of points ninit and cost budget nmax.
%
%  fmin = FUNMIN_G(f,'abstol',abstol,'ninit',ninit,'nmax',nmax) finds
%  minimum value of function f with the guaranteed absolute error abstol,
%  initian number of points ninit and cost budget nmax. All the three
%  field-value pairs are optional and can be supplied in different order.
% 
%  fmin = FUNMIN_G(f,in_param) finds minimum value of function f with a
%  structure input parameters in_param. If a field is not specified, the
%  default value is used.
%  
%    in_param.abstol --- guaranteed absolute error, defoult value is 1e-6.
%  
%    in_param.nint --- initial number of points, default value is 52.
%  
%    in_param.nmax --- cost budget, default value is 1e7.
%  
%  [fmin, out_param] = FUNMIN_G(f,...) returns minimum value fmin of
%  function and an output structure out_param, which has the following
%  fields.
%
%    out_param.abstol --- guaranteed absolute error
% 
%    out_param.ninit --- initial number of points 
% 
%    out_param.nmax --- cost budget
% 
%    out_param.tau --- latest value of tau
% 
%    out_param.exceedbudget --- 0 if the number of points used to find the
%    minimux value is less than the cost budget; 1, otherwise.
% 
%    out_param.npoints --- number of points needed to reach the guaranteed
%    absolute error
%    
%    out_param.error --- estimation of the absolute error bound
%
%    out_param.intervals --- the points of minimum values are located in
%    those intervals
% 
% 
%  Examples  
% 
%  Example 1:
% 
%  >> f=@(x) (x-0.3).^2+1; [fmin,out_param] = funmin_g(f)
%  fmin =
% 
%    0.999999855664237
% 
%  out_param = 
% 
%           abstol: 1.000000000000000e-06
%            ninit: 52
%             nmax: 10000000
%              tau: 101
%     exceedbudget: 0
%          npoints: 6529
%            error: 2.961806568890779e-07
%        intervals: [2x7 double]
%  
%  >> out_param.intervals
% 
%  ans =
% 
%   Columns 1 through 2
% 
%    0.299519512716009   0.299649912671689
%    0.299609467134614   0.299779418954609
% 
%   Columns 3 through 4
% 
%    0.299790957911706   0.299938725490181
%    0.299938725490211   0.300091309687384
% 
%   Columns 5 through 6
% 
%    0.300092543773805   0.300253471873882
%    0.300237843179407   0.300377266854949
% 
%   Column 7
% 
%    0.300427068709191
%    0.300504021795287
% 
% 
%  Example 2:
%
%  >> f=@(x) (x-0.3).^2+1; in_param.abstol = 1e-8;
%  >> in_param.ninit = 10; in_param.nmax = 1e6;
%  >> [fmin,out_param] = funmin_g(f,in_param)
%  
%  fmin =
% 
%    0.999999997487714
% 
%  out_param = 
% 
%           abstol: 1.000000000000000e-08
%                f: @(x)(x-0.3).^2+1
%            ninit: 10
%             nmax: 1000000
%              tau: 17
%     exceedbudget: 0
%          npoints: 18433
%            error: 5.966472205187756e-09
%        intervals: [2x3 double]
%  
%  >> out_param.intervals
%  
%  ans =
% 
%   Columns 1 through 2
% 
%    0.299929032946815   0.299968723949380
%    0.299965645775326   0.300021701388822
% 
%   Column 3
% 
%    0.300021701388667
%    0.300064470566076
%
% 
%  Example 3:
%  
%   >> f=@(x) (x-0.3).^2+1; 
%   >> [fmin,out_param] = funmin_g(f,'ninit',10,'nmax',1e6,'abstol',1e-8)
% 
%   fmin =
% 
%    0.999999997487714
% 
%   out_param = 
% 
%           abstol: 1.000000000000000e-08
%                f: @(x)(x-0.3).^2+1
%            ninit: 10
%             nmax: 1000000
%              tau: 17
%     exceedbudget: 0
%          npoints: 18433
%            error: 5.966472205187756e-09
%        intervals: [2x3 double]
% 
%  >> out_param.intervals
% 
%  ans =
% 
%   Columns 1 through 2
% 
%    0.299929032946815   0.299968723949380
%    0.299965645775326   0.300021701388822
% 
%   Column 3
% 
%    0.300021701388667
%    0.300064470566076
% 
% 
% 
%  Sea also FUNAPPX_G
% 
%  Reference
%  [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang,
%       The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones,
%       Not Balls, Journal of Complexity (2013), to appear, DOI
%       10.1016/j.jco.2013.09.002
%



% Parse and check the validity of input parameters
[f,out_param] = funmin_g_param(varargin{:}); 



%% Main algorithm

% initialize number of points
n = out_param.ninit;
% initialize tau
out_param.tau = ceil((n-1)*2-1);
% cost bugdet flag
out_param.exceedbudget = 1;
% tau change flag
tauchange = 0;


while n < out_param.nmax;
    % Stage 1: estimate weaker and stronger norm
    x = (0:n-1)/(n-1);
    y = f(x);
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = max((n-1)*abs(diff_y-(y(n)-y(1))/(n-1)));
    %approximate the stronger norm of input function
    fn = max((n-1)^2*abs(diff(diff_y)));
    
    % Stage 2: satisfy necessary condition of cone
    if out_param.tau*(gn+fn/(2*n-2))>= fn;
        % Stage 3: check for convergence 
        bn=2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
        cn=2*(n-1)^2*abs(diff_y)./bn; 
        Cn=min(cn,1); % check the conditions for each interval
        ln=(diff_y/2+y(1:n-1))  -abs(diff_y).*(Cn+1./Cn)/4; % minimum values of each interval
        Ln=min(ln); % lower bound
        min_endpoint=min(y); % upper bound
        error=min_endpoint-Ln;
        % satisfy convergence
        if error < out_param.abstol
            out_param.exceedbudget = 0; break;
        end
        % otherwise increase points number
        n=2*(n-1)+1;
    % Stage 2: do not satisfy necessary condition
    else
        % increase tau
        out_param.tau = 2*fn/(gn+fn/(2*n-2));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= ((out_param.tau+1)/2);
            % large enough, go to Stage 3
            bn=2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
            cn=2*(n-1)^2*abs(diff_y)./bn; 
            Cn=min(cn,1); % check the conditions for each interval
            ln=(diff_y/2+y(1:n-1))  -abs(diff_y).*(Cn+1./Cn)/4; % minimum values of each interval
            Ln=min(ln); % lower bound
            min_endpoint=min(y); % upper bound
            error=min_endpoint-Ln;
            % satisfy convergence
            if error < out_param.abstol
                out_param.exceedbudget = 0; break;
            end
            % otherwise increase points number
            n=2*(n-1)+1;
        
        else
        % not large enough, increase points number, and go to Stage 1
            n = 1 + (n-1)*ceil(out_param.tau+1/(2*n-2));
        end;
    end;
end;

% find the intervals containing minimum points 
index=find(cn<1 & ln <= min_endpoint);
m=size(index,2);
ints=zeros(2,m);
delta=zeros(m);
for i=1:m
    delta(i)=(n-1)^2*diff_y(index(i))^2-2*bn*(diff_y(index(i))/2+y(index(i))-bn/8/(n-1)^2-min_endpoint);
    ints(:,i)=x(index(i))+1/2/(n-1)-(n-1)*diff_y(index(i))/bn+[-1,1]*sqrt(delta(i))/bn;
end
   
    
% check tau change flag
if tauchange == 1
    warning('MATLAB:funmin_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similiar functions.')
end;

% check cost budget flag
if out_param.exceedbudget == 1
    warning('MATLAB:funmin_g:exceedbudget','funmin_g attemped to exceed the cost budget. The answer may be unreliable.')
end

fmin = min_endpoint;
out_param.npoints = n;
out_param.error = error;
out_param.intervals = ints;





function [f, out_param] = funmin_g_param(varargin)
% Parse the input to the funmin_g.m function

default.abstol = 1e-6;
default.ninit = 52;
default.nmax = 1e7;

if isempty(varargin)
    help funmin_g
    warning('Function f must be specified. Now funmin_g will use f(x)=(x-0.3)^2+1.')
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
    % There is only one input f or the second input is not satisfied our
    % type. The default parameters are used.
    out_param.abstol = default.abstol;
    out_param.ninit = default.ninit;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in2) %There are more inputs of numerical type. Put them in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'ninit',default.ninit,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) % second input is a structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'ninit',default.ninit,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end});
    out_param = p.Results;
end

% Check whether the tolerance is positive 
if out_param.abstol<=0
    warning(['Error tolerance should be greater that 0.'...
        ' funmin_g will use the default tolerance' num2str(default.abstol)]);
    out_param.abstol = default.abstol;
end

% Check whether the initial number of points is a positive integer
if (~isposint(out_param.ninit))
    if ispositive(out_param.ninit)
        warning(['Initial number of points should be a integer.' ...
        ' funmin_g will use ' num2str(ceil(out_param.ninit))]);
        out_param.ninit = ceil(out_param.ninit);
    else
        warning(['Initial number of points should be a positive integer.' ...
        ' funmin_g will use the default initial number of points ' ...
        num2str(out_param.ninit)]);
        out_param.ninit = default.ninit;
    end
end

% Check whether the cost budget is a positive integer
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning(['Cost budget should be a integer.'' funmin_g will use ' ...
        num2str(ceil(out_param.nmax))]);
        out_param.nmax = ceil(out_param.nmax);
    else
        warning(['Cost budget should be a positive integer.'...
        ' funmin_g will use the default cost budget' num2str(out_param.nmax)]);
        out_param.nmax = default.nmax;
    end
end
        

end

end
