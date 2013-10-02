function [fappx,out_param]=funappx_g(varargin)
%FUNAPPX_G One-dimensional guaranteed function recovery on interval [0,1]
%
%   fappx = FUNAPPX_G(f) recovers function f on the interval [0,1] by a 
%   piecewise linear interpolant fappx to within a guaranteed absolute 
%   error of 1e-6. Default initial number of points is 52 and default cost
%   budget is 1e7.  Input f is a function handle. The statement Y=f(X)
%   should accept a vector argument X and return a vector Y of function
%   values that is the same size as X.
%
%   fappx = FUNAPPX_G(f,abstol,ninit,nmax) for given function f and the ordered 
%   input parameters with the guaranteed absolute error abstol, initial number
%   of points ninit and cost budget nmax.
%
%   fappx = FUNAPPX_G(f,'abstol',abstol,'ninit',ninit,'nmax',nmax) recovers 
%   function f with the guaranteed absolute error abstol, initial number of 
%   points ninit, and cost budget nmax. All three field-value pairs are
%   optional and can be supplied in different order.
%
%   fappx = FUNAPPX_G(f,in_param) recovers function f with the guaranteed
%   absolute error in_param.abstol, initial number of points in_param.ninit,
%   and cost budget in_param.nmax. If a field is not specified, the default
%   value is used.
%
%   in_param.abstol --- guaranteed absolute error, default value is 1e-6.
%
%   in_param.ninit --- initial number of points we used, default value is
%   52
%
%   in_param.nmax --- cost budget, default value is 1e7
%
%   [fappx, out_param] = FUNAPPX_G(f,...) returns an approximated function
%   fappx and an output structure out_param, which has the following fields.
%
%
%   out_param.nmax --- cost budget
%
%   out_param.exceedbudget --- it is 0 if the number of points used in the 
%   construction of fappx is less than cost budget, 1 otherwise.
%
%   out_param.ninit --- initial number of points we used
%
%   out_param.npoints --- number of points we need to reach the guaranteed
%   absolute error
%
%   out_param.errorbound --- estimation of the approximation absolute error
%   bound
%
%   out_param.tau --- this field will give latest value of tau
%
%   out_param.abstol --- guaranteed absolute error
%
%   
%   Examples
%
%   Example 1:
%   
%
%   >> f = @(x) x.^2; [fappx, out_param] = funappx_g(f)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%            ninit: 52
%             nmax: 10000000
%              tau: 101
%     exceedbudget: 0
%          npoints: 7039
%         errbound: 5.0471e-***9
%          
% 
%
%   Example 2:
%
%   >> clear in_param; in_param.abstol = 10^(-8); 
%   >> in_param.ninit = 10; in_param.nmax = 10^6; 
%   >> [fappx, out_param] = funappx_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 1000000
%              tau: 17
%     exceedbudget: 0
%          npoints: 27505
%         errbound: 3.3048e-***10
%
%
%   Example 3:
%
%   >> f = @(x) x.^2; 
%   >> [fappx, out_param] = funappx_g(f,'ninit',10,'nmax',1e6,'abstol',1e-8)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 1000000
%              tau: 17
%     exceedbudget: 0
%          npoints: 27505
%         errbound: 3.3048e-***10
%
%
%   See also INTEGRAL_G, MEANMC_G
%
%   Reference
%   [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%        The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, 
%        Not Balls, submitted for publication, arXiv.org:1303.2412 
%        [math.NA]}, 2013.
%

% check parameter satisfy conditions or not
[f, out_param] = funappx_g_param(varargin{:});

%% main algorithm

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
    gn = max((n-1)*abs(diff_y-(y(n)-y(1))/(n-1)));
    %approximate the stronger norm of input function
    fn = max((n-1)^2*abs(diff(diff_y)));
    
    % Stage 2: satisfy necessary condition
    if out_param.tau*(gn+fn/(2*n-2))>= fn;
        % Stage 3: check for convergence
        errbound = 4*out_param.abstol*(n-1)*(2*n-2-out_param.tau)/out_param.tau;
        % satisfy convergence
        if errbound >= gn;
            out_param.exceedbudget = 0; break;
        end;
        % otherwise increase number of points
        m = max(ceil(1/(n-1)*sqrt(gn*out_param.tau/8/out_param.abstol)),2);
        n = m*(n-1)+1;
    % Stage2: do not satisfy necessary condition
    else
        % increase tau
        out_param.tau = 2*fn/(gn+fn/(2*n-2));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= ((out_param.tau+1)/2);
            % true, go to Stage 3
            errbound = 4*out_param.abstol*(n-1)*(2*n-2-out_param.tau)...
                /out_param.tau;
            if errbound >= gn;
                out_param.exceedbudget = 0; break;
            end;
            m = max(ceil(1/(n-1)*...
                sqrt(gn*out_param.tau/8/out_param.abstol)),2);
            n = m*(n-1)+1;
        else
            % otherwise increase number of points, go to Stage 1
            n = 1 + (n-1)*ceil(out_param.tau+1/(2*n-2));
        end;
    end;
end;

% Check cost budget flag
if out_param.exceedbudget == 1;
    n = 1 + (n-1)/m*floor((out_param.nmax-1)*m/(n-1));
    warning('MATLAB:funappx_g:exceedbudget','funappx_g attemped to exceed the cost bugdet. The answer may be unreliable.')
end;

if tauchange == 1;
    warning('MATLAB:funappx_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similiar functions.')
end;
%out_param.ballradius = 2*out_param.abstol*(out_param.nmax-2)*(out_param.nmax...
%    -2-out_param.tau)/out_param.tau;
out_param.npoints = n;
out_param.errbound = fn/(8*(n-1)^2);
x1 = 0:1/(out_param.npoints-1):1;
y1 = f(x1);
fappx = @(x) interp1(x1,y1,x,'linear');

function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol  = 1e-6;
default.ninit  = 52;
default.nmax  = 1e7;


if isempty(varargin)
    help funappx_g
    warning('Function f must be specified. Now GAIL is using f(x)=x^2.')
    f = @(x) x.^2;
else
    f = varargin{1};
end;

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
    %if only one input f, use all the default parameters
    out_param.abstol = default.abstol;
    out_param.ninit = default.ninit;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'ninit',default.ninit,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'ninit',default.ninit,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let initial number of points be a positive integer
if (~isposint(out_param.ninit) && isposge3(out_param.ninit))
    warning(['Initial number of points should be a positive integer.' ...
             ' Using ', num2str(ceil(out_param.ninit))])
    out_param.ninit = ceil(out_param.ninit);
elseif(~isposint(out_param.ninit) && ~isposge3(out_param.ninit))
    warning(['Initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.ninit)])
    out_param.ninit = default.ninit;
end
% let cost budget be a positive integer
if (~isposint(out_param.nmax) && ispositive(out_param.nmax))
    warning(['Cost budget should be a positive integer.' ...
             ' Using cost budget ', num2str(ceil(out_param.nmax))])
    out_param.nmax = ceil(out_param.nmax);
elseif(~isposint(out_param.nmax) && ~ispositive(out_param.nmax))
    warning(['Cost budget should be a positive integer.' ...
             ' Using default cost budget ' int2str(default.nmax)])
    out_param.nmax = default.nmax;
end
