function [fappx,out_param]=funappxtau_g(f,varargin)
%FUNAPPXTAU_G One dimensional guaranteed function recovery on interval 
%   [0,1] with cone condition tau
%
%   fappx = FUNAPPXTAU_G(f) recovers function f on the interval [0,1] by a 
%   piecewise linear interpolant fappx to within a guaranteed absolute 
%   error of 1e-6. Default cone constant is 10 and default cost budget is
%   1e7.  Input f a function handle. The statement Y=f(X) should accept a 
%   a vector argument X and return a vector Y of function values that is 
%   the same size as X.
%
%   fappx = FUNAPPXTAU_G(f,abstol,tau,nmax) for given function f and the 
%   ordered input parameters with the guaranteed absolute error abstol, 
%   cone condition tau, and cost budget nmax.
%
%   fappx = FUNAPPXTAU_G(f,'abstol',abstol,'tau',tau,'nmax',nmax) recovers
%   function f with the guaranteed absolute error abstol, cone constant tau,
%   and cost budget nmax. All three filed-value pairs are optional and can 
%   be supplied in different order.
%
%   fappx = FUNAPPXTAU_G(f,in_param) recovers function f with the guaranteed
%   absolute error in_param.abstol, cone constant in_param.tau, and cost 
%   budget in_param.nmax. If a field is not specified, the default value is
%   used.
%
%   in_param.abstol --- guaranteed absolute error, default value is 1e-6.
%
%   in_param.tau --- cone constant, default value is 10.
%
%   in_param.nmax --- cost budget, default value is 1e7.
%
%   [fappx, out_param] = FUNAPPXTAU_G(f,...) returns an approximated function
%   fappx and an output structure out_param, which has the following fields.
%
%   out_param.nmax --- cost budget, default value is 1e7.
%
%   out_param.exceedbudget --- it is 0 if the number of points used in the 
%   construction of fappx is less than cost budget, 1 otherwise.
%
%   out_param.npoints --- number of points we need to reach the guaranteed
%   absolute error tol
%
%   out_param.tauchange --- it is 1 if tau is too small, and the algorithm 
%   has used a larger tau.
%
%   out_param.tau --- if input tau is changed, this field will give new 
%   value of tau
%
%   out_param.abtol --- guaranteed absolute error, default value is 1e-6.
%
%   
%   Examples
%
%   Example 1:
%   
%
%   >> format short; f = @(x) x.^2; [fappx, out_param] = funappxtau_g(f)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***6
%              tau: 10
%             nmax: 10000000
%     exceedbudget: 0
%        tauchange: 0
%       ballradius: 2.0000e+***7
%          npoints: 2053
% 
%
%   Example 2:
%
%   >> clear in_param; format short; in_param.abstol = 10^(-8); 
%   >> in_param.tau = 15; in_param.nmax = 10^6; 
%   >> [fappx, out_param] = funappxtau_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%             nmax: 1000000
%              tau: 15
%     exceedbudget: 0
%        tauchange: 0
%       ballradius: 1.3333e+***3
%          npoints: 25633
%
%
%   Example 3:
%
%   >> format short; f = @(x) x.^2; 
%   >> [fappx, out_param] = funappxtau_g(f,'tau',15,'nmax',1e6,'abstol',1e-8)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%             nmax: 1000000
%              tau: 15
%     exceedbudget: 0
%        tauchange: 0
%       ballradius: 1.3333e+***3
%          npoints: 25633
%
%
%
%   See also INTEGRAL_G, MEANMC_G, CUBMC_G
%
%   Reference
%   [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%        The Cost of Deterministic, Adaptive, Automatic Algorithms:  Cones, 
%        Not Balls, submitted for publication, arXiv.org:1303.2412 
%        [math.NA]}, 2013.
%

default.abstol  = 1e-6;
default.tau  = 10;
default.nmax  = 1e7;

if (nargin<1)
    help funappxtau_g
    warning('Function f must be specified. Now GAIL is using f(x)=x^2.')
    f = @(x) x.^2;
end;

if (nargin<2)
   out_param.abstol = default.abstol;
   out_param.tau = default.tau;
   out_param.nmax = default.nmax;
end;

p = inputParser;
addRequired(p,'f',@isfcn);

%% API format: (fcn, struct)
if (nargin == 2 && isstruct(varargin{1}))
 p.StructExpand = true;
 addParamValue(p,'abstol',default.abstol,@isnumeric);
 addParamValue(p,'tau',default.tau,@isnumeric);
 addParamValue(p,'nmax',default.nmax,@isnumeric);
 parse(p,f,varargin{:})
 out_param = p.Results;
end

%% API format---not in order: (fcn, 'input2', inputVal2, 'input3', inputVal3, 'input1', inputVal1)
if (nargin > 2)
  in2 = varargin{1};
  if (ischar(in2)),
    addParamValue(p,'abstol',default.abstol,@isnumeric);
    addParamValue(p,'tau',default.tau,@isnumeric);
    addParamValue(p,'nmax',default.nmax,@isnumeric);
    parse(p,f,varargin{:})
    out_param = p.Results;
  end
end

%% API format---in order: (fcn, inputVal1, inputVal2, inputVal3) 
if (nargin >= 2 && isnumeric(varargin{1}))
    addOptional(p,'abstol',default.abstol,@isnumeric);
    addOptional(p,'tau',default.tau,@isnumeric);
    addOptional(p,'nmax',default.nmax,@isnumeric);
    parse(p,f,varargin{:})
    out_param = p.Results;
end

% check parameter satisfy conditions or not
out_param = funappx_g_param(out_param,default);


%% main algorithm

% initialize number of points
n = ceil((out_param.tau+1)/2)+1;
% cost budget flag
out_param.exceedbudget = 1;
% tau change flag
out_param.tauchange = 0;

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
        out_param.tauchange = 1;
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
    warning('MATLAB:funappxtau_g:exceedbudget','funappx_g attemped to exceed the cost bugdet. The answer may be unreliable.')
end;

if out_param.tauchange == 1;
    warning('MATLAB:funappxtau_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similiar functions.')
end;

out_param.ballradius = 2*out_param.abstol*(out_param.nmax-2)*(out_param.nmax...
    -2-out_param.tau)/out_param.tau;
out_param.npoints = n;
x1 = 0:1/(out_param.npoints-1):1;
y1 = f(x1);
fappx = @(x) interp1(x1,y1,x,'linear');

function out_param = funappx_g_param(out_param,default)
% let error tolerance less than 1 and greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let cone condition greater or equal 2
if (out_param.tau < 2)
    warning(['Cone condition should be greater or equal 2.' ...
             ' Using default cone condition ' num2str(default.tau)])
    out_param.tau = default.tau;
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
