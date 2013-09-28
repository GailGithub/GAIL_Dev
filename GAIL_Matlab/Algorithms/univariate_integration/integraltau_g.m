function [Q,out_param]= integraltau_g(f,varargin)
%  INTEGRALTAU_G 1-D guaranteed function integration using trapezoidal rule
% 
%  Description
%
%   q = INTEGRALTAU_G(f) computes q, the definite integral of function f
%   on the interval [0,1] by trapezoidal rule with 
%   in a guaranteed absolute error of 1e-6. Default cone constant is 100
%   and default cost budget is 1e7. Input f is a 
%   function handle. The function y = f(x) should accept a vector argument 
%   x and return a vector result y, the integrand evaluated at each element
%   of x.
%
%   q = INTEGRALTAU_G(f,in_param) computes q, the definite integral of 
%   function f by trapezoidal rule within a guaranteed absolute error
%   in_param.abstol, cone constant in_param.tau, and cost 
%   budget in_param.nmax. If a field is not specified, the default value is
%   used.
%   
%   in_param.abstol --- absolute error tolerance required by user.
% 
%   in_param.tau --- cone constant specified by user.
%
%   in_param.nmax --- cost budget required by user.
% 
%   q = INTEGRALTAU_G(f,'abstol',abstol,'tau',tau,'nmax',nmax) computes
%   q, the definite integral of function f by trapezoidal rule within a 
%   guaranteed absolute error tolerance abstol, cone constant tau, 
%   and cost budget nmax. All three field-value pairs are optional 
%   and can be supplied.
%
%   q = INTEGRALTAU_G(f,abstol,tau, nmax) computes q, the defintie 
%   integral of function f by trapezoidal rule with the ordered input 
%   parameters, guaranteed absolute error tolerance abstol, cone constant
%   tau, and cost budget nmax.
%
%   [q, out_param] = INTEGRALTAU_G(f,...) returns the approximated 
%   integration q and output structure out_param, which includes the 
%   fileds in_param plus the following fields:
%
%   out_param.exceedbudget --- it is true if the algorithm tries to use 
%   more points than cost budget, false otherwise.
% 
%   out_param.tauchange --- it is true if the cone constant has been
%   changed, false otherwise. See [1] for details. If true, you may wish to
%   change the input in_param.tau to a larger number.
% 
%   out_param.npoints --- number of points we need to 
%   reach the guaranteed absolute error tolernce abstol.
%
%   out_param.errest --- approximation error defined as the differences
%   between the true value and the approximated value of the integral.
%
% 
%   Examples
%
%   Example 1: 
%   >> f = @(x) x; fint = integraltau_g(f)
%   fint = 5.0000e-01
%
%
%   Example 2:
%   >> f = @(x) x.^2; fint = integraltau_g(f,'tol',0.00001,'tau',100,'nmax',1e7)
%   fint = 3.3333e-01
%
%
%   Example 3:
%   >> fint = integraltau_g()
%   Warning: Function f must be specified. Now GAIL is giving you a toy example of
%   f(x)=x^2. >  In ***
%   fint = 0.3333
%
% 
% See Also
% funappx_g, meanMC_g
%
% Reference:
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%      The complexity of guaranteed automatic algorithms: Cones, not
%      balls, preprint, 2013, arXiv:1303.2412 [math.ST].

%%
default.tol  = 1e-6;    % default tolerance
default.tau  = 100;  % default cone constant
default.nmax  = 1e7;    % default maximun allowance of number of points

if (nargin<1)
    help integraltau_g
    warning('Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
end;

if (nargin<2)
   out_param.tol = default.tol;
   out_param.tau = default.tau;
   out_param.nmax = default.nmax;
end;

p = inputParser;
addRequired(p,'f',@isfcn);

%% API format: (fcn, struct)
if (nargin == 2 && isstruct(varargin{1}))
    p.StructExpand = true;
    addParamValue(p,'tol',default.tol,@ispositive);
    addParamValue(p,'tau',default.tau,@ispositive);
    addParamValue(p,'nmax',default.nmax,@ispositive);
    parse(p,f,varargin{:})
    out_param = p.Results;
end

%% API format---not in order: (fcn, 'input2', inputVal2, 'input3', inputVal3, 'input1', inputVal1)
if (nargin > 2)
    in2 = varargin{1};
    if (ischar(in2)),
        addParamValue(p,'tol',default.tol,@ispositive);
        addParamValue(p,'tau',default.tau,@ispositive);
        addParamValue(p,'nmax',default.nmax,@ispositive);
        parse(p,f,varargin{:})
        out_param = p.Results;
    end
end

%% API format---in order: (fcn, inputVal1, inputVal2, inputVal3) 
if (nargin >= 2 && isnumeric(varargin{1}))
    addOptional(p,'tol',default.tol,@ispositive);
    addOptional(p,'tau',default.tau,@ispositive);
    addOptional(p,'nmax',default.nmax,@ispositive);
    parse(p,f,varargin{:})
    out_param = p.Results;
end

% check parameter satisfy conditions or not
out_param = integraltau_g_param(out_param,default);

%% main alg
% out_param = in_param;   % save in_param to out_param
out_param.nmin=ceil((out_param.tau+1)/2)+1; % computes the minimum requirement of number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of fint is less than cost budget
out_param.tauchange=false;  % if the cone constant has been changed
xpts=linspace(0,1,out_param.nmin)'; % generate nmin number of uniformly spaced points in [0,1]
fpts=f(xpts);   % get function values at xpts
sumf=(fpts(1)+fpts(out_param.nmin))/2+sum(fpts(2:out_param.nmin-1));    % computes the sum of trapezoidal rule
ntrap=out_param.nmin-1; % number of trapezoids

while true
    %Compute approximations to the strong and weak norms
    ntrapok=true; %number of trapezoids is large enough for tau
    df=diff(fpts); %first difference of points
    Gf=sum(abs(df-(fpts(ntrap+1)-fpts(1))/ntrap)); %approx weak norm
    Ff=ntrap*(sum(abs(diff(df)))); %approx strong norm
    
    %Check necessary condition for integrand to lie in cone
    if out_param.tau*(Gf+Ff/(2*ntrap)) < Ff %f lies outside cone
        out_param.tau = 2*Ff/(Gf+Ff/(2*ntrap)); %increase tau
        out_param.tauchange=true; %flag the changed tau
        warning('MATLAB:integraltau_g:peaky','This integrand is peaky relative to tau. You may wish to increase tau for similar integrands.');
        if ntrap+1 <= (out_param.tau+1)/2 %the present ntrap is too small for tau
            inflation=ceil((out_param.tau+1)/(2*ntrap)); %prepare to increase ntrap
            ntrapok=false; %flag the number of trapezoids too small for tau
        end
    end
    
    if ntrapok %ntrap large enough for tau
        %compute a reliable error estimate
        errest=out_param.tau*Gf/(4*ntrap*(2*ntrap-out_param.tau));
        if errest <= out_param.tol %tolerance is satisfied
            Q=sumf/ntrap; %compute the integral
            break %exit while loop
        else %need to increase number of trapezoids
            %proposed inflation factor to increase ntrap by
            inflation=max(ceil(1/ntrap*sqrt(out_param.tau*Gf/(8*out_param.tol))),2);
        end
    end
    if ntrap*inflation+1 > out_param.nmax
            %cost budget does not allow intended increase in ntrap
        out_param.exceedbudget=true; %tried to exceed budget
        warning('MATLAB:integraltau_g:exceedbudget','integral_g attempts to exceed the cost budget. The answer may be unreliable.');
        inflation=floor((out_param.nmax-1)/ntrap);
            %max possible increase allowed by cost budget
        if inflation == 1 %cannot increase ntrap at all
            Q=sumf/ntrap; %compute the integral                 
            break %exit while loop
        end
    end
    
    %Increase number of sample points
    expand=repmat(xpts(1:end-1),1,inflation-1);
    addon=repmat((1:inflation-1)'/(inflation*ntrap),1,ntrap)';
    xnew=expand'+addon'; %additional x values
    ynew=f(xnew); %additional f(x) values
    xnew = [xpts(1:end-1)'; xnew];
    ynew = [fpts(1:end-1)'; ynew];
    xpts = [xnew(:); xpts(end)];
    fpts = [ynew(:); fpts(end)];
    ntrap=ntrap*inflation; %new number of trapezoids
    sumf=(fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap));
        %updated weighted sum of function values
    if out_param.exceedbudget %tried to exceed cost budget
        Q=sumf/ntrap; %compute the integral
        break; %exit while loop
    end
    
end

out_param.Q=Q;  % integral of functions
out_param.npoints=ntrap+1;  % number of points finally used
out_param.errest=errest;    % error of integral

function out_param = integraltau_g_param(out_param,default)
% let error tolerence less than 1 and greater than 0
if (out_param.tol <= 0 || out_param.tol >= 1)
    warning('Error tolerence should be less than 1 and greater than 0')
    warning('Use default error tolerence 1e-7')
    out_param.tol = default.tol;
end
% let cone condition greater or equal 2
if (out_param.tau < 2)
    warning('Cone condition should be greater or equal 2')
    warning('Use default cone condition 10')
    out_param.tau = default.tau;
end
% let cost budget be a positive integer
if (~isposint(out_param.nmax))
    warning('Cost budget should be a positive integer')
    warning('Use default cost budget 10000000')
    out_param.nmax = default.nmax;
end
