function [q,out_param] = integralab_g(varargin)
%  INTEGRALAB_G 1-D guaranteed function integration using trapezoidal rule
% 
%  Description
%
%   q = INTEGRALAB_G(f) computes q, the definite integral of function f
%   on the interval [0,1] by trapezoidal rule with 
%   in a guaranteed absolute error of 1e-6. Default starting number of
%   sample points taken is 100 and default cost budget is 1e7. Input f is a 
%   function handle. The function y = f(x) should accept a vector argument 
%   x and return a vector result y, the integrand evaluated at each element
%   of x.
%
%   q = INTEGRALAB_G(f,in_param) computes q, the definite integral of 
%   function f by trapezoidal rule within a guaranteed absolute error
%   in_param.abstol, starting number of points in_param.ninit, and cost 
%   budget in_param.nmax. If a field is not specified, the default value is
%   used.
%   
%   in_param.abstol --- absolute error tolerance required by user.
% 
%   in_param.nlo --- lowest initial number of function values used
%
%   in_param.nhi --- highest initial number of function values used
%
%   in_param.nmax --- cost budget (maximum number of function values)
% 
%   in_param.a --- low end of the integral
%
%   in_param.b --- high end of the integral
%
%   q = INTEGRAL_G(f,'abstol',abstol,'nlo',nlo,'nhi',nhi,,'a',a,'b',b,'nmax',nmax) computes
%   q, the definite integral of function f by trapezoidal rule within a 
%   guaranteed absolute error tolerance abstol, starting number of points 
%   ninit, and cost budget nmax. All three field-value pairs are optional 
%   and can be supplied.
%
%   [q, out_param] = INTEGRAL_G(f,...) returns the approximated 
%   integration q and output structure out_param, which includes the 
%   fileds in_param plus the following fields:
%
%   out_param.exceedbudget --- it is true if the algorithm tries to use 
%   more points than cost budget, false otherwise.
% 
%   out_param.tauchange --- it is true if the cone constant has been
%   changed, false otherwise. See [1] for details. If true, you may wish to
%   change the input in_param.ninit to a larger number.
% 
%   out_param.npoints --- number of points we need to 
%   reach the guaranteed absolute error tolerance abstol.
%
%   out_param.errest --- approximation error defined as the differences
%   between the true value and the approximated value of the integral.
%
%   out_param.nlo --- lowest initial number of function values
%
%   out_param.nhi --- highest initial number of function values
%
%   out_param.nmax --- cost budget (maximum number of function values)
% 
%   out_param.a --- low end of the integral
%
%   out_param.b --- high end of the integral
% 
%
%   Examples
%
%   Example 1: 
%   >> q = integralab_g(@(x) x.^2)
%   q = 0.3333
%
%
%   Example 2:
%   >> f = @(x) exp(-x.^2); q = integralab_g(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.1353
%
%
%   Example 3:
%   >> q = integralab_g()
%   Warning: Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.
%   >  In ***
%   q = 0.3333
%
%
% See also funappxab_g, cubMC_g
%
% Reference:
% [1]  N. Clancy, Y. Ding, C. Hamilton, F. J. Hickernell, and Y. Zhang, 
%      The complexity of guaranteed automatic algorithms: Cones, not
%      balls, Journal of Complexity 2013, to appear, DOI
%      10.1016/j.jco.2013.09.002.

%%
% check parameter satisfy conditions or not
[f,out_param] = integral_g_param(varargin{:});

%% main alg
out_param.tau=ceil((out_param.ninit-1)*2-1); % computes the minimum requirement of number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of q is less than cost budget
out_param.tauchange=false;  % if the cone constant has been changed
xpts=linspace(out_param.a,out_param.b,out_param.ninit)'; % generate ninit number of uniformly spaced points in [0,1]
fpts=f(xpts);   % get function values at xpts
sumf=(out_param.b-out_param.a)*(fpts(1)+fpts(out_param.ninit))/2+sum(fpts(2:out_param.ninit-1));    % computes the sum of trapezoidal rule
ntrap=out_param.ninit-1; % number of trapezoids
intervallen=out_param.b-out_param.a;

while true
    %Compute approximations to the strong and weak norms
    ntrapok=true; %number of trapezoids is large enough for ninit
    df=diff(fpts); %first difference of points
    Gf=sum(abs(df-(fpts(ntrap+1)-fpts(1))/ntrap)); %approx weak norm
    Ff=ntrap*(sum(abs(diff(df))))/intervallen; %approx strong norm
    
    %Check necessary condition for integrand to lie in cone
    if out_param.tau*(Gf+Ff*intervallen/(2*ntrap)) < Ff %f lies outside cone
        out_param.tau = 2*Ff/(Gf+Ff*intervallen/(2*ntrap)); %increase tau
        out_param.tauchange=true; %flag the changed tau
        warning('MATLAB:integral_g:peaky','This integrand is peaky relative to ninit. You may wish to increase ninit for similar integrands.');
        if ntrap+1 <= (out_param.tau+1)/2 %the present ntrap is too small for tau
            inflation=ceil((out_param.tau+1)/(2*ntrap)); %prepare to increase ntrap
            ntrapok=false; %flag the number of trapezoids too small for tau
        end
    end
    
    if ntrapok %ntrap large enough for tau
        %compute a reliable error estimate
        errest=out_param.tau*Gf*intervallen^2/(4*ntrap*(2*ntrap-out_param.tau*intervallen));
        if errest <= out_param.abstol %tolerance is satisfied
            q=sumf/ntrap; %compute the integral
            break %exit while loop
        else %need to increase number of trapezoids
            %proposed inflation factor to increase ntrap by
            inflation=max(ceil(1/ntrap*sqrt(out_param.tau*intervallen*Gf/(8*out_param.abstol))),2);
        end
    end
    if ntrap*inflation+1 > out_param.nmax
            %cost budget does not allow intended increase in ntrap
        out_param.exceedbudget=true; %tried to exceed budget
        warning('MATLAB:integral_g:exceedbudget','integral_g attempts to exceed the cost budget. The answer may be unreliable.');
        inflation=floor((out_param.nmax-1)/ntrap);
            %max possible increase allowed by cost budget
        if inflation == 1 %cannot increase ntrap at all
            q=sumf/ntrap; %compute the integral                 
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
    sumf=(out_param.b-out_param.a)*((fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap)));
        %updated weighted sum of function values
    if out_param.exceedbudget %tried to exceed cost budget
        q=sumf/ntrap; %compute the integral
        break; %exit while loop
    end
end

out_param.q=q;  % integral of functions
out_param.npoints=ntrap+1;  % number of points finally used
out_param.errest=errest;    % error of integral

function [f, out_param] = integral_g_param(varargin)
% parse the input to the integral_g function

% Default parameter values
default.abstol  = 1e-6;
default.nmax  = 1e7;
default.nlo = 10;
default.nhi = 1000;
default.a = 0;
default.b = 1;


if isempty(varargin)
    help integral_g
    warning('Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
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
    out_param.nmax = default.nmax;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.a = default.a;
    out_param.b = default.b;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'nlo',default.nlo,@isnumeric);
        addParamValue(p,'nhi',default.nhi,@isnumeric);
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
if (~isposint(out_param.nlo))
    if isposge3(out_param.nlo)
        warning('MATLAB:integralab_g:lowinitnotint',['Lowest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo))])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:integralab_g:lowinitlt3',['Lowest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~isposint(out_param.nhi))
    if isposge3(out_param.nhi)
        warning('MATLAB:integralab_g:highinitnotint',['Highest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:integralab_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    if isposge3(out_param.nhi)
        warning('MATLAB:integralab_g:nlobtnhi',['Highest initial number of points should be at least equal to to lowest initial number of points.' ...
            ' Using ', num2str(ceil(out_param.nhi)), ' as nlo'])
        out_param.nlo = ceil(out_param.nhi);
    else
        warning('MATLAB:integralab_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+(out_param.b-out_param.a))));
if (~isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('MATLAB:integralab_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:integralab_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end


%% doctest results%%
% doctest integralab_g
% TAP version 13
% 1..3
% ok 1 - "q = integralab_g(@(x) x.^2)"
% ok 2 - "f = @(x) exp(-x.^2); q = integralab_g(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)"
% ok 3 - "q = integralab_g()"
