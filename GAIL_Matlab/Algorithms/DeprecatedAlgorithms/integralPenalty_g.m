function [q,out_param] = integralPenalty_g(varargin)
%INTEGRALPENALTY_G 1-D guaranteed function integration using trapezoidal rule
%
%   q = INTEGRALPENALTY_G(f) computes q, the definite integral of function f on
%   the interval [a,b] by trapezoidal rule with in a guaranteed absolute
%   error of 1e-6. Default starting number of sample points taken is 100
%   and default cost budget is 1e7. Input f is a function handle. The
%   function y = f(x) should accept a vector argument x and return a vector
%   result y, the integrand evaluated at each element of x.
%
%   q = INTEGRALPENALTY_G(f,a,b,abstol) computes q, the definite integral of
%   function f on the finite interval [a,b] by trapezoidal rule with the
%   ordered input parameters, and guaranteed absolute error tolerance
%   abstol.
%
%   q = INTEGRALPENALTY_G(f,'a',a,'b',b,'abstol',abstol) computes q, the definite
%   integral of function f on the finite interval [a,b] by trapezoidal rule
%   within a guaranteed absolute error tolerance abstol. All four
%   field-value pairs are optional and can be supplied.
%
%   q = INTEGRALPENALTY_G(f,in_param) computes q, the definite integral of
%   function f by trapezoidal rule within a guaranteed absolute error
%   in_param.abstol. If a field is not specified, the default value is
%   used.
%
%   [q, out_param] = INTEGRALPENALTY_G(f,...) returns the approximated 
%   integration q and output structure out_param.
%
%   Input Arguments
%
%     f --- input function
%
%     in_param.a --- left end of the integral, default value is 0
%
%     in_param.b --- right end of the integral, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default value
%     is 1e-6
% 
%   Optional Input Arguments
%
%     in_param.nlo --- lowest initial number of function values used, default
%     value is 10
%
%     in_param.nhi --- highest initial number of function values used,
%     default value is 1000
%
%     in_param.nmax --- cost budget (maximum number of function values),
%     default value is 1e7
%
%     in_param.maxiter --- max number of iterations, default value is 1000
% 
%   Output Arguments
%
%     q --- approximated integral
%
%     out_param.f --- input function
%
%     out_param.a --- low end of the integral
%
%     out_param.b --- high end of the integral
%
%     out_param.abstol --- guaranteed absolute error tolerance
% 
%     out_param.nlo --- lowest initial number of function values
%
%     out_param.nhi --- highest initial number of function values
%
%     out_param.nmax --- cost budget (maximum number of function values)
%
%     out_param.maxiter --- max number of iterations
%
%     out_param.ninit --- initial number of points we use, computed by nlo
%     and nhi
%
%     out_param.tauchange --- it is true if the cone constant has been
%     changed, false otherwise. See [1] for details. If true, you may wish to
%     change the input in_param.ninit to a larger number.
% 
%     out_param.tauchange --- it is true if the cone constant has been
%     changed, false otherwise. See [1] for details. If true, you may wish to
%     change the input in_param.ninit to a larger number.
% 
%     out_param.iter --- number of iterations
%
%     out_param.npoints --- number of points we need to 
%     reach the guaranteed absolute error tolerance abstol.
%
%     out_param.errest --- approximation error defined as the differences
%     between the true value and the approximated value of the integral.
%
%     out_param.nstar --- final value of the parameter defining the cone of
%     functions for which this algorithm is guaranteed; nstar = ninit-2
%     initially and is increased as necessary
%
%     out_param.exit --- the state of program when exiting
%            0  Success
%            1  Number of points used is greater than out_param.nmax
%            2  Number of iterations is greater than out_param.maxiter
%
%  Guarantee
%    
%  If the function to be integrated, f, satisfies the cone condition
%                          nstar   ||     f(b)-f(a)  ||
%      ||f''||        <=  -------- ||f'- ----------- ||
%             1           2(b - a) ||       b - a    ||1,
%  then the q output by this algorithm is guaranteed to satisfy
%      |\int_{a}^{b} f(x) dx - q | <= abstol,
%  provided the flag exceedbudget = 0. And the upper bound of the cost is
%          ________________________ 
%         /   nstar*(b-a)^2 Var(f') 
%        / ------------------------ + 2 nstar + 4
%      \/          2 abstol
%
%  Examples
%
%   Example 1: 
%   >> q = integralPenalty_g(@(x) x.^2)
%   q = 0.3333
%
%
%   Example 2:
%   >> f = @(x) exp(-x.^2); q = integralPenalty_g(f,'a',1,'b',2,'abstol',1e-5)
%   q = 0.1353
%
%
%   Example 3:
%   >> q = integralPenalty_g()
%   Warning: Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.
%   >  In ***
%   q = 0.3333
%
%
%   See also INTEGRAL, QUAD, FUNAPPX_G, MEANMC_G, CUBMC_G, FUNMIN_G
%
%  References
%
%   [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%   Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
%   Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
%
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.2)
%   [MATLAB Software], 2017. Available from http://gailgithub.github.io/GAIL_Dev/
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
%   Workshop On Sustainable Software for Science: Practice And Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%


%%
% check parameter satisfy conditions or not
[f,out_param, flip] = integral_g_param(varargin{:});

%% main alg
out_param.tau=ceil((out_param.ninit-1)*2-1); % computes the minimum requirement of number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of q is less than cost budget
out_param.tauchange=false;  % if the cone constant has been changed
xpts=linspace(out_param.a,out_param.b,out_param.ninit)'; % generate ninit number of uniformly spaced points in [0,1]
fpts=f(xpts);   % get function values at xpts
intervallen=out_param.b-out_param.a;
sumf=intervallen*(fpts(1)+fpts(out_param.ninit))/2+sum(fpts(2:out_param.ninit-1));    % computes the sum of trapezoidal rule
ntrap=out_param.ninit-1; % number of trapezoids
iter=0;

if intervallen

    while true
        iter=iter+1;
        %Compute approximations to the strong and weak norms
        ntrapok=true; %number of trapezoids is large enough for ninit
        df=diff(fpts); %first difference of points
        Gf=sum(abs(df-(fpts(ntrap+1)-fpts(1))/ntrap)); %approx weak norm
        Ff=ntrap*(sum(abs(diff(df))))/intervallen; %approx strong norm

        %Check necessary condition for integrand to lie in cone
        if out_param.tau*(Gf/intervallen+Ff/(2*ntrap)) < Ff %f lies outside cone
            out_param.tau = 2*Ff/(Gf/intervallen+Ff/(2*ntrap)); %increase tau
            out_param.tauchange=true; %flag the changed tau
            warning('GAIL:integralPenalty_g:peaky','This integrand is peaky relative to ninit. You may wish to increase ninit for similar integrands.');
            if ntrap+1 <= (out_param.tau+1)/2 %the present ntrap is too small for tau
                inflation=ceil((out_param.tau+1)/(2*ntrap)); %prepare to increase ntrap
                ntrapok=false; %flag the number of trapezoids too small for tau
            end
        end

        if ntrapok %ntrap large enough for tau
            %compute a reliable error estimate
            errest=out_param.tau*Gf*intervallen/(4*ntrap*(2*ntrap-out_param.tau));
            if errest <= out_param.abstol || ... %tolerance is satisfied
                  out_param.exceedbudget %or tried to exceed cost budget
                q=sumf/ntrap; %compute the integral
                %keyboard
                break %exit while loop
           else %need to increase number of trapezoids
                %proposed inflation factor to increase ntrap by
                inflation=max(ceil(1/ntrap*sqrt(out_param.tau*intervallen*Gf/(8*out_param.abstol))),2);
            end
        end
        if ntrap*inflation+1 > out_param.nmax
            %cost budget does not allow intended increase in ntrap
            out_param.exit=1; %tried to exceed budget
            warning('GAIL:integralPenalty_g:exceedbudget','integralPenalty_g attempts to exceed the cost budget. The answer may be unreliable.');
            inflation=floor((out_param.nmax-1)/ntrap);
                %max possible increase allowed by cost budget
            if inflation == 1 %cannot increase ntrap at all
                q=sumf/ntrap; %compute the integral                 
                break %exit while loop
            end
        end

        %Increase number of sample points
        expand=repmat(xpts(1:end-1),1,inflation-1);
        addon=repmat((1:inflation-1)*(intervallen/(inflation*ntrap)),ntrap,1);
        xnew=expand+addon; %additional x values
        xnew=xnew';
        ynew=f(xnew); %additional f(x) values
        xnew = [xpts(1:end-1)'; xnew];
        ynew = [fpts(1:end-1)'; ynew];
        xpts = [xnew(:); xpts(end)];
        fpts = [ynew(:); fpts(end)];
        ntrap=ntrap*inflation; %new number of trapezoids
        sumf=intervallen*((fpts(1)+fpts(ntrap+1))/2+sum(fpts(2:ntrap)));
        if(iter> out_param.maxiter)
           out_param.exit = 2;
           warning('GAIL:integralPenalty_g:exceediter',' Iteration exceeds max iteration ')
           break;
        end;

    end
elseif intervallen == 0
    q = 0;
    errest = 0;
end
if flip==1
    q = -1*q;
end

out_param.iter=iter;
out_param.q=q;  % integral of functions
out_param.npoints=ntrap+1;  % number of points finally used
out_param.errest=errest;    % error of integral

function [f, out_param, flip] = integral_g_param(varargin)
% parse the input to the integralPenalty_g function

% Default parameter values
default.abstol  = 1e-6;
default.nmax  = 1e7;
default.nlo = 10;
default.nhi = 1000;
default.a = 0;
default.b = 1;
default.maxiter = 1000;
% if a<b, flip = 0; if a>b, flip = 1;
flip = 0;

MATLABVERSION= gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;

if isempty(varargin)
    help integralPenalty_g
    warning('GAIL:integralPenalty_g:nofunction','Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
    out_param.f=f;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:integralPenalty_g:notfunction','Function f must be a function handle. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
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
    out_param.nmax = default.nmax;
    out_param.maxiter = default.maxiter;
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
        addOptional(p,'nmax',default.nmax,@isnumeric);
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
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
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('GAIL:integralPenalty_g:anoinfinity',['a can not be infinity nor NaN. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('GAIL:integralPenalty_g:bnoinfinity',['b can not be infinity not Nan. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;
if (out_param.b < out_param.a)
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
    flip=1;
end

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('GAIL:integralPenalty_g:abstolnonpos', ['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let initial number of points be a positive integer
if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('GAIL:integralPenalty_g:lowinitnotint',['Lowest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo))])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('GAIL:integralPenalty_g:lowinitlt3',['Lowest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('GAIL:integralPenalty_g:highinitnotint',['Highest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('GAIL:integralPenalty_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    if gail.isposge3(out_param.nhi)
        warning('GAIL:integralPenalty_g:nlobtnhi',['Highest initial number of points should be at least equal to to lowest initial number of points.' ...
            ' Using ', num2str(ceil(out_param.nlo)), ' as ninit'])
        out_param.nhi = ceil(out_param.nlo);
    else
        warning('GAIL:integralPenalty_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

out_param.ninit = max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+(out_param.b-out_param.a)))),3);

if (~gail.isposint(out_param.nmax))
    if gail.ispositive(out_param.nmax)
        warning('GAIL:integralPenalty_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:integralPenalty_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

if (~gail.isposint(out_param.maxiter))
    if gail.ispositive(out_param.maxiter)
        warning('GAIL:integralPenalty_g:maxiternotint',['Max number of iterations should be a positive integer.' ...
            ' Using max number of iterations as  ', num2str(ceil(out_param.maxiter))])
        out_param.maxiter = ceil(out_param.maxiter);
    else
        warning('GAIL:integralPenalty_g:budgetisneg',['Max number of iterations should be a positive integer.' ...
            ' Using max number of iterations as ' int2str(default.maxiter)])
        out_param.maxiter = default.maxiter;
    end;
end

