function [q,out_param] = integralsim_g(varargin)
%INTEGRALSIM_G 1-D guaranteed function integration using Simpson's rule
% 
%  Description
%
%   q = INTEGRALSIM_G(f) computes q, the definite integral of function f on
%   the interval [0,1] by Simpson's rule with in a guaranteed absolute
%   error of 1e-6. Default starting number of sample points taken is 100
%   and default cost budget is 1e7. Input f is a function handle. The
%   function y = f(x) should accept a vector argument x and return a vector
%   result y, the integrand evaluated at each element of x.
%
%   q = INTEGRAL_G(f,abstol,ninit, nmax,maxiter) computes q, the definite 
%   integral of function f by Simpson's rule with the ordered input 
%   parameters, guaranteed absolute error tolerance abstol, starting number
%   of points ninit, cost budget nmax and max number of iterations maxiter.
%
%   q = INTEGRAL_G(f,'abstol',abstol,'ninit',ninit,'nmax',nmax,'maxiter',maxiter) computes
%   q, the definite integral of function f by Simpson's rule within a 
%   guaranteed absolute error tolerance abstol, starting number of points 
%   ninit, cost budget nmax and max number of iterations maxiter. All three field-value pairs are optional 
%   and can be supplied.
%
%   q = INTEGRAL_G(f,in_param) computes q, the definite integral of 
%   function f by Simpson's rule within a guaranteed absolute error
%   in_param.abstol, starting number of points in_param.ninit, cost 
%   budget in_param.nmax and max number of iterations maxiter. If a field is not specified, the default value is
%   used.
%   
%   [q, out_param] = INTEGRAL_G(f,...) returns the approximated 
%   integration q and output structure out_param.
%
%   Input Arguments
%
%   in_param.abstol --- absolute error tolerance required by user.
% 
%   in_param.ninit --- initial number of function values used
%
%   in_param.nmax --- cost budget (maximum number of function values)
%
%   in_param.maxiter --- max number of iterations, default value is 1000
% 
%   [q, out_param] = INTEGRAL_G(f,...) returns the approximated 
%   integration q and output structure out_param, which includes the 
%   fileds in_param plus the following fields:
%
%   out_param.tauchange --- it is true if the cone constant has been
%   changed, false otherwise. If true, you may wish to
%   change the input in_param.ninit to a larger number.
% 
%   out_param.npoints --- number of points we need to 
%   reach the guaranteed absolute error tolerance abstol.
%
%   out_param.errest --- approximation error defined as the differences
%   between the true value and the approximated value of the integral.
% 
%   out_param.maxiter --- max number of iterations
%
%   out_param.iter --- number of iterations
%
%   out_param.exit --- the state of program when exiting
%            0  Success
%            1  Number of points used is greater than out_param.nmax
%            2  Number of iterations is greater than out_param.maxiter
%
%   Examples
%
%   Example 1: 
%   >> q = integralsim_g(@(x) x.^2)
%   q = 0.3333
%
%
%   Example 2:
%   >> f = @(x) exp(-x.^2); q = integralsim_g(f,'abstol',1e-5,'ninit',53,'nmax',1e7,'maxiter',500)
%   q = 0.7468
%
%
%   Example 3:
%   >> q = integralsim_g()
%   Warning: Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.
%   >  In ***
%   q = 0.3333
%
%
%   See also FUNAPPX_G, MEANMC_G, CUBMC_G, FUNMIN_G
%
%   References
%
%   [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%   Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
%   Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
%            
%   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1) [MATLAB
%   Software], 2015. Available from http://code.google.com/p/gail/
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
%


% check parameter satisfy conditions or not
[f,out_param] = integralsim_g_param(varargin{:});
% MATLABVERSION= gail.matlab_version;

%% main alg
out_param.tau=out_param.ninit-1; % computes the minimum requirement of number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of q is less than cost budget
out_param.tauchange=false;  % if the cone constant has been changed
xpts=linspace(0,1,out_param.ninit)'; % generate ninit number of uniformly spaced points in [0,1]
fpts=f(xpts);   % get function values at xpts
sum1=reshape(fpts(2:out_param.ninit),2,(out_param.ninit-1)/2); %compute the 4-time part of Simpson's rule
sumf=(fpts(1)+fpts(out_param.ninit))+2*sum(fpts(2:out_param.ninit-1))+2*sum(sum1(1,:));    % computes the sum of Simpson's rule
nint=out_param.ninit-1; % number of intevals
iter = 0; % number of iteration used

while true
    iter = iter +1;
    %Compute approximations to the strong and weak norms
    nintok=true; %ninit is large enough for tau
    df=diff(fpts); %first difference of points
    df1=reshape(df,2,length(df)/2); %matrix operation
    df1=df1(2,:)-df1(1,:); %matrix operation
    Gf=sum(abs(2*nint*df1-8*(fpts(nint+1)-2*f(0.5)+fpts(1))/nint)); %approx weak norm   
    Ff=nint^2*(sum(abs(diff(diff(df))))); %approx strong norm
%     Gf=sum(abs(df-(fpts(nint+1)-fpts(1))/nint)); %approx weak norm
%     Ff=nint*(sum(abs(diff(df)))); %approx strong norm
    
    %Check necessary condition for integrand to lie in cone
    if out_param.tau*(Gf+Ff/(2*nint)) < Ff %f lies outside cone
        out_param.tau = 2*Ff/(Gf+Ff/(2*nint)); %increase tau
        out_param.tauchange=true; %flag the changed tau
        warning('GAIL:integralsim_g:peaky','This integrand is peaky relative to ninit. You may wish to increase ninit for similar integrands.');
        if nint+1 <= (out_param.tau+1)/2 %the present ntrap is too small for tau
            inflation=ceil((out_param.tau+1)/(2*nint)); %prepare to increase ntrap
            nintok=false; %flag the number of trapezoids too small for tau
        end
    end
    
    if nintok %ntrap large enough for tau
        %compute a reliable error estimate
        errest=out_param.tau^2*Gf/(36*nint.^3*(2*nint-out_param.tau));
        if errest <= out_param.abstol %tolerance is satisfied
            q=sumf/nint/3; %compute the integral
            break %exit while loop
        else %need to increase number of trapezoids
            %proposed inflation factor to increase ntrap by
            inflation=max(ceil(1/nint*sqrt(out_param.tau*Gf/(8*out_param.abstol))),2);
        end
    end
    if nint*inflation+1 > out_param.nmax
            %cost budget does not allow intended increase in ntrap
        out_param.exit=1; %tried to exceed budget
        warning('GAIL:integralsim_g:exceedbudget','integralsim_g attempts to exceed the cost budget. The answer may be unreliable.');
        inflation=floor((out_param.nmax-1)/nint);
            %max possible increase allowed by cost budget
        if inflation == 1 %cannot increase ntrap at all
            q=sumf/nint/3; %compute the integral                 
            break %exit while loop
        end
    end
    
    %Increase number of sample points
%     expand=repmat(xpts(1:end-1),1,inflation-1);
%     addon=repmat((1:inflation-1)'/(inflation*nint),1,nint)';
%     xnew=expand'+addon'; %additional x values
%     ynew=f(xnew); %additional f(x) values
%     xnew = [xpts(1:end-1)'; xnew];
%     ynew = [fpts(1:end-1)'; ynew];
%     xpts = [xnew(:); xpts(end)];
%     fpts = [ynew(:); fpts(end)];
    nint=nint*inflation; %new number of trapezoids
    if (nint+2)/2-ceil((nint-1)/2) == 0  %check if new number of points is odd
        nint=nint+1;
    end
    xpts=linspace(0,1,nint+1)'; % generate ninit number of uniformly spaced points in [0,1]
    fpts=f(xpts);   % get function values at xpts
    sum1=reshape(fpts(2:nint+1),2,nint/2);
    sumf=(fpts(1)+fpts(nint+1))+2*sum(fpts(2:nint))+2*sum(sum1(1,:));    
        %updated weighted sum of function values
    if out_param.exceedbudget %tried to exceed cost budget
        q=sumf/nint/3; %compute the integral
        break; %exit while loop
    end
    if(iter> out_param.maxiter)
        out_param.exit = 2;
        warning('GAIL:integralsim_g:exceediter',' Iteration exceeds max iteration ')
        break;
    end;

end

out_param.iter = iter;
out_param.q=q;  % integral of functions
out_param.npoints=nint+1;  % number of points finally used
out_param.errest=errest;    % error of integral

function [f, out_param] = integralsim_g_param(varargin)
% parse the input to the integral_g function

% Default parameter values
default.abstol  = 1e-6;
default.ninit  = 53; % must be an odd number
default.nmax  = 1e7;
default.maxiter = 1000;

MATLABVERSION= gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;

if isempty(varargin)
    help integralsim_g
    warning('GAIL:integralsim_g:nofunction','Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:integralsim_g:notfunction','Function f must be a function handle. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = varargin{1};
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
    out_param.abstol = default.abstol;
    out_param.ninit = default.ninit;
    out_param.nmax = default.nmax;
    out_param.maxiter = default.maxiter;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'ninit',default.ninit,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'ninit',default.ninit,@isnumeric);
        f_addParamVal(p,'nmax',default.nmax,@isnumeric);
        f_addParamVal(p,'maxiter',default.maxiter,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('GAIL:integralsim_g:abstolnonpos',['Error tolerance should be greater than 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let ninit be an odd number
if ((out_param.ninit+1)/2-ceil(out_param.ninit/2) ~= 0 )
    warning('GAIL:integralsim_g:initnotodd',['Initial number of points must be an odd number.' ...
            ' Using default number of points ' num2str(default.ninit)])
    out_param.ninit = default.ninit;
end
% let initial number of points be a positive integer
if (~gail.isposint(out_param.ninit))
    if gail.isposge3(out_param.ninit)
        warning('GAIL:integralsim_g:initnotint',['Initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.ninit))])
        out_param.ninit = ceil(out_param.ninit);
    else
        warning('GAIL:integralsim_g:initlt3',['Initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.ninit)])
        out_param.ninit = default.ninit;
    end
end
% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.ispositive(out_param.nmax)
        warning('GAIL:integralsim_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:integralsim_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end
