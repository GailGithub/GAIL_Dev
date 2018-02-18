function [q,out_param] = integral_g(varargin)
%INTEGRAL_G 1-D guaranteed function integration using trapezoidal rule
%
%   q = INTEGRAL_G(f) computes q, the definite integral of function f on
%   the interval [a,b] by trapezoidal rule to within a specified
%   generalized error tolerance, tolfun := max(abstol, reltol*|q|). Default
%   starting number of sample points taken is 100 and default cost budget
%   is 1e7. Input f is a function handle. The function y = f(x) should
%   accept a vector argument x and return a vector result y, the integrand
%   evaluated at each element of x.
%
%   q = INTEGRAL_G(f,a,b,abstol,reltol) computes q, the definite integral
%   of function f on the finite interval [a,b] by trapezoidal rule with the
%   ordered input parameters, and to within the guaranteed generalized
%   error tolerance, tolfun.
%
%   q = INTEGRAL_G(f,'a',a,'b',b,'abstol',abstol,'reltol',reltol) computes
%   q, the definite integral of function f on the finite interval [a,b] by
%   trapezoidal rule to within a guaranteed generalized error tolerance,
%   tolfun. All five field-value pairs are optional and can be supplied by
%   default values.
%
%   q = INTEGRAL_G(f,in_param) computes q, the definite integral of
%   function f by trapezoidal rule to within a guaranteed guaranteed
%   generalized error tolerance, tolfun := max(in_param.abstol,
%   in_param.reltol * |q|). If a field is not specified, the default value
%   is used.
%
%   [q, out_param] = INTEGRAL_G(f,...) returns the approximated integration
%   q and output structure out_param.
%
%   Input Arguments
%
%
%     f --- input function
%
%     in_param.a --- left end of the integral, default value is 0
%
%     in_param.b --- right end of the integral, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6
%
%     in_param.reltol --- guaranteed relative error tolerance, default
%     value is 0
%
%  Optional Input Arguments (Recommended not to change very often)
%
%
%     in_param.nlo --- lowest initial number of function values used,
%     default value is 10
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
%     out_param.reltol --- guaranteed relative error tolerance
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
%     out_param.exceedbudget --- it is true if the algorithm tries to use
%     more points than cost budget, false otherwise.
%
%     out_param.tauchange --- it is true if the cone constant has been
%     changed, false otherwise. See [1] for details. If true, you may wish
%     to change the input in_param.ninit to a larger number.
%
%     out_param.iter --- number of iterations
%
%     out_param.npoints --- number of points we need to reach the
%     generalized error tolerance, tolfun.
%
%     out_param.errest --- approximation error defined as the differences
%     between the true value and the approximated value of the integral.
%
%     out_param.nstar --- final value of the parameter defining the cone of
%     functions for which this algorithm is guaranteed; nstar = ninit-2
%     initially and is increased as necessary
%
%  Guarantee
%
% Please check the details of the guarantee in [1].
%
%
%  Examples
%
%   Example 1:
%   >> q = integral_g(@(x) x.^2)
%   q = 0.3333
%
%
%   Example 2:
%   >> f = @(x) exp(-x.^2); 
%   >> q = integral_g(f,'a',1,'b',2,'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)
%   q = 0.1353
%
%
%   Example 3:
%   >> f = @(x) exp(-x.^2); in_param.abstol = 1e-6; in_param.reltol = 1e-5;
%   >> q = integral_g(f,in_param)
%   q = 0.7468
%
%
%   Example 4:
%   >> q = integral_g()
%   Warning: Function f must be a function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2).
%   >  In ***
%   q = 0.1772
%
%
%   See also INTEGRAL, QUAD, MEANMC_G, CUBMC_G, CUBSOBOL_G, CUBLATTICE_G,
%   FUNAPPX_G,  FUNMIN_G
%
%  References
%
%   [1] Fred J. Hickernell, Martha Razo, and Sunny Yun, "Reliable Adaptive
%   Numerical Integration", 2015+, working.
%
%   [2] Jiazhen Lu, "Adaptive Quadrature with a General Error Criterion,"
%   MS thesis, Illinois Institute of Technology, 2018.
%
%   [3] Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%   Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
%   Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
%
%   [4] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.2) [MATLAB Software], 2017. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%


%%
% check parameter satisfy conditions or not
[f,out_param, flip] = integral_g_param(varargin{:});

%% main alg
out_param.tau=ceil((out_param.ninit-1)*2-1); % computes the minimum number of points to start
out_param.exceedbudget=false;   % if the number of points used in the calculation of q is less than cost budget
out_param.conechange=false;  % if the cone constant has been changed
ntrap=out_param.ninit-1; % number of trapezoids
intervallen=out_param.b-out_param.a; % length of integration interval (>=0)
Varfp=zeros(10,1); %initialize vector of approximations to Var(f')
Varfpup=[inf; zeros(10,1)]; %initialize vector of upper bounds to Var(f')
ii=1; %index

if intervallen>0
    steplen=intervallen/ntrap; % length of subinterval
    steplenvec=zeros(10,1); %vector to record subinterval lengths
    hcut=2*intervallen/(ntrap-1); %minimum interval size
    inflatelim=1.5;
    xpts=(out_param.a:steplen:out_param.b)'; % generate ninit uniformly spaced points in [a,b]
    fpts=f(xpts);   % get function values at xpts
    sumf=(fpts(1)+fpts(out_param.ninit))/2+sum(fpts(2:ntrap)); % computes the sum of trapezoidal rule
    while true
        steplenvec(ii)=steplen;

        %Compute approximation to Var(f')
        Varfp(ii)=sum(abs(diff(fpts,2)))/steplen; %approx Var(f')
        Varfpup(ii+1)=min(Varfpup(ii),Varfp(ii)*inflatelim*hcut/(hcut-2*steplen));
            %update upper bound on Var(f')

        %Check necessary condition for integrand to lie in cone
        if Varfp(ii) > Varfpup(ii+1) %f lies outside cone
            %Decrease hcut
            tempa=1-inflatelim*Varfp(1:ii)/Varfp(ii);
            whpos=find(tempa>0,1,'first');
            hcut=min(steplenvec(whpos:ii)./tempa(whpos:ii));
            out_param.conechange=true; %flag the changed tau
            warning('GAIL:integral_g:spiky','This integrand is spiky relative to ninit. You may wish to increase ninit for similar integrands.');

            %Update Varfpup
            Varfpup(whpos+1:ii+1)=min(Varfp(whpos:ii)*inflatelim*hcut/(hcut-2*steplenvec(whpos:ii)));
        end

        %Check error
        errest=Varfpup(ii+1)*(steplen.^2)/8;
        t=sumf*steplen; %compute the integral
        tminus=t-errest; 
        tplus=t+errest;%compute extreme value
        delta=(max(out_param.abstol,out_param.reltol*abs(tplus))+max(out_param.abstol,out_param.reltol*abs(tminus)))/2;
        if errest <= delta %tolerance is satisfied
            q=(tminus*max(out_param.abstol,out_param.reltol*abs(tplus))+tplus*max(out_param.abstol,out_param.reltol*abs(tminus)))/...
                (max(out_param.abstol,out_param.reltol*abs(tplus))+max(out_param.abstol,out_param.reltol*abs(tminus)));
            %keyboard
            break %exit while loop
        else %need to increase number of trapezoids
            %proposed inflation factor to increase ntrap by
            inflation=2;
        end
        if ntrap*inflation+1 > out_param.nmax
                %cost budget does not allow intended increase in ntrap
            out_param.exceedbudget=true; %tried to exceed budget
            warning('GAIL:integral_g:exceedbudget','integral_g attempts to exceed the cost budget. The answer may be unreliable.');
            inflation=floor((out_param.nmax-1)/ntrap);
                %max possible increase allowed by cost budget
            if inflation == 1 %cannot increase ntrap at all
                q=sumf*steplen; %compute the integral
                break %exit while loop
            end
        end

        %Increase number of sample points
        steplen=steplen/inflation;
        ntrapnew=ntrap*inflation;
        xnew=bsxfun(@plus,(1:inflation-1)'*steplen,xpts(1:ntrap)'); %additional x values
        ynew=f(xnew); %additional f(x) values
        sumf=sumf+sum(ynew(:)); %updated weighted sum of function values
        xpts = [reshape([xpts(1:ntrap)'; xnew],ntrapnew,1); xpts(ntrap+1)];
        fpts = [reshape([fpts(1:ntrap)'; ynew],ntrapnew,1); fpts(ntrap+1)];
        ntrap=ntrapnew; %new number of trapezoids
        ii=ii+1; %increment the counter

    end
else
    q = 0;
    errest = 0;
end
if flip
    q = -1*q;
end
out_param.npoints=ntrap+1;  % number of points finally used
out_param.errest=errest;    % error of integral
out_param.VarfpCI=[Varfp(ii) Varfpup(ii+1)];

% reorder fields in out_param
out_param = orderfields(out_param, ...
           {'f', 'a', 'b', 'abstol', 'reltol', 'nlo', 'nhi', 'nmax', ...
            'ninit', 'tau', 'exceedbudget', 'conechange', 'npoints',...
            'errest', 'VarfpCI'});
            

function [f, out_param, flip] = integral_g_param(varargin)
% parse the input to the integral_g function

% Default parameter values
default.abstol  = 1e-6;
default.reltol = 0;
default.nmax  = 1e7;
default.nlo = 10;
default.nhi = 1000;
default.a = 0;
default.b = 1;
% if a<b, flip = 0; if a>b, flip = 1;
flip = false;

if isempty(varargin)
    help integral_g
    warning('GAIL:integral_g:nofunction','Function f must be specified. Now GAIL is giving you a toy example of f(x)=x^2.')
    f = @(x) x.^2;
    out_param.f=f;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:integral_g:notfunction','Function f must be a function handle. Now GAIL is giving you a toy example of f(x)=x^2.')
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
    out_param.reltol = default.reltol;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'nlo',default.nlo,@isnumeric);
        addParamValue(p,'nhi',default.nhi,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('GAIL:integral_g:anoinfinity',['a can not be infinity nor NaN. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('GAIL:integral_g:bnoinfinity',['b can not be infinity not Nan. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;
if (out_param.b < out_param.a)
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
    flip=1;
end

% let absolute error tolerance be greater than 0
if (out_param.abstol < 0)
    warning('GAIL:integral_g:abstolneg',['Error tolerance should be greater than or equal to 0.' ...
            ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let relative error tolernace be between 0 and 1
if (out_param.reltol < 0)
    warning('GAIL:integral_g:abstolnonpos',['Relative error tolerance should be between 0 and 1.' ...
        ' Using default error tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end
if (out_param.reltol >= 1)
    warning('GAIL:integral_g:abstolgt1',['Relative error tolerance should be between 0 and 1.' ...
        ' Using default error tolerance ' num2str(default.reltol)])
    out_param.reltol = default.reltol;
end
% let generalized error tolernace be greater than 0
tol = max(out_param.abstol, out_param.reltol);
if (tol <= 0)
    warning('GAIL:integral_g:gentolnonpos',['Generalized error tolerance should be greater than 0.' ...
            ' Using default absolute error tolerance ' num2str(default.abstol) ...
            ' and relative error tolerance ' num2str(default.reltol)])
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
end

% let initial number of points be a positive integer
if (~gail.isposint(out_param.nlo))
    if isposge3(out_param.nlo)
        warning('GAIL:integral_g:lowinitnotint',['Lowest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo))])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('GAIL:integral_g:lowinitlt3',['Lowest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~gail.isposint(out_param.nhi))
    if isposge3(out_param.nhi)
        warning('GAIL:integral_g:highinitnotint',['Highest initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('GAIL:integral_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    if isposge3(out_param.nhi)
        warning('GAIL:integral_g:nlobtnhi',['Highest initial number of points should be at least equal to to lowest initial number of points.' ...
            ' Using ', num2str(ceil(out_param.nhi)), ' as nlo'])
        out_param.nlo = ceil(out_param.nhi);
    else
        warning('GAIL:integral_g:highinitlt3',['Highest initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end

out_param.ninit = max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+(out_param.b-out_param.a)))),3);
if (~gail.isposint(out_param.nmax))
    if ispositive(out_param.nmax)
        warning('GAIL:integral_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:integral_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end
