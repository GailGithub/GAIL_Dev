function [pp,out_param]=funappxlocal_g(varargin)
%FUNAPPXLOCAL_G 1-D guaranteed function recovery on a closed interval [a,b]
%
%   pp = FUNAPPXLOCAL_G(f) approximates function f on the default interval [0,1]
%   by a piecewise polynomial structure pp within the guaranteed absolute
%   error tolerance of 1e-6. Input f is a function handle. The
%   statement y = f(x) should accept a vector argument x and return a
%   vector y of function values that is of the same size as x. Output pp
%   may be evaluated via PPVAL.
%   
%   pp = FUNAPPXLOCAL_G(f,a,b,abstol,nlo,nhi) for a given function f and
%   the ordered input parameters that define the finite interval [a,b], a
%   guaranteed absolute error tolerance abstol, a lower bound of initial
%   number of points nlo, and an upper bound of initial number of points
%   nhi.
%
%   pp = FUNAPPXLOCAL_G(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi)
%   recovers function f on the finite interval [a,b], given a guaranteed
%   absolute error tolerance abstol, a lower bound of initial number of
%   points nlo, and an upper bound of initial number of points nhi. All
%   five field-value pairs are optional and can be supplied in different
%   order.
%
%   pp = FUNAPPXLOCAL_G(f,in_param) recovers function f on the finite
%   interval [in_param.a,in_param.b], given a guaranteed absolute error
%   tolerance in_param.abstol, a lower bound of initial number of points
%   in_param.nlo, and an upper bound of initial number of points
%   in_param.nhi. If a field is not specified, the default value is used.
%
%     in_param.a --- left end point of interval, default value is 0
%
%     in_param.b --- right end point of interval, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6
%
%     in_param.nlo --- lower bound of initial number of points we used,
%     default value is 9
%
%     in_param.nhi --- upper bound of initial number of points we used,
%     default value is 100
%
%   [pp, out_param] = FUNAPPXLOCAL_G(f,...) returns a piecewise polynomial
%   structure pp and an output structure out_param, which have the
%   following fields:
%
%     pp.form --- pp means piecewise polynomials
%
%     pp.breaks --- show the location of interpolation points
%
%     pp.coefs --- coefficients for piecewise linear polynomials
%
%     pp.pieces --- number of piecewise linear polynomials
%
%     pp.order --- be 2 as we use piecewise linear polynomials
%
%     pp.dim --- be 1 as we do univariate approximation
%
%     pp.orient --- always be 'first'
%
%     out_param.ninit --- initial number of points we use
%
%     out_param.npoints --- number of points we need to reach the 
%     guaranteed absolute error tolerance
%
%     out_param.errorbound --- an upper bound of the absolute error
%
%     out_param.nstar --- final value of the parameter defining the cone of
%     functions for which this algorithm is guaranteed for each
%     subinterval; nstar = ninit-2 initially
%
%     out_param.a --- left end point of interval
%
%     out_param.b --- right end point of interval
%
%     out_param.abstol --- guaranteed absolute error tolerance
%
%     out_param.nlo --- a lower bound of initial number of points we use
%
%     out_param.nhi --- an upper bound of initial number of points we use
%
%
%   Examples
%
%   Example 1:
%
%
%   >> f = @(x) x.^2; [pp, out_param] = funappxlocal_g(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x1857 double]
%      coefs: [1856x2 double]
%     pieces: 1856
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              f: @(x)x.^2
%              a: 0
%              b: 1
%         abstol: 1.0000e-06
%            nlo: 9
%            nhi: 100
%          ninit: 30
%        npoints: 1857
%     errorbound: 7.7413e-07
%          nstar: [1x64 double]
%
%
%   Example 2:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappxlocal_g(f,-2,2,1e-7,10,20)
% 
% pp = 
% 
%       form: 'pp'
%     breaks: [1x34817 double]
%      coefs: [34816x2 double]
%     pieces: 34816
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-07
%              b: 2
%              f: @(x)x.^2
%            nhi: 20
%            nlo: 10
%          ninit: 18
%        npoints: 34817
%     errorbound: 5.9398e-08
%          nstar: [1x2048 double]
%
%
%   Example 3:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappxlocal_g(f,'a',-2,'b',2,'nhi',100,'nlo',10)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x8065 double]
%      coefs: [8064x2 double]
%     pieces: 8064
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-06
%              b: 2
%              f: @(x)x.^2
%            nhi: 100
%            nlo: 10
%          ninit: 64
%        npoints: 8065
%     errorbound: 6.3562e-07
%          nstar: [1x128 double]
%
%
%   Example 4:
%
%   >> in_param.a = -10; in_param.b = 10; f = @(x) x.^2;
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 100;
%   >> [pp, out_param] = funappxlocal_g(f,in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x45569 double]
%      coefs: [45568x2 double]
%     pieces: 45568
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -10
%         abstol: 1.0000e-06
%              b: 10
%              f: @(x)x.^2
%            nhi: 100
%            nlo: 10
%          ninit: 90
%        npoints: 45569
%     errorbound: 4.7678e-07
%          nstar: [1x512 double]
%
%
%   See also INTEGRAL_G, MEANMC_G, CUBMC_G, FUNMIN_G
%
%
%   References
%
%   [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%        Yizhi Zhang, The Cost of Deterministic, Adaptive, Automatic
%        Algorithms: Cones, Not Balls, Journal of Complexity 30 (2014), 
%        pp. 21-45.
%        
%
%   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%        Lluís Antoni Jiménez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%        "GAIL: Guaranteed Automatic Integration Library (Version 2.0)" 
%        [MATLAB Software], 2014. Available from 
%        http://code.google.com/p/gail/
%
%        If you find GAIL helpful in your work, please support us by citing
%        the above paper and software.
%

% check parameter satisfy conditions or not
[f, out_param] = funappx_g_param(varargin{:});

%%main algorithm

% tau = ceil(out_param.tauhi*(out_param.taulo/out_param.tauhi)^(1/(1+len)));
% n = ceil((tau+1)/2)+1;
% initialize number of points
n = out_param.ninit;
index = [1 n];
% initialize nstar
nstar = n - 2;
% initialize error
err = inf;
% length of interval
len = out_param.b - out_param.a;
x = out_param.a:len/(n-1):out_param.b;
y = f(x);
iter = 0;

while(max(err) >= out_param.abstol)
    iter = iter + 1;
    % Stage 1: Find the maximum error
    tmp = find(err > out_param.abstol);
    
    % Stage 2: Computer the norm
    a = index(min(tmp));
    b = index(max(tmp+1));
    n = b-a+1;
    len = x(b)-x(a);
    diff_y = diff(y(a:b));
    %approximate the weaker norm of input function
    gn = (n-1)/len*max(abs(diff_y-(y(b)-y(a))/(n-1)));
     %approximate the stronger norm of input function
    fn = (n-1)^2/len^2*max(abs(diff(diff_y)));
    
    % Stage 3: satisfy necessary condition
    if nstar(tmp)*(2*gn+fn*len/(n-1)) >= fn*len;
        % Stage 4: check for convergence
        err(tmp) = nstar(tmp)*len*gn/(4*(n-1)*(n-1-nstar(tmp)));
        if err(tmp) >= out_param.abstol;
            % Stage 5:          
            index = [index(1:tmp) a+n-1 index(tmp+1:end)+n-1];
            h = (x(b)-x(a))/(2*(n-1));
            %xnew = repmat(x(a:b-1),1,1)+repmat(h,1,n-1);
            xnew = x(a:b-1)+h;
            ynew = f(xnew);
            xnew1 = [x(a:b-1); xnew];
            xx = xnew1(:)';
            ynew1 = [y(a:b-1); ynew];
            yy = ynew1(:)';
%             xx = x(a):h:x(b);
%             yy = f(xx);
            x = [x(1:a-1) xx x(b:end)];
            y = [y(1:a-1) yy y(b:end)];
            err = [err(1:tmp-1) inf inf err(tmp+1:end)];
            ntemp=max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+h))),3);
            nstar = [nstar(1:tmp-1) ntemp-2 ntemp-2 nstar(tmp+1:end)];
        end;
    else
        % increase nstar
        nstar(tmp) = 2*nstar(tmp);
        % check if number of points large enough
        if n >= nstar(tmp) + 2;
            % true, go to Stage 4
            err(tmp) = nstar(tmp)*len*gn/(4*(n-1)*(n-1-nstar(tmp)));
            if err(tmp) >= out_param.abstol;
                % Stage 5:
                index = [index(1:tmp) a+n-1 index(tmp+1:end)+n-1];
                h = (x(b)-x(a))/(2*(n-1));
                xnew = repmat(x(a:b-1),1,1)+repmat(h,1,n-1);
                ynew = f(xnew);
                xnew = [x(a:b-1); xnew];
                xx = [xnew(:); x(b)]';
                ynew = [y(a:b-1); ynew];
                yy = [ynew(:); y(b)]';
                x = [x(1:a-1) xx x(b+1:end)];
                y = [y(1:a-1) yy y(b+1:end)];
                err = [err(1:tmp-1) inf inf err(tmp+1:end)];
                ntemp=max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+h))),3);
                nstar = [nstar(1:tmp-1) ntemp-2 ntemp-2 nstar(tmp+1:end)];
            end;
        else
            % Stage 5:
            index = [index(1:tmp) a+n-1 index(tmp+1:end)+n-1];
            h = (x(b)-x(a))/(2*(n-1));
            xnew = repmat(x(a:b-1),1,1)+repmat(h,1,n-1);
            ynew = f(xnew);
            xnew = [x(a:b-1); xnew];
            xx = [xnew(:); x(b)]';
            ynew = [y(a:b-1); ynew];
            yy = [ynew(:); y(b)]';
            x = [x(1:a-1) xx x(b+1:end)];
            y = [y(1:a-1) yy y(b+1:end)];
            err = [err(1:tmp-1) inf inf err(tmp+1:end)];
            nstar = [nstar(1:tmp-1) nstar(tmp) nstar(tmp) nstar(tmp+1:end)];
        end;
    end;    
end;
out_param.npoints = index(end);
out_param.errorbound = max(err);
out_param.nstar = nstar;
out_param.iter = iter;
% out_param.err = err;
% x1 = x;
% y1 = f(x1);
pp = interp1(x,y,'linear','pp');
%fappx = @(x) interp1(x1,y,x,'linear');


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.nlo = 9;
default.nhi = 100;

if isempty(varargin)
    warning('MATLAB:funappx_g:nofunction','Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].')
    help funappxlocal_g
    f = @(x) exp(-100*(x-0.5).^2);
    out_param.f = f;
else
    f = varargin{1};
    out_param.f = f;
end;

validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) ...
        || ischar(in2));
end

if ~validvarargin
    %if only one input f, use all the default parameters
    %warning('MATLAB:funappx_g:inputnotcorr','Input can not be recognized. Use default parameters in GAIL.')
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
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
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'nlo',default.nlo,@isnumeric);
        addParamValue(p,'nhi',default.nhi,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let end point of interval not be infinity
if (out_param.a == inf||out_param.a == -inf)
    warning(['a can not be infinity. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf)
    warning(['b can not be infinity. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('MATLAB:funappx_g:blea','b can not be smaller than a; exchange these two. ')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('MATLAB:funappx_g:beqa',['b can not equal a. Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
        ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end


% if (out_param.taulo > out_param.tauhi)
%     out_param.tauhi = out_param.taulo;
% end;
% if (~gail.isposint(out_param.taulo))
%     if gail.isposge2(out_param.taulo)
%         warning('MATLAB:funappx_g:lowtau',['Lower bound of cone condition should be a positive integer.' ...
%             ' Using ', num2str(ceil(out_param.taulo))])
%         out_param.taulo = ceil(out_param.taulo);
%     else
%         warning('MATLAB:funappx_g:lowtault2',[' Lower bound of cone condition of points should be a positive integer.' ...
%             ' Using default number of points ' int2str(default.taulo)])
%         out_param.taulo = default.taulo;
%     end
% end
% if (~gail.isposint(out_param.tauhi))
%     if gail.isposge2(out_param.tauhi)
%         warning('MATLAB:funappx_g:hitau',['Upper bound of cone condition should be a positive integer.' ...
%             ' Using ', num2str(ceil(out_param.tauhi))])
%         out_param.tauhi = ceil(out_param.tauhi);
%     else
%         warning('MATLAB:funappx_g:hitault2',[' Upper bound of cone condition should be a positive integer.' ...
%             ' Using default number of points ' int2str(default.tauhi)])
%         out_param.tauhi = default.tauhi;
%     end
% end
if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('MATLAB:funappx_g:lowinitnotint',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:funappx_g:lowinitlt3',[' Lower bound of initial number of points should be a positive integer greater than 3.' ...
            ' Using 3 as nlo'])
        out_param.nlo = 3;
    end
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('MATLAB:funappx_g:hiinitnotint',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:funappx_g:hiinitlt3',[' Upper bound of points should be a positive integer greater than 3.' ...
            ' Using default number of points ' int2str(default.nhi) ' as nhi' ])
        out_param.nhi = default.nhi;
    end
end

if (out_param.nlo > out_param.nhi)
     warning('MATLAB:funappx_g:logrhi', 'Lower bound of initial number of points is larger than upper bound of initial number of points; Use nhi as nlo')
%     temp = out_param.nlo;
%     out_param.nlo = out_param.nhi;
%     out_param.nhi = temp;
    out_param.nhi = out_param.nlo;
end;

h = out_param.b - out_param.a;
out_param.ninit = max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+h))),3);


