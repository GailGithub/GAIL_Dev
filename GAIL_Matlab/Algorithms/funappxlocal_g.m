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
%   pp = FUNAPPXLOCAL_G(f,a,b,abstol,taulo,tauhi) for a given function f
%   and the ordered input parameters that define the finite interval [a,b],
%   a guaranteed absolute error tolerance abstol, a lower bound of initial
%   cone condition taulo, and an upper bound of initial cone condition
%   tauhi.
%
%   pp =
%   FUNAPPXLOCAL_G(f,'a',a,'b',b,'abstol',abstol,'taulo',taulo,'tauhi',tauhi)
%   recovers function f on the finite interval [a,b], given a guaranteed
%   absolute error tolerance abstol, a lower bound of initial cone
%   condition taulo, and an upper bound of initial cone condition tauhi.
%   All five field-value pairs are optional and can be supplied in
%   different order.
%
%   pp = FUNAPPXLOCAL_G(f,in_param) recovers function f on the finite
%   interval [in_param.a,in_param.b], given a guaranteed absolute error
%   tolerance in_param.abstol, a lower bound of initial cone condition
%   in_param.taulo, an upper bound of initial cone condition
%   in_param.tauhi. If a field is not specified, the default value is used.
%
%     in_param.a --- left end point of interval, default value is 0
%
%     in_param.b --- right end point of interval, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6
%
%     in_param.taulo --- lower bound of cone condition we used,
%     default value is 9
%
%     in_param.tauhi --- upper bound of cone condition we used,
%     default value is 100
%
%   [pp, out_param] = FUNAPPXLOCAL_G(f,...) returns a piecewise polynomial
%   structure pp and an output structure out_param, which has the following
%   fields:
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
%     out_param.tau --- a vector indicate the cone condition of each
%     subinterval
%
%     out_param.a --- left end point of interval
%
%     out_param.b --- right end point of interval
%
%     out_param.abstol --- guaranteed absolute error tolerance
%
%     out_param.taulo --- a lower bound of cone condtion
%
%     out_param.tauhi --- an upper bound of cone condition
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
%     breaks: [1x2049 double]
%      coefs: [2048x2 double]
%     pieces: 2048
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
%          taulo: 9
%          tauhi: 100
%        npoints: 2049
%     errorbound: 4.0640e-07
%            tau: [1x128 double]
%
%
%   Example 2:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappxlocal_g(f,-2,2,1e-7,10,10)
% 
% pp = 
% 
%       form: 'pp'
%     breaks: [1x49153 double]
%      coefs: [49152x2 double]
%     pieces: 49152
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
%          tauhi: 10
%          taulo: 10
%        npoints: 49153
%     errorbound: 4.1392e-08
%            tau: [1x8192 double]
%
%
%   Example 3:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappxlocal_g(f,'a',-2,'b',2,'tauhi',100,'taulo',10)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x8449 double]
%      coefs: [8448x2 double]
%     pieces: 8448
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
%          tauhi: 100
%          taulo: 10
%        npoints: 8449
%     errorbound: 3.5870e-07
%            tau: [1x256 double]
%
%
%   Example 4:
%
%   >> in_param.a = -10; in_param.b = 10; f = @(x) x.^2;
%   >> in_param.abstol = 10^(-7); in_param.taulo = 10; in_param.tauhi = 100;
%   >> [pp, out_param] = funappxlocal_g(f,in_param)
%
% pp = 
%       form: 'pp'
%     breaks: [1x94209 double]
%      coefs: [94208x2 double]
%     pieces: 94208
%      order: 2
%        dim: 1
%     orient: 'first'
%  
% out_param = 
% 
%              a: -10
%         abstol: 1.0000e-07
%              b: 10
%              f: @(x)x.^2
%          tauhi: 100
%          taulo: 10
%        npoints: 94209
%     errorbound: 6.8856e-08
%            tau: [1x2048 double]
%
%
%   See also INTEGRAL_G, MEANMC_G, CUBMC_G
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
%        and Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library
%        (Version 2.0)" [MATLAB Software], 2014. Available from
%        http://code.google.com/p/gail/
%
%        If you find GAIL helpful in your work, please support us by citing
%        the above paper and software.
%
[f, out_param] = funappx_g_param(varargin{:});
len = out_param.b - out_param.a;
tau = ceil(out_param.tauhi*(out_param.taulo/out_param.tauhi)^(1/(1+len)));
n = ceil((tau+1)/2)+1;
index = [1 n];
err = inf;
x = out_param.a:len/(n-1):out_param.b;
y = f(x);

while(max(err) >= out_param.abstol)
    % Stage 1: Find the maximum error
    tmp = find(err==max(err),1);
    
    % Stage 2: Computer the norm
    a = index(tmp); b=index(tmp+1);
    n = b-a+1;
    len = x(b)-x(a);
    diff_y = diff(y(a:b));
    %approximate the weaker norm of input function
    gn = (n-1)/len*max(abs(diff_y-(y(b)-y(a))/(n-1)));
     %approximate the stronger norm of input function
    fn = (n-1)^2/len^2*max(abs(diff(diff_y)));
    
    % Stage 3: satisfy necessary condition
    if tau(tmp)*(gn+fn*len/2/(n-1)) >= fn*len;
        % Stage 4: check for convergence
        err(tmp) = tau(tmp)*len*gn/(4*(n-1)*(2*n-2-tau(tmp)));
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
            tautemp = ceil(out_param.tauhi*(out_param.taulo/out_param.tauhi)^(1/(1+h)));
            tau = [tau(1:tmp-1) tautemp tautemp tau(tmp+1:end)];
        end;
    else
        % increase tau
        tau(tmp) = 2*tau(tmp);
        % check if number of points large enough
        if n > (tau(tmp)+1)/2;
            % true, go to Stage 4
            err(tmp) = tau(tmp)*len*gn/(4*(n-1)*(2*n-2-tau(tmp)));
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
                tautemp = ceil(out_param.tauhi*(out_param.taulo/out_param.tauhi)^(1/(1+h)));
                tau = [tau(1:tmp-1) tautemp tautemp tau(tmp+1:end)];
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
            tau = [tau(1:tmp-1) tau(tmp) tau(tmp) tau(tmp+1:end)];
        end;
    end;    
end;
out_param.npoints = index(end);
out_param.errorbound = max(err);
out_param.tau = tau;
% out_param.err = err;
x1 = x;
y1 = f(x1);
pp = interp1(x1,y1,'linear','pp');
%fappx = @(x) interp1(x1,y,x,'linear');


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.taulo = 9;
default.tauhi = 100;

if isempty(varargin)
    warning('MATLAB:funappx_g:nofunction','Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5).^2).')
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
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
    out_param.taulo = default.taulo;
    out_param.tauhi = default.tauhi;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'taulo',default.taulo,@isnumeric);
        addOptional(p,'tauhi',default.tauhi,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'taulo',default.taulo,@isnumeric);
        addParamValue(p,'tauhi',default.tauhi,@isnumeric);
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


if (out_param.taulo > out_param.tauhi)
    out_param.tauhi = out_param.taulo;
end;
if (~gail.isposint(out_param.taulo))
    if gail.isposge2(out_param.taulo)
        warning('MATLAB:funappx_g:lowtau',['Lower bound of cone condition should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.taulo))])
        out_param.taulo = ceil(out_param.taulo);
    else
        warning('MATLAB:funappx_g:lowtault2',[' Lower bound of cone condition of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.taulo)])
        out_param.taulo = default.taulo;
    end
end
if (~gail.isposint(out_param.tauhi))
    if gail.isposge2(out_param.tauhi)
        warning('MATLAB:funappx_g:hitau',['Upper bound of cone condition should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.tauhi))])
        out_param.tauhi = ceil(out_param.tauhi);
    else
        warning('MATLAB:funappx_g:hitault2',[' Upper bound of cone condition should be a positive integer.' ...
            ' Using default number of points ' int2str(default.tauhi)])
        out_param.tauhi = default.tauhi;
    end
end

