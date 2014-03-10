function [fappx,out_param]=funappxlocal_g(varargin)
[f, out_param] = funappx_g_param(varargin{:});
n = out_param.ninit;
nstar = n - 2;
index = [1 n];
len = out_param.b - out_param.a;
x = out_param.a:len/(n-1):out_param.b;
y = f(x);
diff_y = diff(y);
%approximate the weaker norm of input function
fn = (n-1)^2/len^2*max(abs(diff(diff_y)));
err = fn*len^2/(8*(n-1)^2);

while(max(err) >= out_param.abstol)
    tmp = find(err==max(err),1);
    a = index(tmp); b=index(tmp+1);
    n = b-a+1;
    index = [index(1:tmp) a+n-1 index(tmp+1:end)+n-1];
    h = (x(b)-x(a))/(2*(n-1));
    xnew = repmat(x(a:b-1),1,1)+repmat(h,1,n-1);
    ynew = f(xnew);
    xnew = [x(a:b-1); xnew];
    xx = [xnew(:); x(b)]';
    ynew = [y(a:b-1); ynew];
    yy = [ynew(:); y(b)]';
    diff_y = diff(yy);
    fn = max(abs(diff(diff_y(1:n-1))))/(h^2);
    err1 = fn*h^2/8;
    fn = max(abs(diff(diff_y(n:2*(n-1)))))/(h^2);
    err2 =fn*h^2/8;
    err = [err(1:tmp-1) err1 err2 err(tmp+1:end)];
    x = [x(1:a-1) xx x(b+1:end)];
    y = [y(1:a-1) yy y(b+1:end)];
end;
out_param.npoints = index(end);
out_param.errorbound = max(err);
% out_param.err = err;
x1 = x;
fappx = @(x) interp1(x1,y,x,'linear');


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.ninit = 52;


if isempty(varargin)
    help funappxablocal_g
    warning('Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5).^2).')
    f = @(x) exp(-100*(x-0.5).^2);
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
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
    out_param.minlen = default.minlen;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
    if isnumeric(in2)%if there are multiple inputs with
        %only numeric, they should be put in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'ninit',default.ninit,@isnumeric);
    else
        if isstruct(in2) %parse input structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'ninit',default.ninit,@isnumeric);
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

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
        ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let initial number of points be a positive integer

if (~isposint(out_param.ninit))
    if isposge3(out_param.ninit)
        warning('MATLAB:funappx_g:lowinitnotint',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.ninit))])
        out_param.ninit = ceil(out_param.ninit);
    else
        warning('MATLAB:funappx_g:lowinitlt3',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.ninit)])
        out_param.ninit = default.ninit;
    end
end

