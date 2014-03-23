function [fappx,out_param]=funappxlocal_g(varargin)
[f, out_param] = funappx_g_param(varargin{:});
len = out_param.b - out_param.a;
tau = ceil(out_param.tauhi*(out_param.taulo/out_param.tauhi)^(1/len));
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
            tau = [tau(1:tmp-1) tau(tmp) tau(tmp) tau(tmp+1:end)];
        end;
    else
        % increase tau
        tau(tmp) = 2*tau(tmp);
        % check if number of points large enough
        if n > (tau+1)/2;
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
                tau = [tau(1:tmp-1) tau(tmp) tau(tmp) tau(tmp+1:end)];
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
out_param.x = x;
out_param.tau = tau;
% out_param.err = err;
x1 = x;
fappx = @(x) interp1(x1,y,x,'linear');


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.taulo = 9;
default.tauhi = 100;

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
    out_param.taulo = default.taulo;
    out_param.tauhi = default.tauhi;
else
    p = inputParser;
    addRequired(p,'f',@isfcn);
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

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning(['Error tolerance should be greater than 0.' ...
        ' Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end


if (out_param.taulo > out_param.tauhi)
    out_param.tauhi = out_param.taulo;
end;
if (~isposint(out_param.taulo))
    if isposge2(out_param.taulo)
        warning('MATLAB:funappx_g:lowtau',['Lower bound of cone condition should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.taulo))])
        out_param.taulo = ceil(out_param.taulo);
    else
        warning('MATLAB:funappx_g:lowtault2',[' Lower bound of cone condition of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.taulo)])
        out_param.taulo = default.nlo;
    end
end
if (~isposint(out_param.tauhi))
    if isposge2(out_param.tauhi)
        warning('MATLAB:funappx_g:hitau',['Upper bound of cone condition should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.tauhi))])
        out_param.tauhi = ceil(out_param.tauhi);
    else
        warning('MATLAB:funappx_g:hitault2',[' Upper bound of cone condition should be a positive integer.' ...
            ' Using default number of points ' int2str(default.tauhi)])
        out_param.tauhi = default.tauhi;
    end
end

