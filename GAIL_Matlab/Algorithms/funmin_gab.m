function [fmin,out_param]=funmin_gab(varargin)
%   See also FUNAPPX_G, INTEGRAL_G
%
%  References
%   [1]  Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, 
%   and Yizhi Zhang. The Cost of Deterministic, Adaptive, Automatic 
%   Algorithms: Cones, Not Balls. Journal of Complexity, 30:21-45, 2014
%
%   [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
%   and Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library
%   (Version 1.3.0)" [MATLAB Software], 2014. Available from
%   http://code.google.com/p/gail/
%
%   If you find GAIL helpful in your work, please support us by citing
%   the above paper and software.
%


% Parse and check the validity of input parameters
[f,out_param] = funmin_g_param(varargin{:});


%% Main algorithm

% initialize number of points
n = out_param.ninit;
% initialize tau
out_param.tau = 2*(n-1)-1;
% cost budget flag
out_param.exceedbudget = 1;
% tau change flag
tauchange = 0;
% length of interval
len = out_param.b-out_param.a;
% % add flag
% flag = 0;

while n < out_param.nmax;
    % Stage 1: estimate weaker and stronger norm
    x = out_param.a:len/(n-1):out_param.b;
    y = f(x);
    diff_y = diff(y);
    %approximate the weaker norm of input function
    gn = (n-1)*max(abs(diff_y-(y(n)-y(1))/(n-1)))/len;
    %approximate the stronger norm of input function
    fn = (n-1)^2*max(abs(diff(diff_y)))/len^2;
    
    % Stage 2: satisfy necessary condition of cone
    if out_param.tau*(gn/len+fn/(2*n-2))>= fn;
        % Stage 3: check for convergence
        bn = 2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
        cn = 2*(n-1)^2*abs(diff_y)./bn/len;
        Cn = min(cn,1); % check the conditions for each interval
        ln = (diff_y/2+y(1:n-1))-abs(diff_y).*(Cn+1./Cn)/4;
        % minimum values of each interval
        Ln = min(ln); % lower bound
        Un = min(y); % upper bound
        error = Un-Ln;
        % find the intervals containing minimum points
        index = find(cn<1 & ln < Un);
        m = size(index,2);
        if m > 0
            delta = (n-1)^2/len^2*diff_y(index).^2-2*bn/len*(diff_y(index)./2 ...
                +y(index)-bn*len/8/(n-1)^2-Un);
            ints = zeros(2,m);
            ints(1,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
                -len*sqrt(delta)./bn;
            ints(2,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
                +len*sqrt(delta)./bn;
            leftint = find([1 diff(index)~=1]);
            rightint = find([diff(index)~=1 1]);
            q = size(leftint,2);
            interval = zeros(2,q);
            interval(1,:) = ints(1,leftint);
            interval(2,:) = ints(2,rightint);
        else
            interval = zeros(2,0);
        end
        volumeX = sum(interval(2,:)-interval(1,:));
        % satisfy convergence
        if error < out_param.abstol || volumeX < out_param.TolX
            out_param.exceedbudget = 0; break;
        end
        % otherwise increase points number
        l = n;
        n = 2*(n-1)+1;
        
        % Stage 2: do not satisfy necessary condition
    else
        % increase tau
        out_param.tau = 2*fn/(gn/len+fn/(2*n-2));
        % change tau change flag
        tauchange = 1;
        % check if number of points large enough
        if n >= ((out_param.tau+1)/2);
            % large enough, go to Stage 3
            bn = 2*(n-1)*out_param.tau/(2*(n-1)-out_param.tau)*gn;
            cn = 2*(n-1)^2*abs(diff_y)./bn/len;
            Cn = min(cn,1); % check the conditions for each interval
            ln = (diff_y/2+y(1:n-1))-abs(diff_y).*(Cn+1./Cn)/4;
            % minimum values of each interval
            Ln = min(ln); % lower bound
            Un = min(y); % upper bound
            error = Un-Ln;
            % find the intervals containing minimum points
            index = find(cn<1 & ln < Un);
            m = size(index,2);
            if m > 0
                delta = (n-1)^2/len^2*diff_y(index).^2-2*bn/len*(diff_y(index)./2 ...
                    +y(index)-bn*len/8/(n-1)^2-Un);
                ints = zeros(2,m);
                ints(1,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
                    -len*sqrt(delta)./bn;
                ints(2,:)=x(index)+len/2/(n-1)-(n-1)*diff_y(index)./bn ...
                    +len*sqrt(delta)./bn;
                leftint = find([1 diff(index)~=1]);
                rightint = find([diff(index)~=1 1]);
                q = size(leftint,2);
                interval = zeros(2,q);
                interval(1,:) = ints(1,leftint);
                interval(2,:) = ints(2,rightint);
            else
                interval = zeros(2,0);
            end
            volumeX = sum(interval(2,:)-interval(1,:));
            % satisfy convergence
            if error < out_param.abstol || volumeX < out_param.TolX
                out_param.exceedbudget = 0; break;
            end
            % otherwise increase points number
            l = n;
            n = 2*(n-1)+1;
            
        else
            % not large enough, increase points number, and go to Stage 1
            l = n;
            n = 1 + (n-1)*ceil(out_param.tau+1/(2*n-2));
        end;
    end;
end;


% check tau change flag
if tauchange == 1
    warning('MATLAB:funmin_g:peaky','This function is peaky relative to ninit. You may wish to increase ninit for similar functions.')
end;

% check cost budget flag
if out_param.exceedbudget == 1
    n = l;
    warning('MATLAB:funmin_g:exceedbudget','funmin_g attempted to exceed the cost budget. The answer may be unreliable.')
end

fmin = Un;
out_param.npoints = n;
out_param.error = error;
out_param.volumeX = volumeX;
out_param.tauchange = tauchange;
out_param.intervals = interval;


function [f, out_param] = funmin_g_param(varargin)
% Parse the input to the funmin_g.m function

%% Default parameter values
default.a = 0;
default.b = 1;
default.abstol = 1e-6;
default.TolX = 1e-3;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;
    
if isempty(varargin)
    help funmin_g
    warning('Function f must be specified. Now funmin_g will use f(x)=(x-0.3)^2+1.')
    f = @(x) (x-0.3).^2+1;
else
    f = varargin{1};
end
     
validvarargin=numel(varargin)>1;
if validvarargin
    in2=varargin{2};
    validvarargin=(isnumeric(in2) || isstruct(in2) || ischar(in2));
end
        
if ~validvarargin
% There is only one input f or the second input is not satisfied our type.
% The default parameters are used.
    out_param.a = default.a;
    out_param.b = default.b;
    out_param.abstol = default.abstol;
    out_param.TolX = default.TolX;
    out_param.nlo = default.nlo;
    out_param.nhi = default.nhi;
    out_param.nmax = default.nmax;
else
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    if isnumeric(in2) % more inputs of numerical type. Put them in order.
        addOptional(p,'a',default.a,@isnumeric);
        addOptional(p,'b',default.b,@isnumeric);
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'TolX',default.TolX,@isnumeric);
        addOptional(p,'nlo',default.nlo,@isnumeric);
        addOptional(p,'nhi',default.nhi,@isnumeric);
        addOptional(p,'nmax',default.nmax,@isnumeric);
    else
        if isstruct(in2) % second input is a structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        addParamValue(p,'a',default.a,@isnumeric);
        addParamValue(p,'b',default.b,@isnumeric);
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'TolX',default.TolX,@isnumeric);
        addParamValue(p,'nlo',default.nlo,@isnumeric);
        addParamValue(p,'nhi',default.nhi,@isnumeric);
        addParamValue(p,'nmax',default.nmax,@isnumeric);
    end
    parse(p,f,varargin{2:end});
    out_param = p.Results;
end

% a and b can't be infinity
if (out_param.a == inf||out_param.a == -inf||isnan(out_param.a)==1)
    warning('MATLAB:funmin_g:anoinfinity',['a can not be infinity. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf||isnan(out_param.b)==1)
    warning('MATLAB:funmin_g:bnoinfinity',['b can not be infinity. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

% b is greater than a
if (out_param.b < out_param.a)
    warning('MATLAB:funmin_g:blea','b can not be smaller than a; exchange these two. ')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('MATLAB:funmin_g:beqa',['b can not equal to a. Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% Check whether the error tolerance is nonnegative
if out_param.abstol < 0
    warning(['MATLAB:funmin_g:abstolnegat ','Error tolerance should be greater than or equal to 0.' ...
        ' Using default error tolerance ', num2str(default.abstol)])
    out_param.abstol = default.abstol;
end

% Check whether the length tolerance is nonnegative
if out_param.TolX < 0
    warning(['MATLAB:funmin_g:Xtolernegat ','X tolerance should be greater than or equal to 0.' ...
        ' Using default X tolerance ' num2str(default.TolX)]);
    out_param.abstol = default.TolX;
end

% Check whether the cost budget is a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.ispositive(out_param.nmax)
        warning('MATLAB:funmin_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))]);
        out_param.nmax = ceil(out_param.nmax);
    else
        warning(['Cost budget should be a positive integer.'...
        ' funmin_g will use the default budget' num2str(out_param.nmax)]);
        out_param.nmax = default.nmax;
    end
end

% Check nlo and nhi
if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('MATLAB:funmin_g:lowinitnotint',['Lower bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('MATLAB:funmin_g:lowinitlt3',[' Lower bound of initial number of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nlo)])
        out_param.nlo = default.nlo;
    end
end
if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('MATLAB:funmin_g:hiinitnotint',['Upper bound of initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nhi))])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('MATLAB:funmin_g:hiinitlt3',[' Upper bound of points should be a positive integer.' ...
            ' Using default number of points ' int2str(default.nhi)])
        out_param.nhi = default.nhi;
    end
end
if (out_param.nlo > out_param.nhi)
    out_param.nhi = out_param.nlo;
end;
if (out_param.nlo > out_param.nmax)
    warning('MATLAB:funmin_g:logecost',['Lower bound of initial number of points should be smaller than cost budget.' ...
            ' Using ', num2str(ceil(out_param.nmax/2))])
    out_param.nlo = ceil(out_param.nmax/2);
    out_param.nhi = out_param.nlo;
end;

out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)^(1/(1+out_param.b - out_param.a)));

end

end
