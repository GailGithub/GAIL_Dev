function [fmin,out_param]=funminNoPenalty_g(varargin)
%funminNoPenalty_g 1-D guaranteed locally adaptive function optimization 
%   on [a,b]
%
%   fmin = funminNoPenalty_g(f) 
%
%   fmin = funminNoPenalty_g(f,a,b,abstol) 
%
%   fmin = funminNoPenalty_g(f,'a',a,'b',b,'abstol',abstol) 
%
%   fmin = funminNoPenalty_g(f,in_param) 
%
%   [fmin, out_param] = funminNoPenalty_g(f,...) 
%
%   Properties
%    
%
%   Input Arguments
%
%
%   Optional Input Arguments
%
%
%   Output Arguments
%
%
%  Guarantee
%
%
%   Examples
%
%
%   
%   See also 
%
%
%  References
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%

% check parameter satisfy conditions or not
[f, in_param] = funminNoPenalty_g_param(varargin{:});
MATLABVERSION = gail.matlab_version;
out_param = in_param;
out_param = rmfield(out_param,'memorytest');
out_param = rmfield(out_param,'output_x');

%% main algorithm
a = out_param.a;
b = out_param.b;
abstol = out_param.abstol;
n = out_param.ninit;
out_param.x = a:(b-a)/(n-1):b;
y = f(out_param.x);
Un=min(y);
fh = 4*(b-a)/(n-1);
C0 = 1.2;
max_errest = 1;
iSing = find(isinf(y));
if ~isempty(iSing)
    error('GAIL:funminNoPenalty_g:yInf',['Function f(x) = Inf at x = ', num2str(out_param.x(iSing))]);
end
if length(y) == 1  
    % probably f is a constant function and Matlab would  
    % reutrn only a value fmin 
    fmin = y;
    max_errest = 0;
end
iter = 0;
exit_len = 2;
% we start the algorithm with all warning flags down
out_param.exit = false(1,exit_len); 
C = @(h) C0*fh./(fh-h);

while(max_errest > abstol)
    %% Stage 1: compute length of each subinterval and approximate |f''(t)|
    len = out_param.x(2:end)-out_param.x(1:end-1);
    deltaf = 2*(y(1:end-2)./len(1:end-1)./(len(1:end-1)+len(2:end))-...
                y(2:end-1)./len(1:end-1)./ len(2:end)              +...
                y(3:end  )./len(2:end  )./(len(1:end-1)+len(2:end)));
    deltaf = [0 0 abs(deltaf) 0 0];
    
    %% Stage 2: compute bound of |f''(t)| and estimate error
    h = [out_param.x(2)-a out_param.x(3)-a       ...
         out_param.x(4:end)-out_param.x(1:end-3) ...
         b-out_param.x(end-2)  b-out_param.x(end-1)];
    normbd = C(max(h(1:n-1),h(3:n+1))) .* max(deltaf(1:n-1),deltaf(4:n+2)); % appx f''
    errest = len.^2/8.*normbd;
    % update iterations
    iter = iter + 1;
    max_errest = max(errest);
    if max_errest <= abstol,
        break
    end 
 
    %% Stage 3: find I and update x,y
    diff_y=diff(y);
    ln=diff_y/2+y(1:n-1)-abs(diff_y)/2-errest;
    badinterval = (errest > abstol | ln< Un);
    whichcut = badinterval | [badinterval(2:end) 0] | [0 badinterval(1:end-1)];
    if (out_param.nmax<(n+length(find(whichcut==1))))
        out_param.exit(1) = true;
        warning('GAIL:funminNoPenalty_g:exceedbudget',['funminNoPenalty_g'...
            'attempted to exceed the cost budget. The answer may be '...
            'unreliable.'])
        break;
    end; 
    if(iter==out_param.maxiter)
        out_param.exit(2) = true;
        warning('GAIL:funminNoPenalty_g:exceediter',['Number of iterations has '...
            'reached maximum number of iterations.'])
        break;
    end;
    newx = out_param.x(whichcut) + 0.5 * len(whichcut);
    tt = cumsum(whichcut); 
    out_param.x([1 (2:n)+tt]) = out_param.x;
    y([1 (2:n)+tt]) = y;
    tem = 2 * tt + cumsum(whichcut==0);
    out_param.x(tem(whichcut)) = newx;
    y(tem(whichcut)) = f(newx);
    n = length(out_param.x);

 
end;

%% postprocessing
out_param.iter = iter;
out_param.npoints = n;
out_param.errest = max_errest;
% control the order of out_param
out_param = orderfields(out_param, ...
            {'f', 'a', 'b','abstol','nlo','nhi','ninit','nmax','maxiter',...
             'exit','iter','npoints','errest','x'});
fmin = Un;
if (in_param.memorytest)
  w = whos;
  out_param.bytes = sum([w.bytes]);
end
if (~in_param.output_x)
  out_param = rmfield(out_param,'x');
end


function [f, out_param] = funminNoPenalty_g_param(varargin)
% parse the input to the funminNoPenalty_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;
default.maxiter = 1000;
default.memorytest = false;
default.output_x = false;

MATLABVERSION = gail.matlab_version;
if MATLABVERSION >= 8.3
    f_addParamVal = @addParameter;
else
    f_addParamVal = @addParamValue;
end;

 
if isempty(varargin)
  warning('GAIL:funminNoPenalty_g:nofunction',['Function f must be specified. '...
      'Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval '...
      '[0,1].'])
  help funminNoPenalty_g
  f = @(x) exp(-100*(x-0.5).^2);
  out_param.f = f;
else
  if gail.isfcn(varargin{1})
    f = varargin{1};
    out_param.f = f;
  else
    warning('GAIL:funminNoPenalty_g:notfunction',['Function f must be a '...
        'function handle. Now GAIL is using f(x)=exp(-100*(x-0.5)^2).'])
    f = @(x) exp(-100*(x-0.5).^2);
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
    out_param.nmax = default.nmax ;
    out_param.maxiter = default.maxiter;
    out_param.memorytest = default.memorytest;
    out_param.output_x = default.output_x;
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
        addOptional(p,'nmax',default.nmax,@isnumeric)
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
        addOptional(p,'memorytest',default.memorytest,@logical)
        addOptional(p,'output_x',default.output_x,@logical)
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
        f_addParamVal(p,'memorytest',default.memorytest,@logical);
        f_addParamVal(p,'output_x',default.output_x,@logical);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let end point of interval not be infinity
if (out_param.a == inf||out_param.a == -inf)
    warning('GAIL:funminNoPenalty_g:aisinf',['a cannot be infinity. '...
        'Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf)
    warning(['GAIL:funminNoPenalty_g:bisinf','b cannot be infinity. '...
        'Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('GAIL:funminNoPenalty_g:blea',['b cannot be smaller than a;'...
        ' exchange these two. '])
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('GAIL:funminNoPenalty_g:beqa',['b cannot equal a. '...
        'Use b = ' num2str(out_param.a+1)])
    out_param.b = out_param.a+1;
end;

% let error tolerance greater than 0
if (out_param.abstol <= 0 )
    warning('GAIL:funminNoPenalty_g:tolneg', ['Error tolerance should be greater'...
        ' than 0. Using default error tolerance ' num2str(default.abstol)])
    out_param.abstol = default.abstol;
end
% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.isposintive(out_param.nmax)
        warning('GAIL:funminNoPenalty_g:budgetnotint',['Cost budget should be '...
            'a positive integer. Using cost budget '...
            , num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:funminNoPenalty_g:budgetisneg',['Cost budget should be '...
            'a positive integer. Using default cost budget '...
            int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

if (~gail.isposint(out_param.nlo))
    if gail.isposge3(out_param.nlo)
        warning('GAIL:funminglobal_g:lowinitnotint',['Lower bound of '...
        'initial number of points should be a positive integer.' ...
            ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
    else
        warning('GAIL:funminglobal_g:lowinitlt3',[' Lower bound of '...
        'initial number of points should be a positive integer greater'...
        ' than 3. Using 3 as nlo'])
        out_param.nlo = 3;
   end
        warning('GAIL:funminglobal_g:lowinitnotint',['Lower bound of '...
        'initial nstar should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.nlo)) ' as nlo '])
        out_param.nlo = ceil(out_param.nlo);
end
 if (~gail.isposint(out_param.nhi))
    if gail.isposge3(out_param.nhi)
        warning('GAIL:funminglobal_g:hiinitnotint',['Upper bound of '...
        'initial number of points should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
        out_param.nhi = ceil(out_param.nhi);
    else
        warning('GAIL:funminglobal_g:hiinitlt3',[' Upper bound of '...
        'points should be a positive integer greater than 3. Using '...
        'default number of points ' int2str(default.nhi) ' as nhi' ])
        out_param.nhi = default.nhi;
    end
         warning('GAIL:funminglobal_g:hiinitnotint',['Upper bound of '...
        'initial nstar should be a positive integer.' ...
        ' Using ', num2str(ceil(out_param.nhi)) ' as nhi' ])
        out_param.nhi = ceil(out_param.nhi);
end

if (out_param.nlo > out_param.nhi)
    warning('GAIL:funminNoPenalty_g:logrhi', ['Lower bound of initial number of'...
        ' points is larger than upper bound of initial number of '...
        'points; Use nhi as nlo'])
    out_param.nhi = out_param.nlo;
end;

h = out_param.b - out_param.a;
out_param.ninit = ceil(out_param.nhi*(out_param.nlo/out_param.nhi)...
    ^(1/(1+h)));

if (~gail.isposint(out_param.maxiter))
    if gail.ispositive(out_param.maxiter)
        warning('GAIL:funminNoPenalty_g:maxiternotint',['Max number of '...
            'iterations should be a positive integer. Using max number '...
            'of iterations as  ', num2str(ceil(out_param.maxiter))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('GAIL:funminNoPenalty_g:budgetisneg',['Max number of iterations'...
            ' should be a positive integer. Using max number of '...
            'iterations as ' int2str(default.maxiter)])
        out_param.nmax = default.nmax;
    end;
end
if (out_param.memorytest~=true&&out_param.memorytest~=false)
    warning('GAIL:funminNoPenalty_g:memorytest', ['Input of memorytest'...
        ' can only be true or false; use default value false'])
    out_param.memorytest = false;
end;
if (out_param.output_x~=true&&out_param.output_x~=false)
    warning('GAIL:funminNoPenalty_g:output_x', ['Input of output_x'...
        ' can only be true or false; use default value false'])
    out_param.output_x = false;
end;
