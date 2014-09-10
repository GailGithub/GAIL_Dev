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
%   pp = FUNAPPXLOCAL_G(f,a,b,abstol,nlo,nhi,nmax) for a given function f and
%   the ordered input parameters that define the finite interval [a,b], a
%   guaranteed absolute error tolerance abstol, a lower bound of initial
%   number of points nlo, an upper bound of initial number of points nhi,
%   and a cost budget nmax.
%
%   pp = FUNAPPXLOCAL_G(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax)
%   recovers function f on the finite interval [a,b], given a guaranteed
%   absolute error tolerance abstol, a lower bound of initial number of
%   points nlo, an upper bound of initial number of points nhi, and a cost
%   budget nmax. All six field-value pairs are optional and can be supplied
%   in different order.
%
%   pp = FUNAPPXLOCAL_G(f,in_param) recovers function f on the finite
%   interval [in_param.a,in_param.b], given a guaranteed absolute error
%   tolerance in_param.abstol, a lower bound of initial number of points
%   in_param.nlo, an upper bound of initial number of points in_param.nhi,
%   and a cost budget in_param.nmax. If a field is not specified, the
%   default value is used.
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
%     in_param.nmax --- cost budget, default value is 1e7
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
%     out_param.exceedbudget --- it is 0 if the number of points used in 
%     the construction of pp is less than cost budget, 1 otherwise.
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
%     out_param.nmax --- cost budget
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
%           iter: 7
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
%           iter: 12
%
%
%   Example 3:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappxlocal_g(f,'a',-2,'b',2,'nhi',20,'nlo',10)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x8705 double]
%      coefs: [8704x2 double]
%     pieces: 8704
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
%            nhi: 20
%            nlo: 10
%          ninit: 18
%        npoints: 8705
%     errorbound: 9.5037e-07
%          nstar: [1x512 double]
%           iter: 10
%
%
%   Example 4:
%
%   >> in_param.a = -5; in_param.b = 5; f = @(x) x.^2;
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 20;
%   >> [pp, out_param] = funappxlocal_g(f,in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x36865 double]
%      coefs: [36864x2 double]
%     pieces: 36864
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -5
%         abstol: 1.0000e-06
%              b: 5
%              f: @(x)x.^2
%            nhi: 20
%            nlo: 10
%          ninit: 19
%        npoints: 36865
%     errorbound: 3.1274e-07
%          nstar: [1x2048 double]
%           iter: 12
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
[~,~,~,MATLABVERSION]=GAILstart(0);
if MATLABVERSION >= 8.3
    warning('off', 'MATLAB:interp1:ppGriddedInterpolant');
end;

%%main algorithm
% initialize number of points
ninit = out_param.ninit;
index = [1 ninit];
% initialize nstar
nstar = ninit - 2;
% initialize error
abstol = out_param.abstol;
err = abstol + 1;
len = out_param.b - out_param.a;
x = out_param.a:len/(ninit-1):out_param.b;
y = f(x);
iter = 0;

while(max(err) > abstol)
    iter = iter + 1;
    % length of each subinterval
    len = x(index(2:end))-x(index(1:end-1));
    reshapey = reshape(y(1:end-1),ninit - 1, (index(end) - 1)/(ninit -1));
    diffy = diff([reshapey;y(index(2:end))]);
    %approximate the weaker norm of input function at different subinterval
    gn = (ninit-1)./len.*max(abs(diffy-repmat((y(index(2:end))-y(index(1:end-1)))/(ninit-1),ninit-1,1)),[],1);
    %approximate the stronger norm of input function at different subinterval
    fn = (ninit-1)^2./(len.^2).*max(abs(diff(diffy)),[],1);
    %update cone condition every iteration
    ntemp=max(ceil(out_param.nhi*(out_param.nlo/out_param.nhi).^(1./(1+len))),3);
    nstar = ntemp -2;
    
    if nstar.*(2*gn+fn.*len/(ninit-1)) >= (fn.*len);
        err = nstar.*len.*gn./(4*(ninit-1).*(ninit-1-nstar));
%     else
%         smallconeindex = find(nstar*(2*gn+fn.*len/(ninit-1)) < (fn.*len));
%         err(smllconeindex) = abstol + 1;
    end
    if max(err) > abstol;
        whbad = err > abstol;
        whgood = (whbad ==0);
        %update x and y
        reshapex =  reshape(x(1:end-1),ninit -1,(index(end) - 1)/(ninit -1));
        h = len/2/(ninit-1);
        badind = find(whbad == 1);
        goodind = find(whgood == 1);
        temp = cumsum(whbad);
        cumbad = temp(badind);
        newindex = [badind + [0 cumbad(1:end-1)]; badind + cumbad];
        newindex = newindex(:)';
        newx = reshapex(:,badind) + repmat(h(badind),ninit-1,1);
        newy = f(newx);
        ll = zeros(2*(ninit-1),sum(whbad));
        ll(1:2:end-1,:) = reshapex(:,badind);
        ll(2:2:end,:) = newx;
        llreshapex = reshape(ll, ninit - 1, 2*sum(whbad));
        newreshapex = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        newreshapex(:,newindex) = llreshapex;
        newreshapex(:,goodind + temp(goodind)) = reshapex(:,goodind);
        x = [newreshapex(:)' x(end)];
        ll(1:2:end-1,:) = reshapey(:,badind);
        ll(2:2:end,:) = newy;
        llreshapey = reshape(ll, ninit - 1, 2*sum(whbad));
        newreshapey = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        newreshapey(:,newindex) = llreshapey;
        newreshapey(:,goodind + temp(goodind)) = reshapey(:,goodind);
        y = [newreshapey(:)' y(end)];
        %update err
        cumwhbad = cumsum(whbad+1);
        c = zeros(1,cumwhbad(end));
        newerr = [err(badind); err(badind)];
        c(newindex)=newerr(:)';
        c(goodind + temp(goodind)) = err(goodind);
        err = c;
        %upadte index
        index(2:end) = index(2:end) + cumsum(whbad)*(ninit-1);
        indexbeg = index(1:end-1) + whbad*(ninit-1);
        indexnew = [index(1:end-1); indexbeg];
        indexnew = indexnew(:)';
        index = unique([indexnew index(end)]);
    else
        break;
    end;
    if(iter>=1000)
        warning('MATLAB:funappx_g:exceediter','Iteration exceeds 1000 times.')
        break;
    end;
    if(index(end) >= out_param.nmax)
        warning('MATLAB:funappx_g:exceedbudget','funappx_g attempted to exceed the cost budget. The answer may be unreliable.')
        break;
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
if MATLABVERSION >= 8.3
    warning('on', 'MATLAB:interp1:ppGriddedInterpolant');
end;

%fappx = @(x) interp1(x1,y,x,'linear');


function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.nlo = 9;
default.nhi = 100;
default.nmax = 1e7;

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
    out_param.nmax = default.nmax ;
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
        addParamValue(p,'nmax',default.nmax,@isnumeric);
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
% let cost budget be a positive integer
if (~gail.isposint(out_param.nmax))
    if gail.isposintive(out_param.nmax)
        warning('MATLAB:funappx_g:budgetnotint',['Cost budget should be a positive integer.' ...
            ' Using cost budget ', num2str(ceil(out_param.nmax))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:funappx_g:budgetisneg',['Cost budget should be a positive integer.' ...
            ' Using default cost budget ' int2str(default.nmax)])
        out_param.nmax = default.nmax;
    end;
end

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




