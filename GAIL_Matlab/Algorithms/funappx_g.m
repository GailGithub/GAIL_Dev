function [pp,out_param]=funappx_g(varargin)
%FUNAPPX_G 1-D guaranteed function recovery on a closed interval [a,b]
%
%   pp = FUNAPPX_G(f) approximates function f on the default interval
%   [0,1] by a piecewise polynomial structure pp within the guaranteed
%   absolute error tolerance of 1e-6. Input f is a function handle. The
%   statement y = f(x) should accept a vector argument x and return a
%   vector y of function values that is of the same size as x. Output pp
%   may be evaluated via PPVAL.
%
%   pp = FUNAPPX_G(f,a,b,abstol,nlo,nhi,nmax,maxiter) for a given
%   function f and the ordered input parameters that define the finite
%   interval [a,b], a guaranteed absolute error tolerance abstol, a lower
%   bound of initial number of points nlo, an upper bound of initial number
%   of points nhi, a cost budget nmax and max number of iteration maxiter.
%
%   pp =
%   FUNAPPX_G(f,'a',a,'b',b,'abstol',abstol,'nlo',nlo,'nhi',nhi,'nmax',nmax,'maxiter',maxiter)
%   recovers function f on the finite interval [a,b], given a guaranteed
%   absolute error tolerance abstol, a lower bound of initial number of
%   points nlo, an upper bound of initial number of points nhi, a cost
%   budget nmax and max number of iteration maxiter. All seven field-value
%   pairs are optional and can be supplied in different order.
%
%   pp = FUNAPPX_G(f,in_param) recovers function f on the finite
%   interval [in_param.a,in_param.b], given a guaranteed absolute error
%   tolerance in_param.abstol, a lower bound of initial number of points
%   in_param.nlo, an upper bound of initial number of points in_param.nhi,
%   a cost budget in_param.nmax and max number of iteration
%   in_param.maxiter. If a field is not specified, the default value is
%   used.
%
%     in_param.a --- left end point of interval, default value is 0
%
%     in_param.b --- right end point of interval, default value is 1
%
%     in_param.abstol --- guaranteed absolute error tolerance, default
%     value is 1e-6
%
%     in_param.nlo --- lower bound of initial number of points we used,
%     default value is 10
%
%     in_param.nhi --- upper bound of initial number of points we used,
%     default value is 1000
%
%     in_param.nmax --- when number of points hits the value, iteration
%     will stop, default value is 1e7
%
%     in_param.maxiter --- max number of interation, default value is 1000
%
%   [pp, out_param] = FUNAPPX_G(f,...) returns a piecewise polynomial
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
%     out_param.ninit --- initial number of points we use for each sub
%     interval
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
%     out_param.nmax --- when number of points hits the value, iteration
%     will stop
%
%     out_param.maxiter --- max number of interation
%
%  Guarantee
%
%  For [a,b] there exists a partition, P={[t_0,t_1], [t_1,t_2], ...,
%  [t_{L-1},t_L]}, where a=t_0 < t_1 < ... < t_L=b. If the function to be
%  approximated, f, satisfies the cone condition
%                              2 nstar    ||     f(t_l)-f(t_{l-1})||
%      ||f''||        <=  --------------  ||f'- ----------------- ||
%             \infty       t_l - t_{l-1}  ||        t_l - t_{l-1} ||\infty,
%  for each sub interval [t_{l-1},t_l], where 1 <= l <= L, then the pp
%  output by this algorithm is guaranteed to satisfy
%      ||f-ppval(pp, )||\infty <= abstol.
%
%   Examples
%
%   Example 1:
%
%
%   >> f = @(x) exp(-100*(x-sqrt(2)/2).^2); [pp, out_param] = funappx_g(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x6733 double]
%      coefs: [6732x2 double]
%     pieces: 6732
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% 
% out_param = 
% 
%              f: @(x)exp(-100*(x-1/sqrt(2)).^2)
%              a: 0
%              b: 1
%         abstol: 1.0000e-06
%            nlo: 10
%            nhi: 1000
%           nmax: 10000000
%        maxiter: 1000
%          ninit: 100
%        npoints: 6733
%     errorbound: 9.4644e-07
%          nstar: [1x68 double]
% 
% 
%   Example 2:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappx_g(f,-2,2,1e-7,10,20)
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
%        maxiter: 1000
%            nhi: 20
%            nlo: 10
%           nmax: 10000000
%          ninit: 18
%        npoints: 34817
%     errorbound: 5.9398e-08
%          nstar: [1x2048 double]
%
%
%   Example 3:
%
%   >> f = @(x) x.^2;
%   >> [pp, out_param] = funappx_g(f,'a',-2,'b',2,'nhi',20,'nlo',10)
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
%        maxiter: 1000
%            nhi: 20
%            nlo: 10
%           nmax: 10000000
%          ninit: 18
%        npoints: 8705
%     errorbound: 9.5037e-07
%          nstar: [1x512 double]
%
%
%   Example 4:
%
%   >> in_param.a = -5; in_param.b = 5; f = @(x) x.^2;
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 20;
%   >> [pp, out_param] = funappx_g(f,in_param)
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
%        maxiter: 1000
%            nhi: 20
%            nlo: 10
%           nmax: 10000000
%          ninit: 19
%        npoints: 36865
%     errorbound: 3.1274e-07
%          nstar: [1x2048 double]
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
%        Llu¨ªs Antoni Jim¨¦nez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%        "GAIL: Guaranteed Automatic Integration Library (Version 2.0)"
%        [MATLAB Software], 2014. Available from
%        http://code.google.com/p/gail/
%
%        If you find GAIL helpful in your work, please support us by citing
%        the above paper and software.
%

% check parameter satisfy conditions or not
[f, out_param] = funappx_g_param(varargin{:});
MATLABVERSION= gail.matlab_version;
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
    
%     gn(gn<eps/2)=0;
%     fn(fn<eps/2)=0;
    
    %find nstar not large enough then double it
    smallconeind = find(nstar.*(2*gn+fn.*len/(ninit-1)) <(fn.*len));
    nstar(smallconeind) = 2*nstar(smallconeind);
    
    %check necessary condition if satisfied then compute error
    %otherwise just use the error from last iteration
%     largeconeind = find(nstar.*(2*gn+fn.*len/(ninit-1)) >= (fn.*len));
%     err(largeconeind) = nstar(largeconeind).*len(largeconeind).*gn(largeconeind)./...
%         (4*(ninit-1).*(ninit-1-nstar(largeconeind)));
    err = nstar.*len.*gn./(4*(ninit-1).*(ninit-1-nstar));
    %check if error satisfy the error tolerance 
    if max(err) > abstol;
        %flag sub interval error not satisfy error tolerance 1 in whbad
        whbad = err > abstol;
        %add index for bad sub interval
        badind = find(whbad == 1);
        %flag sub interval error satisfy error tolerance 1 in whgood
        whgood = (whbad ==0);
        %add index for good sub interval
        goodind = find(whgood == 1);   
        %find # of new sub intervals need to be added at each sub
        %interval
        badcumsum = cumsum(whbad);
        %pickup # of new sub intervals at bad intervals
        cumbad = badcumsum(badind);
        %generate new index of sub intervals splitted from bad intervals
        newindex = [badind + [0 cumbad(1:end-1)]; badind + cumbad];
        newindex = newindex(:)';
        %find the length of each sub interval
        h = len/2/(ninit-1);
        %reshape x without end point to a matrix of ninit-1 by # of intervals
        reshapex =  reshape(x(1:end-1),ninit -1,(index(end) - 1)/(ninit -1));
        %generate new points newx need to be added
        newx = reshapex(:,badind) + repmat(h(badind),ninit-1,1);
        %compute value newy of newx
        newy = f(newx);
        %initialize a zero matrix of 2*(ninit-1) by # of bad sub intervals
        %to store all the points after splitting bad sub intervals
        badmatrix = zeros(2*(ninit-1),sum(whbad));
        %insert x at bad sub intervals in badmatrix as the row 1,
        %3,..., end-1
        badmatrix(1:2:end-1,:) = reshapex(:,badind);
        %insert newx at bad sub intervals in badmatrix as the row 2,
        %4,..., end
        badmatrix(2:2:end,:) = newx;
        %reshape badmatrix to the size of ninit -1 by 2*# of bad sub
        %intervals
        badmatreshape = reshape(badmatrix, ninit - 1, 2*sum(whbad));
        %initialize a matrix of ninit - 1 by # of sub intervals after
        %splitting bad sub intervals for x
        newreshapex = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        %insert all the points after splitting bad sub intervals to correct
        %column
        newreshapex(:,newindex) = badmatreshape;
        %insert all the points on good sub intervals to correct column
        newreshapex(:,goodind + badcumsum(goodind)) = reshapex(:,goodind);
        %obtain all the points in vector x
        x = [newreshapex(:)' x(end)];
        %insert y at bad sub intervals in badmatrix as the row 1,
        %3,..., end-1
        badmatrix(1:2:end-1,:) = reshapey(:,badind);
        %insert newy at bad sub intervals in badmatrix as the row 2,
        %4,..., end
        badmatrix(2:2:end,:) = newy;
        %reshape badmatrix to the size of ninit -1 by 2*# of bad sub
        %intervals
        badmatreshape = reshape(badmatrix, ninit - 1, 2*sum(whbad));
        %initialize a matrix of ninit - 1 by # of sub intervals after
        %splitting bad sub intervals for y
        newreshapey = zeros(ninit - 1, 2*sum(whbad)+sum(whgood));
        %insert all the values after splitting bad sub intervals to correct
        %column
        newreshapey(:,newindex) = badmatreshape;
        %insert all the original y on good sub intervals to correct column
        newreshapey(:,goodind + badcumsum(goodind)) = reshapey(:,goodind);
        %obtain all the values in vector y
        y = [newreshapey(:)' y(end)];
        
        %generate error for new sub intervals
        %initialize a vertor of # of sub intervals after splitting
        newerr = zeros(1,2*sum(whbad)+sum(whgood));
        %use the same error for splitted bad interval
        baderr = [err(badind); err(badind)];
        %insert error after splitting bad sub intervals to correct
        %position
        newerr(newindex)=baderr(:)';
        newerr(goodind + badcumsum(goodind)) = err(goodind);
        %obtain error for all sub intervals
        err = newerr;
        
        %upadte index w.p.t x after splitting
        %update index of the original endpoints 
        index(2:end) = index(2:end) + badcumsum*(ninit-1);
        %obtain the index of new endpoins after splitting
        %if one interval not splitted, will get the same index as in
        %previous line
        indexbeg = index(1:end-1) + whbad*(ninit-1);
        %combine two index together and emlinate duplicate indices
        indexnew = [index(1:end-1); indexbeg];
        indexnew = indexnew(:)';
        index = unique([indexnew index(end)]);
    else
        break;
    end;
    if(iter>= out_param.maxiter)
        warning(['MATLAB:funappx_g:exceediter','Iteration exceeds max iteration' num2str(out_param.maxiter)])
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
pp = interp1(x,y,'linear','pp');
if MATLABVERSION >= 8.3
    warning('on', 'MATLAB:interp1:ppGriddedInterpolant');
end;



function [f, out_param] = funappx_g_param(varargin)
% parse the input to the funappx_g function

%% Default parameter values

default.abstol = 1e-6;
default.a = 0;
default.b = 1;
default.nlo = 10;
default.nhi = 1000;
default.nmax = 1e7;
default.maxiter = 1000;

if isempty(varargin)
    warning('MATLAB:funappx_g:nofunction','Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].')
    help funappx_g
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
        addOptional(p,'nmax',default.nmax,@isnumeric)
        addOptional(p,'maxiter',default.maxiter,@isnumeric)
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
        addParamValue(p,'maxiter',default.maxiter,@isnumeric);
    end
    parse(p,f,varargin{2:end})
    out_param = p.Results;
end;

% let end point of interval not be infinity
if (out_param.a == inf||out_param.a == -inf)
    warning('MATLAB:funappx_g:aisinf',['a cannot be infinity. Use default a = ' num2str(default.a)])
    out_param.a = default.a;
end;
if (out_param.b == inf||out_param.b == -inf)
    warning(['MATLAB:funappx_g:bisinf','b cannot be infinity. Use default b = ' num2str(default.b)])
    out_param.b = default.b;
end;

if (out_param.b < out_param.a)
    warning('MATLAB:funappx_g:blea','b cannot be smaller than a; exchange these two. ')
    tmp = out_param.b;
    out_param.b = out_param.a;
    out_param.a = tmp;
elseif(out_param.b == out_param.a)
    warning('MATLAB:funappx_g:beqa',['b cannot equal a. Use b = ' num2str(out_param.a+1)])
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

if (~gail.isposint(out_param.maxiter))
    if gail.isposintive(out_param.maxiter)
        warning('MATLAB:funappx_g:maxiternotint',['Max iteration should be a positive integer.' ...
            ' Using max iteration as  ', num2str(ceil(out_param.maxiter))])
        out_param.nmax = ceil(out_param.nmax);
    else
        warning('MATLAB:funappx_g:budgetisneg',['Max iteration should be a positive integer.' ...
            ' Using default max iteration as ' int2str(default.maxiter)])
        out_param.nmax = default.nmax;
    end;
end



