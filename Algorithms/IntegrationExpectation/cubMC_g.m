function [Q,out_param] = cubMC_g(varargin)
%CUBMC_G Monte Carlo method to evaluate a multidimensional integral
%
%   [Q,out_param] = CUBMC_G(f,hyperbox) estimates the integral of f over
%   hyperbox to within a specified generalized error tolerance, tolfun =
%   max(abstol, reltol*| I |), i.e., | I - Q | <= tolfun with probability
%   at least (1 - alpha), where abstol is the absolute error tolerance, and
%   reltol is the relative error tolerance. Usually the reltol determines
%   the accuracy of the estimation, however, if | I | is rather small,
%   then abstol determines the accuracy of the estimation. Input f is a
%   function handle that accepts an n x d matrix input, where d is the
%   dimension of the hyperbox, and n is the number of points being
%   evaluated simultaneously. 
%
%   When measure is 'uniform', 'uniform box', 'normal' or 'Gaussian', the
%   input hyperbox is a 2 x d matrix, where the first row corresponds to
%   the lower limits and the second row corresponds to the upper limits.
%   When measure is 'uniform ball' or 'uniform sphere', the input hyperbox
%   is a vector with d+1 elements, where the first d values correspond to
%   the center of the ball and the last value corresponds to the radius of
%   the ball. For these last two measures, a user can optionally specify
%   what transformation should be used in order to get a uniform
%   distribution on a ball of sphere. When measure is 'uniform ball_box',
%   the box-to-ball transformation, which gets a set of points uniformly
%   distributed on a ball from a set of points uniformly distributed on a
%   box, will be used. When measure is 'uniform ball_normal', the
%   normal-to-ball transformation, which gets a set of points uniformly
%   distributed on a ball from a set of points normally distributed on the
%   space, will be used. Similarly, the measures 'uniform sphere_box' and
%   'uniform sphere_normal' can be defined. The default transformations are
%   the box-to-ball and the box-to-sphere transformations, depending on the
%   region of integration.
%
%   Q = CUBMC_G(f,hyperbox,measure,abstol,reltol,alpha)
%   estimates the integral of function f over hyperbox to within a
%   specified generalized error tolerance tolfun with guaranteed confidence
%   level 1-alpha using all ordered parsing inputs f, hyperbox, measure,
%   abstol, reltol, alpha, fudge, nSig, n1, tbudget, nbudget, flag. The
%   input f and hyperbox are required and others are optional.
%
%   Q = CUBMC_G(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol,'alpha',alpha)
%   estimates the integral of f over hyperbox to within a specified
%   generalized error tolerance tolfun with guaranteed confidence level
%   1-alpha. All the field-value pairs are optional and can be supplied in
%   different order. If an input is not specified, the default value is used.
%
%   [Q out_param] = CUBMC_G(f,hyperbox,in_param) estimates the integral of
%   f over hyperbox to within a specified generalized error tolerance
%   tolfun with the given parameters in_param and produce output parameters
%   out_param and the integral Q.
%
%   Input Arguments
%
%     f --- the integrand.
%
%     hyperbox --- the integration hyperbox. The default value is
%     [zeros(1,d); ones(1,d)], the default d is 1.
%
%     in_param.measure --- the measure for generating the random variable,
%     the default is 'uniform'. The other measures could be handled are
%     'uniform box', 'normal'/'Gaussian', 'uniform ball'/'uniform
%     ball_box'/'uniform ball_normal' and 'uniform sphere'/'uniform
%     sphere_box'/'uniform sphere_normal'. The input should be
%     a string type, hence with quotes.
%
%     in_param.abstol --- the absolute error tolerance, the default value
%     is 1e-2.
%
%     in_param.reltol --- the relative error tolerance, the default value
%     is 1e-1.
%
%     in_param.alpha --- the uncertainty, the default value is 1%.
%
%   Optional Input Arguments
%
%     in_param.fudge --- the standard deviation inflation factor, the
%     default value is 1.2.
%
%     in_param.nSig --- initial sample size for estimating the sample
%     variance, which should be a moderately large integer at least 30, the
%     default value is 1e4.
%
%     in_param.n1 --- initial sample size for estimating the sample mean,
%     which should be a moderately large positive integer at least 30, the
%     default value is 1e4.
%
%     in_param.tbudget --- the time budget to do the estimation, the
%     default value is 100 seconds.
%
%     in_param.nbudget --- the sample budget to do the estimation, the
%     default value is 1e9.
%
%     in_param.flag --- the value corresponds to parameter checking status.
%
%                         0   not checked
%
%                         1   checked by meanMC_g
%
%                         2   checked by cubMC_g
%
%   Output Arguments
%
%     Q --- the estimated value of the integral.
%
%     out_param.n --- the sample size used in each iteration.
%
%     out_param.ntot --- total sample used, including the sample used to
%     convert time budget to sample budget and the sample in each iteration
%     step.
%
%     out_param.nremain --- the remaining sample budget to estimate I. It was
%     calculated by the sample left and time left.
%
%     out_param.tau --- the total number of iterations.
%
%     out_param.hmu --- estimated integral in each iteration.
%
%     out_param.tol --- a reliable upper bound on error for each iteration.
%
%     out_param.kurtmax --- the upper bound on modified kurtosis.
%
%     out_param.time --- the time elapsed in seconds.
%
%     out_param.var --- the sample variance.
%
%     out_param.exit --- the state of program when exiting.
%
%                       0   success
%
%                       1   Not enough samples to estimate the mean
%
%                       10  hyperbox does not contain numbers
%
%                       11  hyperbox is not 2 x d when measure is 'uniform'
%                           or 'normal'
%
%                       12  hyperbox is only a point in one direction when
%                           measure is 'uniform' or 'normal'
%
%                       13  hyperbox is infinite when measure is 'uniform'
%
%                       14  hyperbox is not doubly infinite when measure
%                           is 'normal'
%
%                       15  hyperbox has an infinite coordinate for the
%                           center of the ball or sphere or a infinite radius
%                           for the ball or sphere
%
%                       16  The radius of the ball or sphere is a non-positive
%                           real number
%
%                       17  The dimension of the ball is zero
%
%                       18  Hyperbox not 1 x (d+1) when measure is 'uniform
%                           ball' or 'uniform sphere'
%
%                       19  The dimension of the sphere is smaller than 2
%
%  Guarantee
% This algorithm attempts to calculate the integral of function f over a
% hyperbox to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level 1-alpha. If the algorithm terminates
% without showing any warning messages and provides an answer Q, then the
% following inequality would be satisfied:
%
% Pr(| Q - I | <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% a function in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And
% the following inequality holds:
%
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
%  Examples
%
% Example 1:
%
% If no parameters are parsed, help text will show up as follows:
% >> cubMC_g
% ***Monte Carlo method to estimate***
%
%
% Example 2:
%
% Estimate the integral with integrand f(x) = sin(x) over the interval
% [1;2] with default parameters.
%
% >> f = @(x) sin(x); interval = [1;2];
% >> Q = cubMC_g(f,interval,'uniform',1e-3,1e-2);
% >> exactsol = 0.9564;
% >> check = double(abs(exactsol-Q) < max(1e-3,1e-2*abs(exactsol)))
% check = 1
%
%
% Example 3:
%
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) over the
% hyperbox [0 0;1 1], where x = [x1 x2] is a vector.
%
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [0 0;1 1];
% >> Q = cubMC_g(f,hyperbox,'uniform',1e-3,0);
% >> exactsol = 0.5577;
% >> check = double(abs(exactsol-Q) < 1e-3)
% check = 1
%
%
% Example 4:
%
% Estimate the integral with integrand f(x) = 2^d*prod(x1*x2*...*xd)+0.555
% over the hyperbox [zeros(1,d);ones(1,d)], where x = [x1 x2... xd] is a
% vector.
%
% >> d = 3; f = @(x) 2^d*prod(x,2)+0.555; hyperbox =[zeros(1,d); ones(1,d)];
% >> in_param.abstol = 1e-3; in_param.reltol = 1e-3;
% >> Q = cubMC_g(f,hyperbox,in_param);
% >> exactsol = 1.555;
% >> check = double(abs(exactsol-Q) < max(1e-3,1e-3*abs(exactsol)))
% check = 1
%
%
% Example 5:
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [-inf -inf;inf inf], where x = [x1 x2] is a vector.
%
% >> f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-inf -inf;inf inf];
% >> Q = cubMC_g(f,hyperbox,'normal',0,1e-2);
% >> exactsol = 1/3;
% >> check = double(abs(exactsol-Q) < max(0,1e-2*abs(exactsol)))
% check = 1
%
%
% Example 6:
% Estimate the integral with integrand f(x) = x1^2+x2^2 in the disk with
% center (0,0) and radius 1, where x = [x1 x2] is a vector.
%
% >> f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1];
% >> Q = cubMC_g(f,hyperbox,'uniform ball','abstol',1e-3,'reltol',1e-3);
% >> exactsol = pi/2;
% >> check = double(abs(exactsol-Q) < max(1e-3,1e-3*abs(exactsol)))
% check = 1
%
%
%   See also FUNAPPX_G, INTEGRAL_G, MEANMC_G, CUBLATTICE_G, CUBSOBOL_G, CUBBAYESLATTICE_G
%
%  References
%
%   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, "Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo
%   sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
%   Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), pp. 105-128,
%   Springer-Verlag, Berlin, 2014. DOI: 10.1007/978-3-642-41095-6_5
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
%   Research via Supportable Scientific Software," Journal of Open Research
%   Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
%   [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
%   Mathematical Software" [Course Slides], Illinois Institute of
%   Technology, Chicago, IL, 2013. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
%   Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
%   James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
%   Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
%   Workshop On Sustainable Software for Science: Practice and Experiences
%   (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
%   pp. 1-21, 2014.
%
%   [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
%   Statistics. London, UK: CHAPMAN & HALL
%
%   [7] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
%   Means of Random Variables, PhD Thesis, Illinois Institute of
%   Technology, 2016.
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%
%   Author:  Lan Jiang

tstart=tic;
[f,hyperbox,out_param] = cubMC_g_param(varargin{:});%check validity of inputs

%------------------------------------------------------------------------------
% TRANSFORMATION

%changing the integrand the hyperbox when measure is uniform ball or
%sphere by applying the appropriate transformation
if strcmpi(out_param.measure,'uniform ball') || strcmpi(out_param.measure,'uniform sphere')% using uniformly distributed samples on a ball or sphere
    if strcmpi(out_param.measure,'uniform ball')% using the formula for the volume of a ball
        volume = ((2.0*pi^(out_param.dim/2.0))/(out_param.dim*gamma(out_param.dim/2.0)))*out_param.radius^out_param.dim; %volume of a d-dimentional ball
    else % using the formula for the volume of a sphere
        volume = ((2.0*pi^(out_param.dim/2.0))/(gamma(out_param.dim/2.0)))*out_param.radius^(out_param.dim - 1); %volume of a d-dimentional sphere
    end

    if out_param.transf == 1 % box-to-ball or box-to-sphere transformation should be used
        if out_param.dim == 1 % It is not necessary to multiply the function f by the volume, since no transformation is being made
            hyperbox = [hyperbox - out_param.radius; hyperbox + out_param.radius];% for one dimension, the ball is actually an interval
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        else
            if strcmpi(out_param.measure,'uniform ball') % box-to-ball transformation
                f = @(t) f(gail.domain_balls_spheres.ball_psi_1(t, out_param.dim, out_param.radius, hyperbox))*volume;% the psi function computes the transformation
            else % box-to-sphere transformation
                f = @(t) f(gail.domain_balls_spheres.sphere_psi_1(t, out_param.dim, out_param.radius, hyperbox))*volume;% the psi function is the transformation
                out_param.dim = out_param.dim - 1;% the box-to-sphere transformation takes points from a (d-1)-dimensional box to a d-dimensional sphere
            end
            hyperbox = [zeros(1, out_param.dim); ones(1, out_param.dim)];% the hyperbox must be the domain of the transformation, which is a unit box
            out_param.measure = 'uniform';% then a uniform distribution on a box can be used
        end
    else % normal-to-ball or normal-to-sphere transformation should be used
        if strcmpi(out_param.measure,'uniform ball') % normal-to-ball transformation
            f = @(t) f(gail.domain_balls_spheres.ball_psi_2(t, out_param.dim, out_param.radius, hyperbox))*volume;% the psi function computes the transformation
        else % normal-to-sphere transformation
            f = @(t) f(gail.domain_balls_spheres.sphere_psi_2(t, out_param.dim, out_param.radius, hyperbox))*volume;% the psi function is the transformation
        end
        hyperbox = bsxfun(@plus, zeros(2, out_param.dim), [-inf; inf]);% the hyperbox must be the domain of the transformation, which is a unit box
        out_param.measure = 'normal';% then a normal distribution can be used
    end
end

%------------------------------------------------------------------------------

%now the integral will be calculated using uniform or normal measure
f=gail.transformIntegrand(f,hyperbox,out_param);
if strcmpi(out_param.measure,'uniform')% using uniformly distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(rand(nfun,out_param.dim)),out_param);
    % using meanMC_g to get the mean
elseif strcmpi(out_param.measure,'normal')% using normally distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(randn(nfun,out_param.dim)),out_param);
    % using meanMC_g to get the mean
end
out_param.time=toc(tstart); %elapsed time

% control the order of out_param
% out_param = orderfields(out_param, ...
% {'Yrand','abstol','reltol','tol','alpha','fudge', 'tau','hmu','time',...
% 'n1','nSig', 'n','nremain','nbudget','ntot','tbudget','var','kurtmax','exitflag'});

end

function [f,hyperbox,out_param] = cubMC_g_param(varargin)
% Parameter checking and parsing
default.measure = 'uniform';% default measure
default.transf = 1;% default transformation (box-to-ball or box-to-sphere, depending on the region of integration)
default.dim = 1;% default dimension
default.hyperbox = [zeros(1,default.dim);ones(1,default.dim)];% default hyperbox
default.abstol  = 1e-2;% default absolute error tolerance
default.reltol  = 1e-1;% default absolute error tolerance
default.alpha = 0.01;% default uncertainty
default.fudge = 1.2; % default variance inflation factor
default.nSig = 1e4; % default nSig initial sample size to estimate sigma
default.n1 = 1e4; % default n1 initial sample size to estimate Q
default.tbudget = 100;% default time budget in seconds
default.nbudget = 1e9; % default sample budget
default.flag = 0; % default value of parameter checking status
if isempty(varargin) % if no input, print error message and use the default setting
    help cubMC_g
    warning('GAIL:cubMC_g:fnotgiven',['f must be specified.'...
        'Now GAIL is using f = @(x) x.^2. '...
        'Integration hyperbox must be specified.'...
        'Now GAIL is using interval [0;1] with dimension 1.'])
    f = @(x) x.^2;
    hyperbox = default.hyperbox;
end
if numel(varargin)==1
    % if there is only function but no hyperbox input. Use default hyperbox.
    help cubMC_g
    warning('GAIL:cubMC_g:hyperboxnotgiven',...
        'the hyperbox must be specified, Now GAIL is using interval [0;1] with dimension 1')
    f = varargin{1};
    hyperbox = default.hyperbox;
end
if numel(varargin) >= 2
    f = varargin{1};
    hyperbox = varargin{2}; % the first input is function, the second input is hyperbox.
end

validvarargin=numel(varargin)>2;% check if there is any optional parameter input
if validvarargin
    in3=varargin{3}; % check the third input
    validvarargin=(isnumeric(in3) || isstruct(in3) || ischar(in3));
    % to see if it is numeric structure or character.
end
MATLABVERSION = gail.matlab_version;

if MATLABVERSION >= 8.3
  f_addParamVal= @addParameter;
else
  f_addParamVal = @addParamValue;
end

if ~validvarargin
% if there is no optional input, use default settings.
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.reltol = default.reltol;
    out_param.alpha = default.alpha;
    out_param.fudge = default.fudge;
    out_param.nSig = default.nSig;
    out_param.n1 = default.n1;
    out_param.tbudget = default.tbudget;
    out_param.nbudget = default.nbudget;
    out_param.flag = default.flag;
else % if there is some optional input
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'hyperbox',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        %if there are multiple inputs with only numeric, they should be put
        %in order.
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','uniform box',...
            'normal','Gaussian','uniform ball', 'uniform ball_box',...
            'uniform ball_normal','uniform sphere','uniform sphere_box',...
            'uniform sphere_normal'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.n1,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
        addOptional(p,'flag',default.flag,@isnumeric);
    else
        if isstruct(in3) %the input is structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end
        % if there are multiple inputs with name and numeric, they should
        % be put in order.
        f_addParamVal(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','uniform box',...
            'normal','Gaussian','uniform ball', 'uniform ball_box',...
            'uniform ball_normal','uniform sphere','uniform sphere_box',...
            'uniform sphere_normal'})));
        f_addParamVal(p,'abstol',default.abstol,@isnumeric);
        f_addParamVal(p,'reltol',default.reltol,@isnumeric);
        f_addParamVal(p,'alpha',default.alpha,@isnumeric);
        f_addParamVal(p,'fudge',default.fudge,@isnumeric);
        f_addParamVal(p,'nSig',default.nSig,@isnumeric);
        f_addParamVal(p,'n1',default.n1,@isnumeric);
        f_addParamVal(p,'tbudget',default.tbudget,@isnumeric);
        f_addParamVal(p,'nbudget',default.nbudget,@isnumeric);
        f_addParamVal(p,'flag',default.flag,@isnumeric);
    end
    parse(p,f,hyperbox,varargin{3:end})
    out_param = p.Results;
end
if (~gail.isfcn(f))
    warning('GAIL:cubMC_g:ynotfcn',...
        ['f must be a function handle.'...
        ' Now GAIL is using default f = @(x) x.^2 .'])
    %print warning message
    f = @(x) x.^2;
end

if any(isnan(hyperbox(:))); %check hyperbox for not a number
    out_param.exit=10; out_param = cubMC_g_err(out_param); return;
end

%hyperbox validation should depend on the measure.
if isfield(out_param,'measure'); % the sample measure
    out_param.measure=validatestring(out_param.measure,{'uniform','uniform box',...
        'normal','Gaussian','uniform ball', 'uniform ball_box',...
        'uniform ball_normal', 'uniform sphere', 'uniform sphere_box',...
            'uniform sphere_normal'});
    % simplifying out_param.measure and storing which transformation should
    % be applied
    if strcmpi(out_param.measure,'Gaussian')
        out_param.measure='normal';
    elseif strcmpi(out_param.measure,'uniform box')
        out_param.measure='uniform';
    elseif strcmpi(out_param.measure,'uniform ball') || strcmpi(out_param.measure,'uniform sphere')
        out_param.transf = default.transf;
    elseif strcmpi(out_param.measure,'uniform ball_box')
        out_param.measure='uniform ball';
        out_param.transf = 1;
    elseif strcmpi(out_param.measure,'uniform ball_normal')
        out_param.measure='uniform ball';
        out_param.transf = 2;
    elseif strcmpi(out_param.measure,'uniform sphere_box')
        out_param.measure='uniform sphere';
        out_param.transf = 1;
    elseif strcmpi(out_param.measure,'uniform sphere_normal')
        out_param.measure='uniform sphere';
        out_param.transf = 2;
    end
else
    out_param.measure=default.measure;
end

if strcmpi(out_param.measure,'uniform')||strcmpi(out_param.measure,'normal')
    [two, out_param.dim]=size(hyperbox); %hyperbox should be 2 x dimension
    if two==0 && isfield(out_param,'hyperbox');
        %if hyperbox specified through out_param structure
        hyperbox=out_param.hyperbox; %then get it from there
        [two, out_param.dim]=size(hyperbox); %and get the dimension
    end
    if two~=2 %if hyperbox is given as row vector for dimension 1, fix that
        if out_param.dim==2; out_param.dim=two; hyperbox=hyperbox';
        else out_param.exit=11; out_param = cubMC_g_err(out_param); return;
            %else, return an error
        end
    end

    hyperbox=[min(hyperbox,[],1); max(hyperbox,[],1)];
    %ensure left and right endpoints are in order
    if any(hyperbox(1,:)==hyperbox(2,:)); %hyperbox is a point in one direction
        out_param.exit=12; out_param = cubMC_g_err(out_param); return;
    end

    if strcmpi(out_param.measure,'uniform')&&~all(isfinite(hyperbox(:)))
        %cannot integrate on an infinite hyperbox with the uniform distribution
        out_param.exit=13; out_param = cubMC_g_err(out_param); return;
    end
    if strcmpi(out_param.measure,'normal')&&any(isfinite(hyperbox(:)))
        %must integrate on an infinite hyperbox with the normal distribution
        out_param.exit=14; out_param = cubMC_g_err(out_param); return;
    end
elseif strcmpi(out_param.measure,'uniform ball') || strcmpi(out_param.measure,'uniform sphere')
    [one, out_param.dim]=size(hyperbox); %hyperbox should be 1 x dimension
    if one==0 && isfield(out_param,'hyperbox');
        %if hyperbox specified through out_param structure
        hyperbox=out_param.hyperbox; %then get it from there
        [one, out_param.dim]=size(hyperbox); %and get the dimension
    end
    if one~=1 %if hyperbox is given as column vector for dimension 1, fix that
        if out_param.dim==1; out_param.dim=one; hyperbox=hyperbox';
        else out_param.exit=18; out_param = cubMC_g_err(out_param); return;
            %else, return an error
        end
    end

    if ~all(isfinite(hyperbox(:)))
        %the coordinates of the center of the ball and the radius must be
        %finite
        out_param.exit=15; out_param = cubMC_g_err(out_param); return;
    end
    out_param.radius = hyperbox(out_param.dim);% the last value of the hyperbox is the radius
    out_param.dim = out_param.dim - 1;% the last values doesn't count as coordinate
    if out_param.radius <= 0
       %the radius must be a positive real number
       out_param.exit=16; out_param = cubMC_g_err(out_param); return;
    end

    if strcmpi(out_param.measure,'uniform ball')
        if out_param.dim <= 0
            %the dimension must be a positive integer number
            out_param.exit=17; out_param = cubMC_g_err(out_param); return;
        end
    else
       if out_param.dim <= 1
            %the dimension must be >= 2 when measure is uniform sphere
            out_param.exit=19; out_param = cubMC_g_err(out_param); return;
       end
    end

    hyperbox = hyperbox(1,1:out_param.dim);
end

out_param.hyperbox=hyperbox; %copy hyperbox into the out_param structure

if out_param.flag == 0
    if (out_param.abstol < 0)
        %absolute error tolerance should be positive
        warning('GAIL:cubMC_g:abstolneg',...
            ['The absolute error tolerance should be positive; '...
            'We will take the absolute value of the absolute error tolerance provided.'])
        out_param.abstol = abs(out_param.abstol);
    end
if (out_param.reltol < 0 || out_param.reltol > 1)
    % relative error tolerance should be in [0,1]
    warning('GAIL:cubMC_g:reltolneg',...
        ['Relative error tolerance should be in [0,1]; ' ...
        'We will use the default value of the error tolerance 1e-1.'])
    out_param.reltol = default.reltol;
end
    if (out_param.alpha <= 0 ||out_param.alpha >= 1)
        %uncertainty should be 1 (0,1)
        warning('GAIL:cubMC_g:alphanot01',...
            ['The uncertainty should be in (0,1); '...
            'We will use the default value 1e-2.'])
        out_param.alpha = default.alpha;
    end
    if (out_param.fudge <= 1)
     %standard deviation inflation factor should be a number bigger than 1
        warning('GAIL:cubMC_g:fudgelessthan1',...
            ['The fudge factor should be bigger than 1; '...
            'We will use the default value 1.2.'])
        out_param.fudge = default.fudge;
    end
    if (~gail.isposge30(out_param.nSig))
        %the sample to estimate sigma should be a positive integer
        warning('GAIL:cubMC_g:nsignotposint',...
            ['The number nSig should a positive integer greater than 30; '
            'We will take the default value 1e4.'])
        out_param.nSig = default.nSig;
    end
    if (~gail.isposge30(out_param.n1))
    %initial sample size to estimate Q should be a positive integer
    warning('GAIL:cubMC_g:n1notposint',...
        ['The number n1 should a positive integer greater than 30; '...
        'We will use the default value 1e4.'])
    out_param.n1 = default.n1;
    end
    if (out_param.tbudget < 0)
        %the time budget in seconds should be positive
        warning('GAIL:cubMC_g:tbudgetneg',...
            ['Time budget should be positive; '...
            'We will use the absolute value of the time budget'])
        out_param.tbudget = abs(out_param.tbudget);
    end
    if (out_param.tbudget == 0)
        %the time budget in seconds should be positive
        warning('GAIL:cubMC_g:tbudget0',...
            ['Time budget should be positive rather than 0; '...
            'We will use the default value of the time budget 100 seconds.'])
        out_param.tbudget = default.tbudget;
    end
    if (~gail.isposge30(out_param.nbudget))
        %the sample budget should be a positive integer
        warning('GAIL:cubMC_g:nbudgetnotposint',...
            ['The number of sample budget should be a large positive integer;'...
            'We will take the default value of 1e9.'])
        out_param.nbudget = default.nbudget;
    end
    out_param.flag=2;
end
end

function  out_param =cubMC_g_err(out_param)
%Handles errors in cubMC_g and cubMC_g_param
%to give an exit with information
%out_param.exit = 0   success
%                 10  hyperbox does not contain numbers
%                 11  hyperbox not 2 x d when measure is uniform or normal
%                 12  hyperbox is only a point in one direction
%                 13  hyperbox is infinite when measure is uniform
%                 14  hyperbox is not doubly infinite when measure is normal
%                 15  hyperbox has an infinite coordinate for the center of
%                       the ball or a infinite radius for the ball
%                 16  the radius of the ball is a nonpositive real number
%                 17  the dimension of the ball is zero
%                 18  hyperbox not 1 x (d+1) when measure is uniform ball
%                 19  the dimension of the sphere is smaller than 2

if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 10; error('GAIL:cubMC_g:hyperboxnotnum',...
            'hyperbox must contain numbers.');
    case 11; error('GAIL:cubMC_g:hyperboxnot2d',...
            'hyperbox must be 2 x d when measure is ''uniform'' or ''normal''.');
    case 12; error('GAIL:cubMC_g:hyperboxnotlessthan2',...
            ['hyperbox must be more than a point in any coordinate'...
            ' direction when measure is ''uniform'' or ''normal''']);
    case 13; error('GAIL:cubMC_g:hyperboxnotfiniteforuniform',...
            'hyperbox must be finite when measure is uniform.');
    case 14; error('GAIL:cubMC_g:hyperboxnotinffornormal',...
            ['hyperbox must be infinite in both directions' ...
            ' when measure is normal']);
    case 15; error('GAIL:cubMC_g:hyperboxnotfiniteforuniformonaball',...
            'hyperbox must have only finite values for a distribution on a ball.');
    case 16; error('GAIL:cubMC_g:radiusnonpositive',...
            'the radius of the ball must be a positive real number.');
    case 17; error('GAIL:cubMC_g:dimensionequalszero',...
            'The number of dimentions of the ball must be a positive integer.');
    case 18; error('GAIL:cubMC_g:hyperboxnot1d',...
            'hyperbox must be 1 x (d+1) when measure is ''uniform ball''.');
    case 19; error('GAIL:cubMC_g:dimensionsmallerthan2',...
            'The number of dimentions of the sphere must be greater than or equal to 2.');
end
end
