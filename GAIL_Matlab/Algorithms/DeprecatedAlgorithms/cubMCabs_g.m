function [Q,out_param] = cubMCabs_g(varargin)
%CUBMCABS_G Monte Carlo method to evaluate a multidimensional integral to
%within a specified absolute error tolerance with guaranteed confidence
%level 1-alpha.
%
%   [Q,out_param] = CUBMCABS_G(f,hyperbox) estimates the integral of f over
%   hyperbox to within an specified absolute error tolerance 1e-2 and with
%   guaranteed confidence level 99%. Input f is a function handle. The
%   function f should accept an n x d matrix input, where d is the dimension
%   of the hyperbox, and n is the number of points being evaluated
%   simultaneously. The input hyperbox is a 2 x d matrix, where the first row
%   corresponds to the lower limits and the second row corresponds to the
%   upper limits.
% 
%   Q = CUBMCABS_G(f,hyperbox,measure,abstol,alpha,n_sigma,fudge,tbudget,nbudget,npcmax,checked)
%   estimates the integral of f over hyperbox with respect to a given
%   measure. The answer is given to within an specified absolute error
%   tolerance abstol with guaranteed confidence level 1-alpha. All parameters
%   should be input in the order specified above. If an input is not
%   specified, the default value is used.
% 
%   Q = CUBMCABS_G(f,hyperbox,'measure','uniform','abstol',abstol,'alpha',alpha,
%   'n_sigma',n_sigma,'fudge',fudge,'tbudget',tbudget,'nbudget',nbudget,
%   'npcmax',npcmax,'checked',checked) estimates the integral of f over
%   hyperbox to within an specified absolute error tolerance abstol with
%   guaranteed confidence level 1-alpha. All the field-value pairs are
%   optional and can be supplied in different order. If an input is not
%   specified, the default value is used.
% 
%   Q = CUBMCABS_G(f,hyperbox,in_param) estimates the integral of f over
%   hyperbox to within an specified absolute error tolerance in_param.abstol
%   with guaranteed confidence level 1-in_param.alpha. If a field is not
%   specified, the default value is used.
% 
%   Input Arguments
%
%     f --- the integrand.
% 
%     hyperbox --- the integration hyperbox. The default value is
%     [zeros(1,d); ones(1,d)], the default d is 1.
% 
%     in_param.measure --- the measure for generating the random variable, the
%     default is uniform. The other measure we could handle is normal/Gaussian.
% 
%     in_param.abstol --- the absolute error tolerance, the default value is
%     1e-2.
% 
%     in_param.alpha --- the uncertainty, the default value is 1%.
% 
%     in_param.n_sigma --- initial sample size for estimating the sample
%     variance, the default value is 1e4.
% 
%     in_param.fudge --- the standard deviation inflation factor, the default
%     value is 1.2.
% 
%     in_param.tbudget --- the time budget to do the two-stage estimation,
%     the default value is 100 seconds.
% 
%     in_param.nbudget --- the sample budget to do the two-stage estimation,
%     the default value is 1e8.
% 
%     in_param.npcmax --- number of elements in an array of optimal size to
%     calculate the mean, the default value is 1e6.
% 
%     in_param.checked --- the value corresponds to parameter checking status.
%                         0   not checked
%                         1   checked by meanMCabs_g
%                         2   checked by cubMCabs_g
%
%   Output Arguments
%
%     Q --- the estimated value of the integral.
% 
%     out_param.time_n_sigma_predict --- the estimated time to get n_sigma
%     samples.
% 
%     out_param.n_left_predict --- using the time left to predict the number
%     of samples left.
% 
%     out_param.nmax --- the maximum sample budget to estimate the mean, it
%     comes from both the sample budget and the time budget.
% 
%     out_param.var --- the sample variance.
% 
%     out_param.kurtmax --- the upper bound on modified kurtosis.
% 
%     out_param.time --- the time elapsed.
% 
%     out_param.n_mu --- the sample size that needed to estimate the mean,
%     which comes from Berry-Esseen inequality and Chebyshev inequality.
% 
%     out_param.n --- the total sample size needed to do the two stage
%     estimation.
%
%     out_param.exit --- the state of program when exiting.
%                       0   success
%                       1   No enough samples to estimate the mean.
%                       2   Initial try out time costs more than
%                           10% of time budget. 
%                       3   The estimated time for estimating variance 
%                           is bigger than half of the time budget.
%                       10  hyperbox does not contain numbers.
%                       11  hyperbox not 2 x d.
%                       12  hyperbox is only a point in one direction.
%                       13  hyperbox is infinite when measure is uniform.
%                       14  hyperbox is not doubly infinite when measure
%                           is normal.
% 
%  Guarantee
% 
% If the modified kurtosis of the integrand, f, is less than the kurtmax,
% which is defined in terms of the uncertainty, alpha, the sample size to
% estimate variance, n_sigma, and the standard deviation inflation factor,
% fudge, then the inequality
% 
% Pr(|I-Q| <= abstol) >= 1-alpha 
% 
% holds. Here I is the true integral (or mean) of f, and Q is the output
% of CUBMCABS_G. 
% 
% The cost of the two-stage algorithm also satisfies the inequality
% 
% Pr (N_tot <= N_up) >= 1-beta
% 
% where N_tot is the total cost of samples, N_up is the upper bound on the
% cost, which is roughly proportional to sigma^2/abstol^2, beta is the
% level of uncertainty on the cost. For details, please refer to [1].
% 
%  Examples
% 
% Example 1:
% Estimate the integral with integrand f(x) = sin(x) in the interval [1;2]
% 
% >> f=@(x) sin(x);interval = [1;2];
% >> Q = cubMCabs_g(f,interval,'uniform',1e-3)
% Q = 0.95***
% 
% 
% Example 2: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [0 0;1 1], where x is a vector x = [x1 x2].
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2);hyperbox = [0 0;1 1];
% >> Q = cubMCabs_g(f,hyperbox,'uniform',1e-3)
% Q = 0.55***
% 
% 
% Example 3: 
% Estimate the integral with integrand f(x) = 2^d*prod(x1*x2*...*xd)+0.555 in the
% hyperbox [zeros(1,d);ones(1,d)], where x is a vector x = [x1 x2 ... xd].
% 
% >> d=3;f=@(x) 2^d*prod(x,2)+0.555;hyperbox = [zeros(1,d);ones(1,d)];
% >> Q = cubMCabs_g(f,hyperbox,'uniform',1e-3)
% Q = 1.5***
% 
%
% Example 4: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [-inf -inf;inf inf], where x is a vector x = [x1 x2].
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2);hyperbox = [-inf -inf;inf inf];
% >> Q = cubMCabs_g(f,hyperbox,'normal',1e-3)
% Q = 0.33***
% 
% 
%   See also FUNAPPX_G, INTEGRAL_G, MEANMCABS_G
% 
%  References
%
%   [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
%   Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
%   arXiv:1208.4318 [math.ST]
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
%   Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
%   1.3.0)" [MATLAB Software], 2014. Available from
%   http://gailgithub.github.io/GAIL_Dev/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.


tstart=tic;
[f,hyperbox,out_param] = cubMCabs_g_param(varargin{:});%check validity of inputs
f=transformIntegrand(f,hyperbox,out_param); 
if strcmpi(out_param.measure,'uniform')% the using uniformly distributed samples
    [Q,out_param] = meanMCabs_g(@(nfun)f(rand(nfun,out_param.dim)),out_param);
else strcmpi(out_param.measure,'normal')% using normally distributed samples
    [Q,out_param] = meanMCabs_g(@(nfun)f(randn(nfun,out_param.dim)),out_param);
end
out_param.time=toc(tstart); %elapsed time
end

function newf=transformIntegrand(oldf,hyperbox,out_param)
% Transform integrand linearly so that the hyperbox would not be changed
    if strcmpi(out_param.measure,'uniform') %uniform measure
        a=hyperbox(1,:); %left endpoint
        b=hyperbox(2,:); %right endpoint
    if all(a==0) && all(b==1) %no change needed
        newf=oldf; 
    else %transform points and integrand
        bmina=b-a; %hyperbox width
        volbox=prod(bmina); %volume of the hyperbox
        newf=@(x) oldf(x.*repmat(bmina,size(x,1),1)+repmat(a,size(x,1),1))...
            .*volbox;
       %stretch and shift, then multiply by volume
    end
    elseif strcmpi(out_param.measure,'normal')
        newf=oldf;% no change if it is normal measure.
    end   
end

function [f,hyperbox,out_param] = cubMCabs_g_param(varargin)
% Parameter checking and parsing
default.measure = 'uniform';% default measure
default.dim = 1;% default dimension
default.hyperbox = [zeros(1,default.dim);ones(1,default.dim)];% default hyperbox
default.abstol  = 1e-2;% default absolute error tolerance
default.alpha = 0.01;% default uncertainty
default.n_sigma = 1e4; % default n_sigma
default.fudge = 1.2; % default variance inflation factor
default.tbudget = 100;% default time budget
default.nbudget = 1e8; % default sample budget
default.npcmax = 1e6;% default n piece maximum
default.checked = 0; % default value of parameter checking status
if isempty(varargin) % if no input print error message and use the default setting
    help cubMCabs_g
    warning('GAIL:cubMCabs_g:fnotgiven',['f must be specified.'...
        'Now GAIL is using f = @(x) x.^2. '...
        'Integration hyperbox must be specified.'...
        'Now GAIL is using interval [0;1] with dimension 1.'])
    f = @(x) x.^2;
    hyperbox = default.hyperbox;
elseif numel(varargin)==1
    % if there is only function but no hyperbox input. Use default hyperbox.
    help cubMCabs_g
    warning('GAIL:cubMCabs_g:hyperboxnotgiven',...
        'the hyperbox must be specified, Now GAIL is using interval [0;1] with dimension 1')
    f = varargin{1};
    hyperbox = default.hyperbox;    
else
    f = varargin{1};
    hyperbox = varargin{2}; % the first input is function, the second input is hyperbox.
end
    
validvarargin=numel(varargin)>2;% check if there is any optional parameter input
if validvarargin
    in3=varargin{3}; % check the third input
    validvarargin=(isnumeric(in3) || isstruct(in3) || ischar(in3));
    % to see if it is numeric structure or character.
end

if ~validvarargin
% if there is no optional input, use default settings.
    out_param.measure = default.measure;
    out_param.abstol = default.abstol;
    out_param.alpha = default.alpha;
    out_param.n_sigma = default.n_sigma;
    out_param.fudge = default.fudge;
    out_param.tbudget = default.tbudget;
    out_param.nbudget = default.nbudget;
    out_param.npcmax = default.npcmax;
    out_param.checked = default.checked;
else % if there is some optional input 
    p = inputParser;
    addRequired(p,'f',@gail.isfcn);
    addRequired(p,'hyperbox',@isnumeric);
    if isnumeric(in3) || ischar(in3)
        %if there are multiple inputs with only numeric, they should be put
        %in order.
        addOptional(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
        addOptional(p,'abstol',default.abstol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'n_sigma',default.n_sigma,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
        addOptional(p,'npcmax',default.npcmax,@isnumeric);
        addOptional(p,'checked',default.checked,@isnumeric);
    else
        if isstruct(in3) %the input is structure
            p.StructExpand = true;
            p.KeepUnmatched = true;
        end % if there are multiple inputs with name and numeric, they
    % could be put not in order
        addParamValue(p,'measure',default.measure,...
            @(x) any(validatestring(x, {'uniform','normal','Gaussian'})));
        addParamValue(p,'abstol',default.abstol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'n_sigma',default.n_sigma,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
        addParamValue(p,'npcmax',default.npcmax,@isnumeric); 
        addParamValue(p,'checked',default.checked,@isnumeric); 
    end
    parse(p,f,hyperbox,varargin{3:end})
    out_param = p.Results;
end
if any(isnan(hyperbox(:))); %check hyperbox for not a number
    out_param.exit=10; out_param = cubMCabs_g_err(out_param); return; 
end
[two, out_param.dim]=size(hyperbox); %hyperbox should be 2 x dimension
if two==0 && isfield(out_param,'hyperbox'); 
    %if hyperbox specified through out_param structure
    hyperbox=out_param.hyperbox; %then get it from there
    [two, out_param.dim]=size(hyperbox); %and get the dimension
end
if two~=2 %if hyperbox is given as row vector for dimension 1, fix that
    if out_param.dim==2; out_param.dim=two; hyperbox=hyperbox';
    else out_param.exit=11; out_param = cubMCabs_g_err(out_param); return; 
        %else, return an error
    end
end
hyperbox=[min(hyperbox,[],1); max(hyperbox,[],1)]; 
%ensure left and right endpoints are in order
if any(hyperbox(1,:)==hyperbox(2,:)); %hyperbox is a point in one direction
    out_param.exit=12; out_param = cubMCabs_g_err(out_param); return;
end
out_param.hyperbox=hyperbox; %copy hyperbox into the out_param structure

if isfield(out_param,'measure'); % the sample measure
    out_param.measure=validatestring(out_param.measure,{'uniform','normal','Gaussian'});
    if strcmpi(out_param.measure,'Gaussian')
        out_param.measure='normal'; 
    end
else
    out_param.measure=default.measure;
end
if strcmpi(out_param.measure,'uniform')&&~all(isfinite(hyperbox(:)))
    %cannot integrate on an infinite hyperbox with the uniform distribution
    out_param.exit=13; out_param = cubMCabs_g_err(out_param); return;
end
if strcmpi(out_param.measure,'normal')&&any(isfinite(hyperbox(:)))
    %must integrate on an infinite hyperbox with the normal distribution
    out_param.exit=14; out_param = cubMCabs_g_err(out_param); return;
end
if out_param.checked == 0
    if (out_param.abstol <= 0) %absolute error tolerance
        warning('GAIL:cubMCabs_g:abstolneg',...
            'the absolute error tolerance should be larger than 0, use the absolute value.')
        out_param.abstol = abs(out_param.abstol);
    end
    if (out_param.alpha <= 0 ||out_param.alpha >= 1) %uncertainty
        warning('GAIL:cubMCabs_g:alphanot01',...
            ['the uncertainty should be less than 1 and bigger than 0, '...
            'use the the default value.'])
        out_param.alpha = default.alpha;
    end
    if (~gail.isposint(out_param.n_sigma)) %the sample to estimate sigma
        warning('GAIL:cubMCabs_g:nsignotposint',...
            ['the number n_sigma should a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.n_sigma = ceil(abs(out_param.n_sigma));
    end
    if (out_param.fudge <= 1) %standard deviation inflation factor/fudge factor
        warning('GAIL:cubMCabs_g:fudgelessthan1',...
            'the fudge factor should be bigger than 1, use the default value.')
        out_param.fudge = default.fudge;
    end
    if (out_param.tbudget <= 0) % time budget
        warning('GAIL:cubMCabs_g:tbudgetlneg',...
            ['Time budget should be bigger than 0, '...
            'use the absolute value of time budget'])
        out_param.tbudget = abs(out_param.tbudget);
    end
    if (~gail.isposint(out_param.nbudget)) % sample budget should be a postitive integer
        warning('GAIL:cubMCabs_g:nbudgetnotposint',...
            ['the number of sample budget should be a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.nbudget = ceil(abs(out_param.nbudget));
    end
    if (~gail.isposint(out_param.npcmax))
        % maximum number of scalar values of x per vector should be a positive integer
        warning('GAIL:cubMCabs_g:npcmaxnotposint',...
            ['the number of each piece of the samples should be' ...
            'a positive integer, take the absolute value and ceil.'])
        out_param.npcmax = ceil(abs(out_param.npcmax));
    end
    out_param.checked=2;
end
end

function  out_param =cubMCabs_g_err(out_param)
%Handles errors in cubMCabs_g and cubMCabs_g_param
%to give an exit with information
%out_param.exit = 0   success
%                 10  hyperbox does not contain numbers
%                 11  hyperbox not 2 x d
%                 12  hyperbox is only a point in one direction
%                 13  hyperbox is infinite when measure is uniform
%                 14  hyperbox is not doubly infinite when measure is normal
if ~isfield(out_param,'exit'); return; end
if out_param.exit==0; return; end
switch out_param.exit
    case 10; error('GAIL:cubMCabs_g:hyperboxnotnum',...
            'hyperbox must contain numbers.');
    case 11; error('GAIL:cubMCabs_g:hyperboxnot2d',...
            'hyperbox must be 2 x d.');
    case 12; error('GAIL:cubMCabs_g:hyperboxnotlessthan2',...
            'hyperbox must be more than a point in any coordinate direction.');
    case 13; error('GAIL:cubMCabs_g:hyperboxnotfiniteforuniform',...
            'hyperbox must be finite when measure is uniform.');
    case 14; error('GAIL:cubMCabs_g:hyperboxnotinffornormal',...
            ['hyperbox must be infinite in both directions' ...
        ' when measure is normal']);
end
end
