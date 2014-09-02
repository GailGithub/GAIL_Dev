function [Q,out_param] = cubMC_g(varargin)
%CUBMC_G Monte Carlo method to evaluate a multidimensional integral to
%within a specified generalized error tolerance tolfun = max(abstol,
%reltol|I|) with guaranteed confidence level 1-alpha.
%
%   [Q,out_param] = CUBMC_G(f,hyperbox) estimates the integral of f over
%   hyperbox to within an specified generalized error tolerance tolfun =
%   max(abstol, reltol|I|) and with guaranteed confidence level 99%. Input
%   f is a function handle. The function f should accept an n x d matrix
%   input, where d is the dimension of the hyperbox, and n is the number of
%   points being evaluated simultaneously. The input hyperbox is a 2 x d
%   matrix, where the first row corresponds to the lower limits and the
%   second row corresponds to the upper limits.
% 
%   Q = CUBMC_G(f,hyperbox,measure,abstol,reltol,alpha,fudge,nSig,n1,
%   tbudget,nbudget,checked) estimates the integral of f over hyperbox with
%   respect to a given measure. The answer is given to within an specified
%   generalized error tolerance tolfun with guaranteed confidence level
%   1-alpha. All parameters should be input in the order specified above.
%   If an input is not specified, the default value is used.
% 
%   Q = CUBMC_G(f,hyperbox,'measure','uniform','abstol',abstol,'reltol',
%   reltol,'alpha',alpha,'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',
%   tbudget,'nbudget',nbudget,'checked',checked) estimates the integral of
%   f over hyperbox to within a specified generalized error tolerance
%   tolfun with guaranteed confidence level 1-alpha. All the field-value
%   pairs are optional and can be supplied in different order. If an input
%   is not specified, the default value is used.
% 
%   [Q out_param] = CUBMC_G(f,hyperbox,in_param) estimates the integral of
%   f over hyperbox to within a specified generalized error tolerance
%   tolfun with the given parameters in_param and produce output parameters
%   out_param.
% 
%   Input Arguments
%
%     f --- the integrand.
% 
%     hyperbox --- the integration hyperbox. The default value is
%     [zeros(1,d); ones(1,d)], the default d is 1.
% 
%     in_param.measure --- the measure for generating the random variable, the
%     default is uniform. The other measure could be handled is normal/Gaussian.
% 
%     in_param.abstol --- the absolute error tolerance, the default value is
%     1e-2.
%
%     in_param.reltol --- the relative error tolerance, the default value is
%     1e-1.
% 
%     in_param.alpha --- the uncertainty, the default value is 1%.
% 
%     in_param.fudge --- the standard deviation inflation factor, the default
%     value is 1.1.
%
%     in_param.nSig --- initial sample size for estimating the sample
%     variance, the default value is 1e3.
% 
%     in_param.n1 --- initial sample size for estimating the sample
%     mean, the default value is 1e4.
% 
%     in_param.tbudget --- the time budget to do the estimation, the
%     default value is 100 seconds.
% 
%     in_param.nbudget --- the sample budget to do the estimation, the
%     default value is 1e8.
% 
%     in_param.checked --- the value corresponds to parameter checking status.
%                         0   not checked
%                         1   checked by meanMC_g
%                         2   checked by cubMC_g
%
%   Output Arguments
%
%     Q --- the estimated value of the integral.
% 
%     out_param.n --- sample used in each iteration.
%
%     out_param.ntot --- total sample used.
%
%     out_param.tau --- the iteration step.
%
%     out_param.mu --- estimated mean in each iteration
%
%     out_param.tol --- the tolerance for each iteration
%  
%     out_param.kurtmax --- the upper bound on modified kurtosis.
% 
%     out_param.time --- the time elapsed.
%
%     out_param.exit --- the state of program when exiting.
%                       0   success
%                       1   Not enough samples to estimate the mean.
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
% Guarantee
% This algorithm attampts to calculate the mean of a random variable to a
% certain tolerance with guaranteed confidence level 1-alpha. If the
% algorithm terminated without showing any warning messages and provide an
% answer \hat{mu}, then the follow inequality would be satisfied:
% 
% Pr(|mu-\hat{mu}| <= max(abstol,reltol|mu|)) >= 1-alpha
%
% where abstol is the absolute error tolerance and reltol is the relative
% error tolerance, if the true mean mu is rather small as well as the
% reltol, then the abstol would be satisfied, and vice versa. 
%
% The cost of the algorithms is also bounded above by N_up, which is in
% terms of abstol, reltol, nSig, n_1, fudge, alpha_sigma, kmax, beta.
% And the following inequality would hold:
% 
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
% 
%  Examples
% 
% Example 1:
% If no parsing any parameter, help text will show up as following
% >> cubMC_g
% ***Monte Carlo method to estimate***
%
%
% Example 2:
% Estimate the integral with integrand f(x) = sin(x) in the interval [1;2]
% 
% >> f=@(x) sin(x);interval = [1;2];
% >> Q = cubMC_g(f,interval,'uniform',1e-3,1e-2)
% Q = 0.95***
% 
% 
% Example 3: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [0 0;1 1], where x is a vector x = [x1 x2].
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2);hyperbox = [0 0;1 1];
% >> Q = cubMC_g(f,hyperbox,'measure','uniform','abstol',1e-3,'reltol',1e-13)
% Q = 0.55***
% 
% 
% Example 4: 
% Estimate the integral with integrand f(x) = 2^d*prod(x1*x2*...*xd)+0.555 in the
% hyperbox [zeros(1,d);ones(1,d)], where x is a vector x = [x1 x2 ... xd].
% 
% >> d=3;f=@(x) 2^d*prod(x,2)+0.555;hyperbox =[zeros(1,d);ones(1,d)];
% >> in_param.abstol = 1e-3;in_param.reltol=1e-3;
% >> Q = cubMC_g(f,hyperbox,in_param)
% Q = 1.5***
% 
%
% Example 5: 
% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [-inf -inf;inf inf], where x is a vector x = [x1 x2].
% 
% >> f=@(x) exp(-x(:,1).^2-x(:,2).^2);hyperbox = [-inf -inf;inf inf];
% >> Q = cubMC_g(f,hyperbox,'normal',1e-3,1e-2)
% Q = 0.33***
% 
% 
%   See also FUNAPPX_G, INTEGRAL_G, MEANMC_G
% 
%  References
%
%   [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, Guaranteed
%   conservative fixed width confidence intervals via Monte Carlo sampling,
%   Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F. Y. Kuo, G. W.
%   Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin, 2014, to appear,
%   arXiv:1208.4318 [math.ST]
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
%   Yizhi Zhang, "GAIL: Guaranteed Automatic Integration Library (Version
%   1.3.0)" [MATLAB Software], 2014. Available from
%   http://code.google.com/p/gail/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above paper and software.


tstart=tic;
[f,hyperbox,out_param] = cubMC_g_param(varargin{:});%check validity of inputs
f=gail.transformIntegrand(f,hyperbox,out_param); 
if strcmp(out_param.measure,'uniform')% the using uniformly distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(rand(nfun,out_param.dim)),out_param);
   % out_param.Q=Q;% using meanMC_g to get the mean 
else strcmp(out_param.measure,'normal')% using normally distributed samples
    [Q,out_param] = meanMC_g(@(nfun)f(randn(nfun,out_param.dim)),out_param);
    %out_param.Q=Q;% using meanMC_g to get the mean
end
out_param.time=toc(tstart); %elapsed time
end

function [f,hyperbox,out_param] = cubMC_g_param(varargin)
% Parameter checking and parsing
default.measure = 'uniform';% default measure
default.dim = 1;% default dimension
default.hyperbox = [zeros(1,default.dim);ones(1,default.dim)];% default hyperbox
default.abstol  = 1e-2;% default absolute error tolerance
default.reltol  = 1e-1;% default absolute error tolerance
default.alpha = 0.01;% default uncertainty
default.fudge = 1.1; % default variance inflation factor
default.nSig = 1e3; % default nSig initial sample size to estimate sigma
default.n1 = 1e4; % default n1 initial sample size to estimate Q
default.tbudget = 100;% default time budget
default.nbudget = 1e8; % default sample budget
default.checked = 0; % default value of parameter checking status
if isempty(varargin) % if no input, print error message and use the default setting
    help cubMC_g
    warning('MATLAB:cubMC_g:fnotgiven',['f must be specified.'...
        'Now GAIL is using f = @(x) x.^2. '...
        'Integration hyperbox must be specified.'...
        'Now GAIL is using interval [0;1] with dimension 1.'])
    f = @(x) x.^2;
    hyperbox = default.hyperbox;
elseif numel(varargin)==1
    % if there is only function but no hyperbox input. Use default hyperbox.
    help cubMC_g
    warning('MATLAB:cubMC_g:hyperboxnotgiven',...
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
    out_param.reltol = default.reltol;
    out_param.alpha = default.alpha;    
    out_param.fudge = default.fudge;
    out_param.nSig = default.nSig;
    out_param.n1 = default.n1;
    out_param.tbudget = default.tbudget;
    out_param.nbudget = default.nbudget;
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
        addOptional(p,'reltol',default.reltol,@isnumeric);
        addOptional(p,'alpha',default.alpha,@isnumeric);
        addOptional(p,'fudge',default.fudge,@isnumeric);
        addOptional(p,'nSig',default.nSig,@isnumeric);
        addOptional(p,'n1',default.nSig,@isnumeric);
        addOptional(p,'tbudget',default.tbudget,@isnumeric);
        addOptional(p,'nbudget',default.nbudget,@isnumeric);
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
        addParamValue(p,'reltol',default.reltol,@isnumeric);
        addParamValue(p,'alpha',default.alpha,@isnumeric);
        addParamValue(p,'fudge',default.fudge,@isnumeric);
        addParamValue(p,'nSig',default.nSig,@isnumeric);
        addParamValue(p,'n1',default.n1,@isnumeric);
        addParamValue(p,'tbudget',default.tbudget,@isnumeric);
        addParamValue(p,'nbudget',default.nbudget,@isnumeric);
        addParamValue(p,'checked',default.checked,@isnumeric); 
    end
    parse(p,f,hyperbox,varargin{3:end})
    out_param = p.Results;
end
if any(isnan(hyperbox(:))); %check hyperbox for not a number
    out_param.exit=10; out_param = cubMC_g_err(out_param); return; 
end
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
out_param.hyperbox=hyperbox; %copy hyperbox into the out_param structure

if isfield(out_param,'measure'); % the sample measure
    out_param.measure=validatestring(out_param.measure,{'uniform','normal','Gaussian'});
    if strcmpi(out_param.measure,'Gaussian')
        out_param.measure='normal'; 
    end
else
    out_param.measure=default.measure;
end
if strcmp(out_param.measure,'uniform')&&~all(isfinite(hyperbox(:)))
    %cannot integrate on an infinite hyperbox with the uniform distribution
    out_param.exit=13; out_param = cubMC_g_err(out_param); return;
end
if strcmp(out_param.measure,'normal')&&any(isfinite(hyperbox(:)))
    %must integrate on an infinite hyperbox with the normal distribution
    out_param.exit=14; out_param = cubMC_g_err(out_param); return;
end
if out_param.checked == 0
    if (out_param.abstol <= 0) 
        %absolute error tolerance should be positive
        warning('MATLAB:cubMC_g:abstolneg',...
            'the absolute error tolerance should be larger than 0, use the absolute value.')
        out_param.abstol = abs(out_param.abstol);
    end
if (out_param.reltol <= 0 || out_param.reltol >= 1)
    % relative error tolerance should be in (0,1)
    warning('MATLAB:cubMC_g:reltolneg',...
        ['Relative error tolerance should be between 0 and 1, ' ...
        'use the default value of the error tolerance'])
    out_param.abstol = abs(out_param.abstol);
end    
    if (out_param.alpha <= 0 ||out_param.alpha >= 1) 
        %uncertainty should be 1 (0,1)
        warning('MATLAB:cubMC_g:alphanot01',...
            ['the uncertainty should be less than 1 and bigger than 0, '...
            'use the the default value.'])
        out_param.alpha = default.alpha;
    end
    if (out_param.fudge <= 1) 
        %standard deviation inflation factor should be a number bigger than
        %1
        warning('MATLAB:cubMC_g:fudgelessthan1',...
            'the fudge factor should be bigger than 1, use the default value.')
        out_param.fudge = default.fudge;
    end    
    if (~gail.isposint(out_param.nSig))
        %the sample to estimate sigma should be a positive integer
        warning('MATLAB:cubMC_g:nsignotposint',...
            ['the number nSig should a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.nSig = ceil(abs(out_param.nSig));
    end
    if (~gail.isposint(out_param.n1)) 
    %initial sample size to estimate Q should be a posotive integer
    warning('MATLAB:cubMC_g:n1notposint',...
        ['the number n1 should a positive integer, '...
        'take the absolute value and ceil.'])
    out_param.n1 = ceil(abs(out_param.n1));
    end
    if (out_param.tbudget <= 0) % time budget should be positive
        warning('MATLAB:cubMC_g:tbudgetlneg',...
            ['Time budget should be bigger than 0, '...
            'use the absolute value of time budget'])
        out_param.tbudget = abs(out_param.tbudget);
    end
    if (~gail.isposint(out_param.nbudget)) 
        % sample budget should be a postitive integer
        warning('MATLAB:cubMC_g:nbudgetnotposint',...
            ['the number of sample budget should be a positive integer,'...
            'take the absolute value and ceil.'])
        out_param.nbudget = ceil(abs(out_param.nbudget));
    end
    out_param.checked=2;
end
end

function  out_param =cubMC_g_err(out_param)
%Handles errors in cubMC_g and cubMC_g_param
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
    case 10; error('MATLAB:cubMC_g:hyperboxnotnum',...
            'hyperbox must contain numbers.');
    case 11; error('MATLAB:cubMC_g:hyperboxnot2d',...
            'hyperbox must be 2 x d.');
    case 12; error('MATLAB:cubMC_g:hyperboxnotlessthan2',...
            'hyperbox must be more than a point in any coordinate direction.');
    case 13; error('MATLAB:cubMC_g:hyperboxnotfiniteforuniform',...
            'hyperbox must be finite when measure is uniform.');
    case 14; error('MATLAB:cubMC_g:hyperboxnotinffornormal',...
            ['hyperbox must be infinite in both directions' ...
        ' when measure is normal']);
end
end
