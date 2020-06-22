%CUBBAYESNET_G Bayesian cubature method to estimate the integral of a
% random variable using digital nets over a d-dimensional region within a
% specified generalized error tolerance with guarantees under Bayesian
% assumptions. Currently, only Sobol points are supported.
%
%   [OBJ,Q] = CUBBAYESNET_G(f,dim,'absTol',absTol,'relTol',relTol,
%   'order',order,'arbMean',arbMean) initializes the object with the given
%   parameters and also returns an estimate of integral Q.
%
%   [Q,OutP] = COMPINTEG(OBJ) estimates the integral of f over hyperbox
%   [0,1]^d using digital nets (Sobol points) to within a specified
%   generalized error tolerance, tolfun = max(abstol, reltol*| I |), i.e.,
%   | I - Q | <= tolfun with confidence of at least 99%, where I is the
%   true integral value, Q is the estimated integral value, abstol is the
%   absolute error tolerance, and reltol is the relative error tolerance.
%   Usually the reltol determines the accuracy of the estimation; however,
%   if | I | is rather small, then abstol determines the accuracy of the
%   estimation.
%
%   It is recommended to use COMPINTEG for estimating the integral
%   repeatedly after the object initialization.
%
%   OutP is the structure holding additional output params, more details
%   provided below. Input f is a function handle that accepts an n x d
%   matrix input, where d is the dimension of the hyperbox, and n is the
%   number of points being evaluated simultaneously.
%
%   The following additional input parameter passing styles also supported:
%
%   [OBJ,Q] = CUBBAYESNET_G(f,dim); estimates the integral of f over
%   hyperbox [0,1]^dim using digital nets (Sobol points). All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.
%
%   [OBJ,Q] = CUBBAYESNET_G(f,dim,absTol,relTol); estimates the integral
%   of f over hyperbox [0,1]^dim using digital nets (Sobol points). All parameters
%   should be input in the order specified above. The answer is given within
%   the generalized error tolerance tolfun. All other input parameters
%   are initialized with default values as given below.
%
%   [OBJ,Q] = CUBBAYESNET_G(f,dim,inParms); estimates the integral
%   of f over hyperbox [0,1]^dim digital nets (Sobol points).
%   The structure inParams shall hold the optional input parameters.
%
%   Input Arguments
%
%     f --- the integrand.
%     dim --- number of dimensions of the integrand.
%
%   Optional Input Arguments
%
%     absTol --- absolute error tolerance | I - Q | <= absTol. Default is 0.01
%     relTol --- relative error tolerance | I - Q | <= I*relTol. Default is 0
%     arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean. Default is 'true'
%     alpha --- confidence level for a credible interval of Q. Default is 0.01
%     mmin --- min number of samples to start with: 2^mmin. Default is 8
%     mmax --- max number of samples allowed: 2^mmax. Default is 20
%
%   Output Arguments
%
%    n --- number of samples used to compute the integral of f.
%    time --- time to compute the integral in seconds.
%    exitFlag --- indicates the exit condition of the algorithm:
%                      1 - integral computed within the error tolerance and
%                      without exceeding max sample limit 2^mmax
%                      2 - used max number of samples and yet not met the
%                      error tolerance
%    ErrBd  --- estimated integral error | I - Q |
%    optParams --- optional parameters useful to debug and get better
%                  understanding of the algorithm
%    optParams.aMLEAll ---- returns the shape parameters computed
%
%
%  Guarantee
%
%   This algorithm attempts to calculate the integral of function f over the
%   hyperbox [0,1]^dim to a prescribed error tolerance tolfun:=
%   max(abstol,reltol*| I |) with guaranteed confidence level, e.g., 99% when
%   alpha=0.5%. If the algorithm terminates without showing any warning
%   messages and provides an answer Q, then the following inequality would be
%   satisfied:
%
%   Pr(| Q - I | <= tolfun) = 99%
%
%   Please refer to our paper [1] for detailed arguments and proofs.
%
%
%  Examples
%
%   Example 1:
%   If no parameters are parsed, help text will show up as follows:
%   >> help cubBayesNet_g
%   ***Bayesian cubature method to estimate the integral ***
%
%
%   Example 2: Quadratic
%   Estimate the integral with integrand f(x) = x.^2 over the interval [0,1]
%   with default parameters: order=1, abstol=0.01, relTol=0
%
%   >> warning('off','GAIL:cubBayesNet_g:fdnotgiven')
%   >> [~,muhat] = cubBayesNet_g;
%   >> exactInteg = 1.0/3;
%   >> warning('on','GAIL:cubBayesNet_g:fdnotgiven')
%   >> check = double(abs(exactInteg-muhat) < 0.01)
%   check = 1
%
%
%   Example 3: ExpCos
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over
%   the interval [0,1] with parameters: order=1, abstol=0.001, relTol=0.01
%
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; absTol=1e-3; relTol=1e-2;
%   >> exactInteg = besseli(0,1)^dim;
%   >> inputArgs = {'absTol',absTol,'relTol',relTol};
%   >> [~,muhat]=cubBayesNet_g(fun, dim, inputArgs{:});
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%
%
%  Please refer to dt_cubBayesNet_g for more examples
%
%  See also CUBBAYESLATTICE_G, CUBSOBOL_G, CUBLATTICE_G, CUBMC_G, MEANMC_G,
%  INTEGRAL_G
%
%
%  References
%
%   [1] Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using
%   matching kernels and designs," PhD thesis, Illinois Institute of
%   Technology, 2019.
%
%   [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%
%
classdef cubBayesNet_g < handle
  
  properties
    f = @(x) x.^2; %function to integrate
    dim = 1; %dimension of the integrand
    mmin = 8; %min number of points to start with = 2^mmin
    mmax = 20; %max number of points allowed = 2^mmax
    absTol = 0.01; %absolute tolerance
    relTol = 0; %relative tolerance
    order = 1; %order of the kernel
    digitalNetAlpha = 1; %order of the digital nets
    stopCriterion = 'MLE'; %Available options {'MLE', 'GCV', 'full'}
    
    alpha = 0.01; % p-value, default 0.1%.
    stopAtTol = true; %automatice mode: stop after meeting the error tolerance
    arbMean = true; %by default use zero mean algorithm
    fName = 'None'; %name of the integrand
    figSavePath = ''; %path where to save he figures
    visiblePlot = true; %make plots visible
    debugEnable = false; %enable debug prints
    gaussianCheckEnable = false; %enable plot to check Guassian pdf
    verify_ftilde = false; % enable to verrify iterative computation of ftilde
    net_type = 'Sobol'; % type of nets. Options: NX or Sobol
  end
  
  properties (SetAccess = private)
    avoidCancelError = true;
    
    gen_h = 0;
    gen_h_un = 0;
    
    mvec = [];
    fullBayes = false; % Full Bayes - assumes m and s^2 as hyperparameters,
    GCV = false; % Generalized cross validation
    uncert = 0;  % quantile value for the error bound
    
    % variables to save debug info in each iteration, indexed in sequence
    errorBdAll = [];  % computed error bound
    muhatAll = []; % computed muhat, i.e, estimated integral
    aMLEAll = []; % shape parameter
    lossMLEAll = []; % MLE objective function loss
    timeAll = []; % time to to finish
    dscAll = []; % discriminant
    s_All = []; % scale parameter
  end
  
  methods
    function [obj,muhat] = cubBayesNet_g(varargin)  %Constructor
      if nargin > 0
        iStart = 1;
        if isa(varargin{1},'cubBayesNet_g')
          obj = copy(varargin{1});
          iStart = 2;
        end
        if nargin >= iStart
          % parse and set the input arguments to obj
          obj.parse_input_args(nargin-iStart+1, varargin{:});
        end
      else
        obj.warn_fd();
      end
      
      if strcmp(obj.stopCriterion, 'full')
        % Full Bayes : The posterior error is a Student-t distribution
        obj.fullBayes = true;
      elseif strcmp(obj.stopCriterion, 'GCV')
        % use Generalized cross validation
        obj.GCV = true;
      else
        % use empirical Bayes : Maximum likelihood estimation
        obj.fullBayes = false;
        obj.GCV = false;
      end
      
      % Credible interval : two-sided confidence
      % i.e., 1-alpha percent quantile
      if obj.fullBayes
        % degrees of freedom = 2^mmin - 1
        obj.uncert = -tinv(obj.alpha/2, (2^obj.mmin) - 1);
      else
        obj.uncert = -norminv(obj.alpha/2);
      end
      
      obj.mvec = obj.mmin:obj.mmax;
      length_mvec = length(obj.mvec);
      
      obj.errorBdAll = zeros(length_mvec,1);
      obj.muhatAll = zeros(length_mvec,1);
      obj.aMLEAll = zeros(length_mvec,1);
      obj.lossMLEAll = zeros(length_mvec,1);
      obj.timeAll = zeros(length_mvec,1);
      obj.dscAll = zeros(length_mvec,1);
      obj.s_All = zeros(length_mvec,1);
      
      if obj.verify_ftilde==true || obj.debugEnable==true
        warning('Caution: debugEnable is set, this will increase computation time !')
      end
      
      if nargout > 1
        muhat = compInteg(obj);
      end
    end
    
    % parses each input argument passed and assigns to obj.
    % Warns any unsupported options passed
    function parse_input_args(obj, nargin, varargin)
      if nargin > 0
        iStart = 1;
        
        % need atleast two input arguments
        if nargin < 2
          obj.warn_fd();
        else
          if isa(varargin{iStart}, 'function_handle') && isnumeric(varargin{1+iStart})
            obj.f = varargin{iStart};  % must be a function_handle
            obj.dim = varargin{1+iStart};  % must be an int
          else
            % input passing style not supported
            error('GAIL:cubBayesNet_g:input_invalid',...
              'Invalid input aruguments.\n');
          end
        end
        
        if nargin == 4 && isa(varargin{2+iStart}, 'float') && isa(varargin{3+iStart}, 'float')
          obj.absTol = varargin{2+iStart};
          obj.relTol = varargin{3+iStart};
        else
          
          if nargin == 3
            inParams = varargin{2+iStart};
            if isa(inParams, 'struct')
              % convert struct to name value pairs
              inputArgs = [fieldnames(inParams) struct2cell(inParams)];
              inputArgs = reshape(inputArgs', 1, []);
              varargin = inputArgs;
            else
              % input passing style not supported
              error('GAIL:cubBayesNet_g:input_invalid',...
                'Invalid input aruguments.\n');
            end
          end
          
          % parse each input argument passed
          wh = find(strcmp(varargin(iStart:end),'absTol'));
          if ~isempty(wh), obj.absTol = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'relTol'));
          if ~isempty(wh), obj.relTol = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'order'));
          if ~isempty(wh), obj.order = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'arbMean'));
          if ~isempty(wh), obj.arbMean = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'stopAtTol'));
          if ~isempty(wh), obj.stopAtTol = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'figSavePath'));
          if ~isempty(wh), obj.figSavePath = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'fName'));
          if ~isempty(wh), obj.fName = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'alpha'));
          if ~isempty(wh), obj.alpha = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'stopCriterion'));
          if ~isempty(wh), obj.stopCriterion = varargin{wh+iStart}; end
          
        end
      end
      
      function validate_input_args(obj)
        if ~gail.isfcn(obj.f)
          warning('GAIL:cubBayesNet_g:fnotfcn',...
            'The given input f should be a function handle.\n' );
        end
        
        if ~(obj.order==1)
          warning('GAIL:cubBayesNet_g:r_invalid',...
            'Kernel order, r=%d, is not supported; it must be 1. The algorithm is using default value r=1.\n', ...
            obj.order);
          obj.order = 1;
        end
        
        if obj.dim>20
          error('GAIL:cubBayesNet_g:dim_invalid',...
            'Integrand dimension=%d, is not implemented; max allowed is 20.\n', ...
            obj.dim);
        end
        
        stopCriterions = {'full','GCV','MLE'};
        if ~ismember(obj.stopCriterion, stopCriterions)
          str_stopCriterions = strjoin(stopCriterions, ', ');
          warning('GAIL:cubBayesNet_g:stop_crit_invalid',...
            'Stop criterion = "%s" is not supported; it must be from "%s". The algorithm is using default value "MLE".\n', ...
            obj.stopCriterion, str_stopCriterions);
          obj.stopCriterion = 'MLE';
        end
        
        if obj.absTol <= 0
          warning('GAIL:cubBayesNet_g:absTol_invalid',...
            'absTol %f is invalid, must be positive and non-zero, changing to default value 0.01', obj.absTol)
          obj.absTol = 0.01;
        end
        
        if obj.relTol < 0
          warning('GAIL:cubBayesNet_g:relTol_invalid',...
            'relTol %f is invalid, must be positive, changing to default value 0', obj.relTol)
          obj.relTol = 0;
        end
        
      end
      
      if ~license('test', 'Signal_Toolbox')
        error('GAIL:cubBayesNet_g:fwht_not_supported',...
          'Signal Processing Toolbox license required to continue !\n');
      end
      
      validate_input_args(obj);
    end
    
    % computes the integral
    function [muhat,out] = compInteg(obj)
      % comment this line of code to use GPU for computations
      gpuArray = @(x) x;   gather = @(x) x;
      obj = gpuArray(obj);
      successFlag = false;
      tstart = tic; %start the clock
      
      numM = length(obj.mvec);
      
      % Initialize the points generator
      gen_digital_nets(obj,true);
      
      % no periodization transfrom required for this algorithm
      ff = obj.f;
      
      %% Iteratively find the number of points required for the cubature to meet
      % the error threshold
      for iter = 1:numM
        tstart_iter = tic;
        m = obj.mvec(iter);
        n = 2^m;
        
        if iter == 1
          % xpts = sobstr(n0:nnext,1:obj.dim);  % grab Sobol' points
          % xpts = digitalseq_b2g(obj.dim, nnext-n0+1)';
          [xpts,xpts_un] = gen_digital_nets(obj,false,1,n);
          
          fx = gpuArray(ff(xpts)); %evaluate integrand
          ftilde = cubBayesNet_g.fwht_hs(fx);
        else
          % xptsnext = sobstr(n0:nnext,1:obj.dim);
          % xptsnext = digitalseq_b2g(obj.dim, nnext-n0+1)';
          [xptsnext,xptsnext_un] = gen_digital_nets(obj,false,1+n/2,n);
          
          xpts = [xpts;xptsnext];
          xpts_un = [xpts_un;xptsnext_un];
          
          fx = gpuArray(ff(xptsnext));  % initialize for inplace computation
          ftildeNext = cubBayesNet_g.fwht_hs(fx);
          ftilde = [(ftilde+ftildeNext)/2; (ftilde-ftildeNext)/2];
          % ftilde = [(ftilde+ftildeNext); (ftilde-ftildeNext)];
          
          if obj.verify_ftilde == true
            fx = gpuArray(ff(xpts)); %evaluate integrand
            ftilde_direct = cubBayesNet_g.fwht_hs(fx);
            if std(abs(ftilde-ftilde_direct)) > 0.1
              fprintf('fwht_hs computation is wrong');
            end
          end
        end
        
        
        %br_xpts = bitrevorder(gpuArray(xpts));
        br_xpts = gpuArray(xpts_un);
        
        %Compute MLE parameter
        lnaRange = [-8,5];
        lnaMLE = fminbnd(@(lna) ...
          ObjectiveFunction(obj, exp(lna),br_xpts,ftilde), ...
          lnaRange(1),lnaRange(2),optimset('TolX',1e-2));
        
        thetaOpt = exp(lnaMLE);
        [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj, thetaOpt,br_xpts,ftilde);
        
        %Check error criterion
        
        % out.ErrBd = 2.58*sqrt(DSC*RKHSnorm/n);
        if obj.fullBayes==true
          % full bayes
          if obj.avoidCancelError
            DSC = abs(Lambda_ring(1)/1);
          else
            DSC = abs((Lambda(1)/1) - 1);
          end
          % 1-alpha two sided confidence interval
          ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/(n-1));
        elseif obj.GCV==true
          % GCV based stopping criterion
          if obj.avoidCancelError
            DSC = abs(Lambda_ring(1)/(1 + Lambda_ring(1)));
          else
            DSC = abs(1 - (1/Lambda(1)));
          end
          temp = Lambda;
          temp(1) = n*(1 + Lambda_ring(1));
          C_Inv_trace = sum(1./temp(temp~=0));
          ErrBd = obj.uncert*sqrt(DSC * (RKHSnorm) /C_Inv_trace);
        else
          % empirical bayes
          if obj.avoidCancelError
            DSC = abs(Lambda_ring(1)/(1 + Lambda_ring(1)));
          else
            DSC = abs(1 - (1/Lambda(1)));
          end
          ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/n);
        end
        out.ErrBd = ErrBd;
        
        % store the debug information
        obj.dscAll(iter) = sqrt(DSC);
        obj.s_All(iter) = sqrt(RKHSnorm/n);
        
        if obj.arbMean == true % zero mean case
          % muhat is just the sample mean
          muhat = ftilde(1);
          
        else % non zero mean case
          muhat = ftilde(1)/Lambda(1);
        end
        muminus = muhat - out.ErrBd;
        muplus = muhat + out.ErrBd;
        obj.muhatAll(iter) = muhat;
        obj.errorBdAll(iter) = out.ErrBd;
        obj.aMLEAll(iter) = thetaOpt;
        obj.lossMLEAll(iter) = loss;
        
        if 2*out.ErrBd <= ...
            max(obj.absTol,obj.relTol*abs(muminus)) + max(obj.absTol,obj.relTol*abs(muplus))
          % if stopAtTol true, exit the loop
          % else, run for for all 'n' values.
          % Used to compute errorValues for 'n' vs error plotting
          if obj.stopAtTol==true
            successFlag = true;
            break
          end
        end
        
        obj.timeAll(iter) = toc(tstart_iter);
        
        % fprintf('thetaOpt=%1.3f, n=%d, ErrBd=%1.2e, time=%1.2f, absTol=%1.1e\n', ...
        %  thetaOpt,n,out.ErrBd,obj.timeAll(iter),obj.absTol);
      end %iteration loop
      
      out.n = n;
      out.time = toc(tstart);
      
      optParams.ErrBdAll = obj.errorBdAll;
      optParams.muhatAll = obj.muhatAll;
      optParams.mvec = obj.mvec;
      optParams.aMLEAll = obj.aMLEAll;
      optParams.timeAll = obj.timeAll;
      optParams.s_All = obj.s_All;
      optParams.dscAll = obj.dscAll;
      optParams.absTol = obj.absTol;
      optParams.relTol = obj.relTol;
      optParams.stopAtTol = obj.stopAtTol;
      out.optParams = optParams;
      
      if obj.stopAtTol==true
        % fprintf('Success n %d err_n %1.2e time %1.3f where absTol %1.2e relTol %1.2e\n', ...
        %    n, out.ErrBd, out.time, obj.absTol, obj.relTol)
      end
      
      if successFlag==true
        out.exitflag = 1;
      else
        out.exitflag = 2;  % error tolerance may not be met
      end
      
      if successFlag==false
        warning('GAIL:cubBayesNet_g:maxreached',...
          ['In order to achieve the guaranteed accuracy, ', ...
          sprintf('used maximum allowed sample size %d. \n', n)] );
      end
      % convert from gpu memory to local
      muhat = gather(muhat);
      out = gather(out);
      
    end  %function
    
    
    % generates higher order interlaced digital nets or sobol points
    function [xpts,xpts_un] = gen_digital_nets(obj,init,nStart,nEnd)
      if init==true
        if obj.digitalNetAlpha==1
          obj.gen_h = scramble(sobolset(obj.dim),'MatousekAffineOwen'); %generate a Sobol' sequence
          obj.gen_h_un = sobolset(obj.dim); %generate a Sobol' sequence
        else
          if strcmp(obj.net_type,'NX')
            switch obj.digitalNetAlpha
              case 2
                gen_name = 'nxmats\nx_s5_alpha2_m32.col';
              case 3
                gen_name = 'nxmats\nx_s5_alpha3_m32.col';
              otherwise
                error('Higher order=%d NX net not supported', obj.digitalNetAlpha)
            end
          else
            switch obj.digitalNetAlpha
              case 2
                gen_name = 'sobolmats\sobol_alpha2_Bs53.col';
              case 3
                gen_name = 'sobolmats\sobol_alpha3_Bs53.col';
              otherwise
                error('Higher order=%d sobol net not supported', obj.digitalNetAlpha)
            end
          end
          genmat = load(['dirk_nuyens\qmc-generators\DIGSEQ\' gen_name]);
          digitalseq_b2g('init0', genmat)
        end
      else
        if obj.digitalNetAlpha==1
          xpts = obj.gen_h(nStart:nEnd,1:obj.dim); %grab Sobol' points
          xpts_un = obj.gen_h_un(nStart:nEnd,1:obj.dim); %grab unscrambled Sobol' points
        else
          xpts_un = digitalseq_b2g(obj.dim, nEnd-nStart+1)';
          xpts = xpts_un;  % not yet implemented
        end
      end
    end
    
    
    % Objective function to find the optimal shape parmaeter
    function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,theta,xpts,ftilde)
      
      n = length(ftilde);
      % Compute the eigenvalues of the covariance matrix
      [Lambda,Lambda_ring] = cubBayesNet_g.kernel(xpts,obj.order,theta,obj.avoidCancelError);
      
      % Not required ftilde = abs(ftilde);  % remove any negative/imaginary values
      temp = ((ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)));
      
      % compute loss
      % Note: unlike FFT, FWT divides the output by 'n'
      if obj.GCV==true
        % GCV
        temp_gcv = abs(ftilde(Lambda~=0)./(Lambda(Lambda~=0))).^2 ;
        loss1 = 2*log(sum(1./Lambda(Lambda~=0)));
        loss2 = log(sum(temp_gcv(2:end)));
        % ignore all zero eigenvalues
        loss = loss2 - loss1;
        
        if obj.arbMean==true
          RKHSnorm = sum(temp_gcv(2:end));
        else
          RKHSnorm = sum(temp_gcv);
        end
      else
        % default: MLE
        if obj.arbMean==true
          RKHSnorm = sum(temp(2:end));  % already divided by 'n'
          temp_1 = sum(temp(2:end));
        else
          RKHSnorm = sum(temp);  % already divided by 'n'
          temp_1 = sum(temp);
        end
        if obj.debugEnable==true
          cubBayesNet_g.alertMsg(temp_1, 'Imag');
        end
        
        % loss = mean(log(Lambda)) + log(RKHSnorm);
        loss1 = sum(log(abs(Lambda(Lambda~=0))));
        loss2 = n*log(abs(temp_1 + eps)); % add eps to avoid zero
        % ignore all zero val eigenvalues
        loss = loss1 + loss2;
      end
      
      if obj.debugEnable==true
        cubBayesNet_g.alertMsg(loss1, 'Inf');
        cubBayesNet_g.alertMsg(loss2, 'Inf');
        cubBayesNet_g.alertMsg(loss, 'Inf', 'Imag', 'Nan');
        cubBayesNet_g.alertMsg(Lambda, 'Imag');
      end
    end
    
  end % end of methods
  
  % static methods are the ones that do not need to access obj
  methods(Static)
    % compute fast walsh transform in 'hadamard' ordering.
    % Matlab by default uses 'sequency' ordering, thus the need to be spcific
    function t = fwht_hs(fx)
      [n, ~] = size(fx);
      t = fwht(fx,n,'hadamard');
      
      % Note: Unlike fft, fwht normalizes the output, i.e. divideds by 'n'
    end
    
    % prints debug message if the given variable is Inf, Nan or
    % complex, etc
    function alertMsg(varargin)
      %varname = @(x) inputname(1);
      if nargin > 1
        iStart = 1;
        varTocheck = varargin{iStart};
        iStart = iStart + 1;
        inpvarname = inputname(1);
        
        while iStart <= nargin
          type = varargin{iStart};
          iStart = iStart + 1;
          
          switch type
            case 'Nan'
              if isnan(varTocheck)
                warning('%s has NaN values', inpvarname);
              end
            case 'Inf'
              if isinf(varTocheck)
                warning('%s has Inf values', inpvarname);
              end
            case 'Imag'
              if ~isreal(varTocheck)
                warning('%s has complex values', inpvarname)
              end
            otherwise
              warning('%s : unknown type check requested !', inpvarname)
          end
        end
      end
      
    end
    
    % Builds the difference matrix to compute kernel matrix, useful for
    % debugging or verification purpose
    % x : input points of size [n,d]
    function dm = diffMatrix(x)
      [n, d] = size(x);
      x_int = x*(2^64);  % convert to integer values
      A = reshape(x_int, n,1,d);
      dmA = repmat(A, [1,n,1]);  % stack into n columns
      B = reshape(x_int, 1,n,d);
      dmB = repmat(B, [n,1,1]);  % stack into n rows
      
      dm = bitxor(dmA, dmB)*(2^(-64));  % bitwise subtraction
    end
    
    
    % Computes modified kernel Km1 = K - 1
    % Useful to avoid cancellation error in the computation of (1 - n/\lambda_1)
    function [Km1, K] = kernel_t(aconst, Bern)
      d = size(Bern, 2);
      if length(aconst) == 1
        theta = ones(1,d)*aconst;
      else
        theta = aconst;  % Allow theta vary per dimension
      end
      
      Kjm1 = theta(1)*Bern(:,1);  % Kernel at j-dim minus One
      Kj = 1 + Kjm1;  % Kernel at j-dim
      
      for j=2:d
        Kjm1_prev = Kjm1; Kj_prev = Kj;  % save the Kernel at the prev dim
        
        Kjm1 = theta(j)*Bern(:,j).*Kj_prev + Kjm1_prev;
        Kj = 1 + Kjm1;
      end
      
      Km1 = Kjm1; K = Kj;
    end
    
    % Interface for all the kernels and to call cancellation error fix
    function [Lambda,Lambda_ring] = kernel(xpts,order,theta,avoidCancelError)
      kernelFunc = cubBayesNet_g.BuildKernelFunc(order);
      if avoidCancelError
        [C1m1] = cubBayesNet_g.kernel_t(theta, kernelFunc(xpts));
        % eigenvalues must be real : Symmetric pos definite Kernel
        Lambda_ring = real(cubBayesNet_g.fwht_hs(C1m1));
        
        Lambda = Lambda_ring;
        Lambda(1) = Lambda_ring(1) + 1;  %n;
      else
        C1 = prod(1 + theta*kernelFunc(xpts),2);
        Lambda = real(cubBayesNet_g.fwht_hs(C1));
        Lambda_ring = Lambda;
        Lambda_ring(1) = Lambda_ring(1) - 1;  %n;
      end
      
      if any(Lambda < 0)
        % warning('eigen values are negative!\n')
        Lambda_ring(2:end) = abs(Lambda_ring(2:end));
        Lambda(2:end) = Lambda_ring(2:end);
      end
      
      if 0 % Create the full kernel matrix to compare
        dm = cubBayesNet_g.diffMatrix(xpts);
        [n,dim] = size(xpts);
        K = 1 + theta*kernelFunc(dm);
        if dim > 1
          K = prod(K, 3);  % for d>1 reduce the 3rd dim by multiplication
        end
        lambda = eig(K)/n;
        temp = lambda./sort(Lambda);
        if std(abs(temp)) > 0.1 || max(abs(temp-1)) > 0.1
          warning('estimated eigenvalues wrong')
          figure; loglog(sort(abs(temp),'descend'))
        end
      end
      
    end
    
    % Builds High order walsh kernel function
    function [kernFunc] = BuildKernelFunc(order)
      
      %a1 = @(x)(-floor(log2(x)));
      function out = a1(x)
        out = -floor(log2(x));
        out(x==0) = 0;  % a1 is zero when x is zero
      end
      
      %t1 = @(x)(2.^(-a1(x)));
      function out = t1(x)
        out = (2.^(-a1(x)));
        out(x==0) = 0;  % t1 is zero when x is zero
      end
      
      %t2 = @(x)(2.^(-2*a1(x)));
      function out = t2(x)
        out = (2.^(-2*a1(x)));
        out(x==0) = 0;  % t2 is zero when x is zero
      end
      
      s1 = @(x)(1-2*x);
      s2 = @(x)(1/3 - 2*(1-x).*x);
      ts2 = @(x)((1-5*t1(x))/2 - (a1(x)-2).*x);
      ts3 = @(x)((1-43*t2(x))/18 + (5*t1(x)-1).*x +(a1(x)-2).*x.^2);
      
      if order==1
        kernFunc = @(x)(6*( (1/6) - 2.^(floor(log2(x))-1) ));
      elseif order==2
        omega2_1D = @(x)(s1(x) + ts2(x));
        kernFunc = omega2_1D;
      elseif order==3
        omega3_1D = @(x)(s1(x) + s2(x) + ts3(x));
        kernFunc = omega3_1D;
      else
        error('kernel order not yet supported');
      end
      
    end
    
    
    % plots the objective function for the MLE of theta
    function [minTheta, hFigCost] = plotObjectiveFunc(obj)
      
      n = 2^obj.mmax;
      gen_digital_nets(obj,true); % initialize digital net sequence
      [xpts, xpts_un] = gen_digital_nets(obj,false,1,n); % grab points
      % use scrambled points only to compute f
      fx = obj.f(xpts);  % No periodization transform required
      numM = length(obj.mvec);
      
      %% plot MLEKernel cost function
      lnTheta = -8:0.2:5;
      
      costMLE = zeros(numM,numel(lnTheta));
      tstart=tic;
      
      for iter = 1:numM
        nii = 2^obj.mvec(iter);
        
        eigvalK = zeros(numel(lnTheta),nii);
        %br_xpts = bitrevorder(xpts(1:nii, :));
        %ftilde = cubBayesNet_g.fwht_hs(obj.f(br_xpts)); %/nii;
        ftilde = cubBayesNet_g.fwht_hs(fx(1:nii)); %, nii, 'hadamard'
        br_xpts = xpts_un(1:nii,:);
        
        tic
        %par
        for k=1:numel(lnTheta)
          [costMLE(iter,k),eigvalK(k,:)] = ObjectiveFunction(obj, exp(lnTheta(k)),...
            br_xpts,ftilde);
          %[costMLE(iter,k),eigvalK(k,:)] = MLEKernel(obj, exp(lnTheta(k)),...
          %  br_xpts,ftilde);
        end
        toc
      end
      
      toc(tstart)
      
      if obj.visiblePlot==false
        hFigCost = figure('visible','off');
      else
        hFigCost = figure();
      end
      
      if ~isreal(costMLE)
        fprintf('costMLE has complex values !! \n')
      end
      
      % semilogx : loglog
      semilogx(exp(lnTheta),real(costMLE));
      set(hFigCost, 'units', 'inches', 'Position', [0 0 13.5 11.5])
      xlabel('Shape param, \(\theta\)')
      ylabel('MLE Cost, \( \log \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
      % ylabel('Log MLE Obj. fun.')
      axis tight;
      if obj.arbMean
        mType = '\(m \neq 0\)';  % arb mean
      else
        mType = '\(m = 0\)';  % zero mean
      end
      %title(lgd,'Sample Size, \(n\)'); legend boxoff
      title(sprintf('%s d=%d r=%d %s', obj.fName, obj.dim, ...
        obj.order, mType));
      [minVal,Index] = min(real(costMLE),[],2);
      
      % mark the min theta values found using fminbnd
      minTheta = exp(lnTheta(Index));
      hold on;
      semilogx(minTheta,minVal, '.');
      if exist('aMLEAll', 'var')
        semilogx(obj.aMLEAll,obj.lossMLEAll, '+');
      end
      temp = string(obj.mvec); temp=strcat('\(2^{',temp,'}\)');
      temp(end+1) = '\(\theta_{min_{true}}\)';
      if exist('aMLEAll', 'var')
        temp(end+1) = '\(\theta_{min_{est}}\)';
      end
      
      legend(temp,'location','best'); axis tight
      % gail.save_image(hFigCost, 'Paper_cubBayesLattice_g', plotFileName)
    end
    
    
    
    % plots the kernel with different params
    % useful for visualization and debugging
    % For Ex:
    %   obj = cubBayesNet_g(); obj.demoKernel(1024,2,2,1)
    function demoKernel(npts,ndims,order,theta)
      
      sobstr=sobolset(ndims); %generate a Sobol' sequence
      xpts = sobstr(1:npts,1:ndims); %grab Sobol' points
      
      kernelFunc = cubBayesNet_g.BuildKernelFunc(order);
      C1 = prod(1 + theta*kernelFunc(xpts),2);
      % Lambda = cubBayesNet_g.kernel(xpts,order,theta,avoidCancelError);
      
      if ndims==1
        hFig = figure();
        set(hFig, 'units', 'inches', 'Position', [4 4 6.5 5.5])
        plot(xpts, C1, '.', 'MarkerSize', 10); grid on; axis([0 1 -1 3])
      elseif ndims==2
        ns = sqrt(npts);
        if floor(ns) ~= ns
          error('number of points should be n = 2^m');
        end
        x = sobstr(1:ns,1:1); %grab Sobol' points
        x = sort(x);
        [X,Y] = meshgrid(x,x);
        xpts = [X(:) Y(:)];
        C1 = prod(1 + theta*kernelFunc(xpts),2);
        Z = reshape(C1, ns,ns);
        figure(); surf(X,Y, Z); axis tight
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$\omega_2$')
      elseif ndims==2 && false
        figure(); plot3(xpts(:,1), xpts(:,2), K1,'.'); axis tight; grid on
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$\omega_2$')
        
      else
        error('demoKernel: ndims > 2 not implemented');
      end
    end
    function warn_fd()
      warning('GAIL:cubBayesNet_g:fdnotgiven',...
        'At least, function f and dimension need to be specified');
    end
  end % end of static functions
end