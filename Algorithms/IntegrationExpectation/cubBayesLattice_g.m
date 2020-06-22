%CUBBAYESLATTICE_G Bayesian cubature method to estimate the integral
% of a random variable using rank-1 Lattices over a d-dimensional region
% within a specified generalized error tolerance with guarantees under Bayesian
% assumptions.
%
%   [OBJ,Q] = CUBBAYESLATTICE_G(f, dim, 'absTol',absTol, ...
%         'relTol',relTol, 'order',order, ...
%         'ptransform',ptransform, 'arbMean',arbMean);
%   initializes the object with the given parameters and also returns an
%   estimate of integral Q.
%
%   [Q,OutP] = COMPINTEG(OBJ); estimates the integral of f over hyperbox
%   [0,1]^dim using rank-1 Lattice sampling to within a specified generalized
%   error tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <=
%   tolfun with confidence of at least 99%, where I is the true integral
%   value, Q is the estimated integral value, abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation; however, if | I | is
%   rather small, then abstol determines the accuracy of the estimation.
%   Given the construction of our Lattices, d must be a positive integer
%   with 1 <= dim <= 600. For higher dimensions, it is recommended to use 
%   simpler periodization transformation like 'Baker'.
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
%   [OBJ,Q] = CUBBAYESLATTICE_G(f,dim); estimates the integral of f over
%   hyperbox [0,1]^dim using rank-1 Lattice sampling. All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.
%
%   [OBJ,Q] = CUBBAYESLATTICE_G(f,dim,absTol,relTol); estimates the integral
%   of f over hyperbox [0,1]^dim using rank-1 Lattice sampling. All parameters
%   should be input in the order specified above. The answer is given within
%   the generalized error tolerance tolfun. All other input parameters
%   are initialized with default values as given below.
%
%   [OBJ,Q] = CUBBAYESLATTICE_G(f,dim,inParms); estimates the integral
%   of f over hyperbox [0,1]^dim using rank-1 Lattice sampling.
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
%     order --- order of the Bernoulli polynomial of the kernel r=1,2.
%             If r == 0, algorithm automatically chooses the kernel order
%             which can be a non-integer value.
%             Default is 2
%     ptransform --- periodization variable transform to use: 'Baker',
%                    'C0', 'C1', 'C1sin', or 'C2sin'. Default is 'C1sin'
%     arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean. Default is 'true'
%     alpha --- confidence level for a credible interval of Q. Default is 0.01
%     mmin --- min number of samples to start with: 2^mmin. Default is 10
%     mmax --- max number of samples allowed: 2^mmax. Default is 22
%     stopCriterion -- stopping criterion to use. Supports three options 
%                     1) MLE: Empirical Bayes, 
%                     2) GCV: Generalized Cross Validation
%                     3) full: Full Bayes
%                     Default is MLE: Empirical Bayes
%     useGradient -- If true uses gradient descent in parameter search.
%                   Default is false
%     oneTheta -- If true uses common shape parameter for all dimensions,
%                 else allow shape parameter vary across dimensions.
%                 Default is true
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
%   This algorithm attempts to calculate the integral of function f over the
%   hyperbox [0,1]^dim to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
%   with guaranteed confidence level, e.g., 99% when alpha=0.5%. If the
%   algorithm terminates without showing any warning messages and provides
%   an answer Q, then the following inequality would be satisfied:
%
%   Pr(| Q - I | <= tolfun) = 99%
%
%   Please refer to our paper [1] for detailed arguments and proofs.
%
%  Examples
%
%   Example 1:
%   If no parameters are parsed, help text will show up as follows:
%   >> help cubBayesLattice_g
%   ***Bayesian cubature method to estimate the integral ***
%
%
%   Example 2: Quadratic
%   Estimate the integral with integrand f(x) = x.^2 over the interval [0,1]
%   with default parameters: order=2, ptransform=C1sin, abstol=0.01, relTol=0
%
%   >> warning('off','GAIL:cubBayesLattice_g:fdnotgiven')
%   >> [~,muhat] = cubBayesLattice_g;
%   >> exactInteg = 1.0/3;
%   >> warning('on','GAIL:cubBayesLattice_g:fdnotgiven')
%   >> check = double(abs(exactInteg-muhat) < 0.01)
%   check = 1
%
%
%   Example 3: ExpCos
%   Shape parameter independently chosen for each dimension
%
%   Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over the
%   interval [0,1] with parameters: order=2, ptransform=C1sin, abstol=0.001,
%   relTol=0.01
%
%   >> fun = @(x) exp(sum(cos(2*pi*x), 2));
%   >> dim=2; absTol=1e-2; relTol=1e-2;
%   >> exactInteg = besseli(0,1)^dim;
%   >> inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
%   >> inputArgs = [{fun, dim, 'absTol',absTol,'oneTheta',false} inputArgs];
%   >> obj=cubBayesLattice_g(inputArgs{:});
%   >> [muhat,outParams]=compInteg(obj);
%   >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
%   check = 1
%   >> etaDim = size(outParams.optParams.aMLEAll, 2)
%   etaDim = 2
%
%
%  Please refer to dt_cubBayesLattice_g for more examples
%
%
%  See also CUBBAYESNET_G, CUBSOBOL_G, CUBLATTICE_G, CUBMC_G, MEANMC_G,
%  INTEGRAL_G
%
%  References
%
%   [1] Jagadeeswaran Rathinavel, Fred J. Hickernell, Fast automatic
%   Bayesian cubature using lattice sampling.  Stat Comput 29, 1215-1229
%   (2019). https://doi.org/10.1007/s11222-019-09895-9
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
classdef cubBayesLattice_g < handle
  properties
    f = @(x) x.^2; %function to integrate
    dim = 1; %dimension of the integrand
    absTol = 0.01; %absolute tolerance
    relTol = 0; %relative tolerance
    order = 2; %Bernoulli order of the kernel
    %If zero, choose order automatically
    alpha = 0.01; % p-value, default 0.1%.
    ptransform = 'C1sin'; %periodization transform
    stopAtTol = true; %automatic mode: stop after meeting the error tolerance
    arbMean = true; %by default use zero mean algorithm
    stopCriterion = 'MLE'; % Available options {'MLE', 'GCV', 'full'}
    mmin = 8; %min number of samples to start with = 2^mmin
    mmax = 22; %max number of samples allowed = 2^mmax
    useGradient = false; %If true uses gradient descent in parameter search
    oneTheta = true; %If true use common shape parameter for all dimensions
    %else allow shape parameter vary across dimensions
  end
  
  properties (SetAccess = private)
    ff = []; %integrand after the periodization transform
    mvec = []; % n = 2^m
    uncert = 0;  % quantile value for the error bound
    gen_vec = []; % generator for the Lattice points
    
    fullBayes = false; % Full Bayes - assumes m and s^2 as hyperparameters,
    GCV = false; % Generalized cross validation
    vdc_order = false; % use Lattice points generated in vdc order
    kernType = 1; % Type-1: Bernoulli polynomaial based algebraic convergence
    % Type-2: Truncated series
    
    % only for developers use
    fName = 'None'; %name of the integrand
    figSavePath = ''; %path where to save he figures
    visiblePlot = false; %make plots visible
    debugEnable = false; %enable debug prints
    gaussianCheckEnable = false; %enable plot to check Gaussian pdf
    avoidCancelError = true; % avoid cancellation error in stopping criterion
    
    % variables to save debug info in each iteration
    errorBdAll = [];
    muhatAll = [];
    aMLEAll = [];
    lossMLEAll = [];
    timeAll = [];
    dscAll = [];
    s_All = [];
  end
  
  methods
    function [obj,muhat] = cubBayesLattice_g(varargin)  %Constructor
      
      if nargin > 0
        iStart = 1;
        if isa(varargin{1},'cubBayesLattice_g')
          obj = copy(varargin{1});
          iStart = 2;
        end
        % parse each input argument passed
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
      
      % apply periodization transformation to the function
      obj.ff = obj.doPeriodTx(obj.f, obj.ptransform);
      obj.gen_vec = obj.get_lattice_gen_vec(obj.dim);
      
      obj.mvec = obj.mmin:obj.mmax;
      length_mvec = length(obj.mvec);
      
      % to store debug info
      temp = cubBayesLattice_g.gpuArray_(zeros(length_mvec,1));
      obj.errorBdAll = temp;
      obj.muhatAll = temp;
      if obj.oneTheta==true
        obj.aMLEAll = temp;
      else
        obj.aMLEAll = cubBayesLattice_g.gpuArray_(zeros(length_mvec,obj.dim));
      end
      obj.lossMLEAll = temp;
      obj.timeAll = temp;
      obj.dscAll = temp;
      obj.s_All = temp;
      
      if nargout > 1
        muhat = compInteg(obj);
      end
      
    end
    
    % computes the integral
    function [muhat,out] = compInteg(obj)
      
      tstart = tic; %start the clock
      numM = length(obj.mvec);
      
      % pick a random value to apply as shift
      shift = rand(1,obj.dim);
      
      xun_=[];xpts_=[];ftilde_=[]; % temporary storage between iterations
      %% Iteratively find the number of points required for the cubature to meet
      % the error threshold
      for iter = 1:numM
        tstart_iter = tic;
        m = obj.mvec(iter);
        n = 2^m;
        
        %Update function values
        if obj.vdc_order
          [ftilde_,xun_,xpts_] = iter_fft_vdc(obj,iter,shift,xun_,xpts_,ftilde_);
        else
          [ftilde_,xun_,xpts_] = iter_fft(obj,iter,shift,xun_,xpts_,ftilde_);
        end
        [stop_flag,muhat,order_] = stopping_criterion(obj, xun_, ftilde_, iter, m);
        
        obj.timeAll(iter) = toc(tstart_iter);  % store per iteration time
        
        % if stopAtTol true, exit the loop
        % else, run for for all 'n' values.
        % Used to compute error values for 'n' vs error plotting
        if obj.stopAtTol==true && stop_flag==true
          break
        end
      end
      out.n = n;
      out.time = toc(tstart);
      out.ErrBd = obj.errorBdAll(end);
      
      optParams.ErrBdAll = obj.errorBdAll;
      optParams.muhatAll = obj.muhatAll;
      optParams.mvec = obj.mvec;
      optParams.aMLEAll = obj.aMLEAll;
      optParams.timeAll = obj.timeAll;
      optParams.s_All = obj.s_All;
      optParams.dscAll = obj.dscAll;
      optParams.absTol = obj.absTol;
      optParams.relTol = obj.relTol;
      optParams.shift = shift;
      optParams.stopAtTol = obj.stopAtTol;
      optParams.r = order_;
      out.optParams = optParams;
      if stop_flag==true
        out.exitflag = 1;
      else
        out.exitflag = 2;  % error tolerance may not be met
      end
      
      if stop_flag==false
        warning('GAIL:cubBayesLattice_g:maxreached',...
          ['In order to achieve the guaranteed accuracy, ', ...
          sprintf('used maximum allowed sample size %d. \n', n)] );
      end
      
      % convert from gpu memory to local
      muhat=obj.gather_(muhat);
      out=obj.gather_(out);
      
    end
    
    %% Efficient FFT computation algorithm, avoids recomputing the full fft
    function [ftilde_,xun_,xpts_] = iter_fft(obj,iter,shift,xun,xpts,ftildePrev)
      m = obj.mvec(iter);
      n = 2^m;
      
      % In every iteration except the first one, "n" number_of_points is doubled,
      % but FFT is only computed for the newly added points.
      % Previously computed FFT is reused.
      if iter == 1
        % In the first iteration compute full FFT
        xun_ = mod(bsxfun(@times,(0:1/n:1-1/n)',obj.gen_vec),1);
        xpts_ = mod(bsxfun(@plus,xun_,shift),1);  % shifted
        
        % Compute initial FFT
        ftilde_ = fft(obj.gpuArray_(obj.ff(xpts_))); %evaluate integrand's fft
      else
        xunnew = mod(bsxfun(@times,(1/n:2/n:1-1/n)',obj.gen_vec),1);
        xnew = mod(bsxfun(@plus,xunnew,shift),1);
        [xun_, xpts_] = obj.merge_pts(xun, xunnew, xpts, xnew, n, obj.dim);
        mnext=m-1;
        
        % Compute FFT on next set of new points
        ftildeNextNew = fft(obj.gpuArray_(obj.ff(xnew)));
        if obj.debugEnable
          cubBayesLattice_g.alertMsg(ftildeNextNew, 'Nan', 'Inf');
        end
        
        % combine the previous batch and new batch to get FFT on all points
        ftilde_ = obj.merge_fft(ftildePrev, ftildeNextNew, mnext);
      end
    end
    
    % Lattice points are ordered in van der Corput sequence, so we cannot use
    % Matlab's built-in fft routine. We use a custom one instead.
    function [ftilde_,xun_,xpts_] = iter_fft_vdc(obj,iter,shift,xun,xpts,ftildePrev)
      m = obj.mvec(iter);
      n = 2^m;
      
      % In every iteration except the first one, "n" number_of_points is doubled,
      % but FFT is only computed for the newly added points.
      % Previously computed FFT is reused.
      if iter == 1
        % in the first iteration compute the full FFT
        [xpts_,xun_] = obj.simple_lattice_gen(n,obj.dim,shift,true);
        
        % Compute initial FFT
        ftilde_ = obj.fft_DIT(obj.gpuArray_(obj.ff(xpts_)),m); %evaluate integrand's fft
      else
        [xnew,xunnew] = obj.simple_lattice_gen(n,obj.dim,shift,false);
        mnext=m-1;
        
        % Compute FFT on next set of new points
        ftildeNextNew = obj.fft_DIT(obj.gpuArray_(obj.ff(xnew)),mnext);
        if obj.debugEnable
          cubBayesLattice_g.alertMsg(ftildeNextNew, 'Nan', 'Inf');
        end
        
        xpts_ = [xpts; xnew];
        temp = zeros(n,obj.dim);
        temp(1:2:n-1,:) = xun;
        temp(2:2:n,:) = xunnew;
        xun_ = temp;
        % combine the previous batch and new batch to get FFT on all points
        ftilde_ = obj.merge_fft(ftildePrev, ftildeNextNew, mnext);
      end
      if obj.debugEnable
        cubBayesLattice_g.alertMsg(ftilde_, 'Inf', 'Nan');
      end
    end
    
    % decides if the user-defined error threshold is met
    function [success,muhat,r] = stopping_criterion(obj, xpts, ftilde, iter, m)
      
      tstart=tic;
      n=2^m;
      success = false;
      r = obj.order;
      
      % search for optimal shape parameter
      if obj.kernType==2
        % Truncated series kernel
        options = optimset('TolX',1e-2, 'PlotFcns',@optimplotfval);
        fLoss = @(param)ObjectiveFunctionFmin(obj,exp(param(1)),...
          1 + exp(param(2)),xpts,ftilde);
        % theta0 = ones(1,obj.dim+1)*0.1;
        theta0 = [0.1,0.1];
        
        [aOPT] = fminsearch(fLoss, ...
          theta0,optimset('TolX',1e-2));
        
        thetaOpt = exp(aOPT(1));
        bOpt = 1 + exp(aOPT(2));
        r = bOpt;
        [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunctionFmin(obj, ...
          thetaOpt,bOpt,xpts,ftilde);
        
      elseif obj.oneTheta==false
        % Use d length shape parameters, one per dimension
        if obj.useGradient==true
          % using Matlab fminunc
          theta0 = ones(1,obj.dim)*0.1;
          options = optimoptions('fminunc',...
            'Display','off',...
            'TolX',1e-2, ...
            'Algorithm','trust-region',...
            'SpecifyObjectiveGradient',true);
          try
            [lnThetaOpt] = fminunc(...
              @(phi)dObjectiveFunction(obj, exp(phi),xpts,ftilde), ...
              theta0,options);
          catch mErr
            error('*** cubBayesLattice_g: Error when evaluating fminunc.');
          end
          thetaOpt = exp(lnThetaOpt);
        else
          % Nelder Mead: gradient not required
          nm_start = tic;
          theta0 = zeros(1,obj.dim)*0.1; %
          nmOptions = optimset('TolX',1e-2);
          [lnaMLE_] = fminsearch(@(lna) ...
            ObjectiveFunction(obj, exp(lna),xpts,ftilde), ...
            theta0,nmOptions);
          thetaOpt = exp(lnaMLE_);
          time_nm = toc(nm_start);
        end
        [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj, ...
          thetaOpt,xpts,ftilde);
      else
        % Use one shape parameter for all dimensions
        % bounded search
        lnaRange = [-5,5];
        lnaMLE = fminbnd(@(lna) ...
          ObjectiveFunction(obj, exp(lna),xpts,ftilde), ...
          lnaRange(1),lnaRange(2),optimset('TolX',1e-2));
        thetaOpt = exp(lnaMLE);
        [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj, ...
          thetaOpt,xpts,ftilde);
      end
      
      %Check error criterion
      % compute DSC
      if obj.fullBayes==true
        % full bayes
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/n);
        else
          DSC = abs((Lambda(1)/n) - 1);
        end
        % 1-alpha two sided confidence interval
        ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/(n-1));
      elseif obj.GCV==true
        % GCV based stopping criterion
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/(n + Lambda_ring(1)));
        else
          DSC = abs(1 - (n/Lambda(1)));
        end
        temp = Lambda;
        temp(1) = n + Lambda_ring(1);
        C_Inv_trace = sum(1./temp(temp~=0));
        ErrBd = obj.uncert*sqrt(DSC * (RKHSnorm) /C_Inv_trace);
      else
        % empirical bayes
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/(n + Lambda_ring(1)));
        else
          DSC = abs(1 - (n/Lambda(1)));
        end
        ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/n);
      end
      
      if obj.arbMean==true % zero mean case
        muhat = ftilde(1)/n;
      else % non zero mean case
        muhat = ftilde(1)/Lambda(1);
      end
      muminus = muhat - ErrBd;
      muplus = muhat + ErrBd;
      
      % store intermediate values for post analysis
      % store the debug information
      obj.dscAll(iter) = sqrt(DSC);
      obj.s_All(iter) = sqrt(RKHSnorm/n);
      obj.muhatAll(iter) = muhat;
      obj.errorBdAll(iter) = ErrBd;
      obj.aMLEAll(iter, :) = thetaOpt;
      obj.lossMLEAll(iter) = loss;
      
      if obj.gaussianCheckEnable == true
        % plots the transformed and scaled integrand values as normal plot
        % Useful to verify the assumption, integrand was an instance of a Gaussian process
        CheckGaussianDensity(obj, ftilde, Lambda)
      end
      
      if 2*ErrBd <= ...
          max(obj.absTol,obj.relTol*abs(muminus)) + max(obj.absTol,obj.relTol*abs(muplus))
        if obj.errorBdAll(iter)==0
          obj.errorBdAll(iter) = eps;
        end
        % stopping criterion achieved
        success = true;
      end
      
    end
    
    % objective function to estimate shape parameter eta and order r
    % MLE : Maximum likelihood estimation
    % GCV : Generalized cross validation
    function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunctionFmin(obj,...
        a,b,xun,ftilde)
      
      if length(a) > 1
        if size(a,1) > size(a,2)
          a = a';  % we need row vector
        end
      end
      n = length(ftilde);
      [Lambda,Lambda_ring] = obj.dKernel(xun,b,a,obj.avoidCancelError, ...
        obj.kernType,false,obj.debugEnable);
      
      % compute RKHSnorm
      temp = abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)) ;
      
      % compute loss
      if obj.GCV==true
        % GCV
        temp_gcv = abs(ftilde(Lambda~=0)./(Lambda(Lambda~=0))).^2 ;
        loss1 = 2*log(sum(1./Lambda(Lambda~=0)));
        loss2 = log(sum(temp_gcv(2:end)));
        % ignore all zero eigenvalues
        loss = loss2 - loss1;
        
        if obj.arbMean==true
          RKHSnorm = sum(temp_gcv(2:end))/n;
        else
          RKHSnorm = sum(temp_gcv)/n;
        end
      else
        % default: MLE
        if obj.arbMean==true
          RKHSnorm = sum(temp(2:end))/n;
          temp_1 = sum(temp(2:end));
        else
          RKHSnorm = sum(temp)/n;
          temp_1 = sum(temp);
        end
        
        % ignore all zero eigenvalues
        loss1 = sum(log(abs(Lambda(Lambda~=0))));
        loss2 = n*log(abs(temp_1 + eps)); % add eps to avoid zero
        loss = loss1 + loss2;
      end
      
      if obj.debugEnable
        cubBayesLattice_g.alertMsg(RKHSnorm, 'Imag');
        cubBayesLattice_g.alertMsg(loss1, 'Inf');
        cubBayesLattice_g.alertMsg(loss2, 'Inf');
        cubBayesLattice_g.alertMsg(loss, 'Inf', 'Imag', 'Nan');
        cubBayesLattice_g.alertMsg(Lambda, 'Imag');
      end
    end
    
    
    function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,a,xun,ftilde)
      [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunctionFmin(obj,a,...
        obj.order,xun,ftilde);
    end
    
    
    % computes Gradient of the objective function
    function [lossObj,lossD] = dObjectiveFunction(obj,a,xun,ftilde)
      
      if length(a) > 1
        if size(a,1) > size(a,2)
          a = a';  % we need row vector
        end
      end
      n = length(ftilde);
      [dLambda, Lambda] = obj.dKernel(xun,obj.order,a,obj.avoidCancelError, ...
        obj.kernType,true,obj.debugEnable);
      
      % ignore all zero eigenvalues
      if obj.GCV==true
        % GCV
        temp_gcv = abs(ftilde(Lambda~=0)./(Lambda(Lambda~=0))).^2 ;
        deriv_part_num = abs(...
          bsxfun(@times, ...
          (temp_gcv)./(Lambda(Lambda~=0)), dLambda(Lambda~=0, :)));
        % (ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0).^3), dLambda(Lambda~=0, :)));
        if obj.arbMean==true
          deriv_part_num = sum(deriv_part_num(2:end, :), 1);
          deriv_part_den = sum(temp_gcv(2:end));
        else
          deriv_part_num = sum(deriv_part_num(1:end, :), 1);
          deriv_part_den = sum(temp_gcv);
        end
        deriv_part = deriv_part_num/deriv_part_den;
        det_part_num = sum(...
          bsxfun(@times, ...
          dLambda(Lambda~=0, :), (Lambda(Lambda~=0).^2)), 1);
        det_part = det_part_num/sum(1./Lambda(Lambda~=0));
        lossD = -2*deriv_part + 2*det_part;
        lossGCV = log(deriv_part_den) - 2*log(sum(1./Lambda(Lambda~=0)));
        lossObj = lossGCV;
      else
        % trace_part = sum(dLambda(Lambda~=0)./Lambda(Lambda~=0))/n;
        trace_part = sum(bsxfun(@rdivide, dLambda(Lambda~=0,:), Lambda(Lambda~=0)), 1)/n;
        
        % default: MLE
        temp = abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)) ;
        
        deriv_part_num = abs(...
          bsxfun(@times, ...
          (abs(ftilde(Lambda~=0))./Lambda(Lambda~=0)).^2, dLambda(Lambda~=0, :)));

        if obj.arbMean==true
          deriv_part_num = sum(deriv_part_num(2:end, :), 1);
          deriv_part_den = sum(temp(2:end));
        else
          deriv_part_num = sum(deriv_part_num(1:end, :), 1);
          deriv_part_den = sum(temp);
        end
        
        deriv_part = deriv_part_num/deriv_part_den;
        lossD = trace_part - deriv_part;
        lossDet = sum(log(abs(Lambda(Lambda~=0))))/n;
        lossMLE = lossDet + log(abs(deriv_part_den));
        lossObj = lossMLE;
        
        if obj.debugEnable
          cubBayesLattice_g.alertMsg(lossD, 'Inf', 'Imag', 'Nan');
          cubBayesLattice_g.alertMsg(lossObj, 'Inf', 'Imag', 'Nan');
        end
      end
      
    end
    
    
    % parses each input argument passed and assigns to obj.
    % Warns any unsupported options passed
    function parse_input_args(obj, nargin, varargin)
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
          error('GAIL:cubBayesLattice_g:input_invalid',...
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
            error('GAIL:cubBayesLattice_g:input_invalid',...
              'Invalid input aruguments.\n');
          end
        end
        
        % parse additional input arguments
        wh = find(strcmp(varargin(iStart:end),'absTol'));
        if ~isempty(wh), obj.absTol = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'relTol'));
        if ~isempty(wh), obj.relTol = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'order'));
        if ~isempty(wh), obj.order = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'ptransform'));
        if ~isempty(wh), obj.ptransform = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'arbMean'));
        if ~isempty(wh), obj.arbMean = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'stopAtTol'));
        if ~isempty(wh), obj.stopAtTol = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'figSavePath'));
        if ~isempty(wh), obj.figSavePath = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'fName'));
        if ~isempty(wh), obj.fName = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'stopCriterion'));
        if ~isempty(wh), obj.stopCriterion = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'alpha'));
        if ~isempty(wh), obj.alpha = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'useGradient'));
        if ~isempty(wh), obj.useGradient = varargin{wh+iStart}; end
        wh = find(strcmp(varargin(iStart:end),'oneTheta'));
        if ~isempty(wh), obj.oneTheta = varargin{wh+iStart}; end
        
      end
      
      if obj.order <= 0
        % find optimal kernel order
        obj.kernType = 2;
      end
      
      function validate_input_args(obj)
        if ~gail.isfcn(obj.f)
          warning('GAIL:cubBayesLattice_g:fnotfcn',...
            'The given input f should be a function handle.\n' );
        end
        
        if obj.dim>22
          error('GAIL:cubBayesLattice_g:dim_invalid',...
            'Integrand dimension=%d, is not implemented; max allowed is 22.\n', ...
            obj.dim);
        end
        
        if obj.kernType==1
          if ~(obj.order==0 || obj.order==1 || obj.order==2)
            warning('GAIL:cubBayesLattice_g:r_invalid',...
              'Kernel order, r=%d, is not supported; it must be 1 or 2. The algorithm is using default value r=1.\n', ...
              obj.order);
            obj.order = 1;
          end
        else
          if obj.oneTheta==false
            warning('GAIL:cubBayesLattice_g:r_invalid',...
              'Kernel order, r, must be specified when choosing oneTheta = false. The algorithm is using default value r=1.\n');
          end
          obj.order = 1;
        end
        
        stopCriterions = {'full','GCV','MLE'};
        if ~ismember(obj.stopCriterion, stopCriterions)
          str_stopCriterions = strjoin(stopCriterions, ',');
          warning('GAIL:cubBayesLattice_g:stop_crit_invalid',...
            'Stop criterion = "%s" is not supported; it must be from "%s". The algorithm is using default value "MLE".\n', ...
            obj.stopCriterion, str_stopCriterions);
          obj.stopCriterion = 'MLE';
        end
        
        var_txs = {'Baker','C0','C1','C1sin','C2sin','C3sin','none'};
        if ~ismember(obj.ptransform, var_txs)
          str_var_txs = strjoin(var_txs, ',');
          warning('GAIL:cubBayesLattice_g:var_transform_invalid',...
            'Periodizing transform = "%s" is not supported; the value must be from "%s". The algorithm is using default value "Baker".\n', ...
            obj.ptransform, str_var_txs);
          obj.ptransform = 'Baker';
        end
        
        if obj.absTol <= 0
          warning('GAIL:cubBayesLattice_g:absTol_invalid',...
            'absTol %f is invalid, must be positive and non-zero, changing to default value 0.01', obj.absTol)
          obj.absTol = 0.01;
        end
        
        if obj.relTol < 0
          warning('GAIL:cubBayesLattice_g:relTol_invalid',...
            'relTol %f is invalid, must be positive, changing to default value 0', obj.relTol)
          obj.relTol = 0;
        end
      end
      
      validate_input_args(obj);
    end
    
    
  end
  
  methods(Static)
    % prints debug message if the given variable is Inf, Nan or
    % complex, etc
    % Example: alertMsg(x, 'Inf', 'Imag')
    %          prints if variable 'x' is either Infinite or Imaginary
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
              if any(isnan(varTocheck))
                warning('%s has NaN values', inpvarname);
              end
            case 'Inf'
              if any(isinf(varTocheck))
                warning('%s has Inf values', inpvarname);
              end
            case 'Imag'
              if ~all(isreal(varTocheck))
                warning('%s has complex values', inpvarname)
              end
            otherwise
              warning('%s : unknown type check requested !', inpvarname)
          end
        end
      end
      
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
    
    function g = truncated_series(N, r)
      tilde_g_0 = 0;
      m = 1:(-1 + N/2);
      tilde_g_h1 = N./abs(m).^(r);
      m = (N/2):(-1 + N);
      tilde_g_h2 = N./abs(N-m).^(r);
      tilde_g = [tilde_g_0 tilde_g_h1 tilde_g_h2];
      g = ifft(tilde_g)';
    end
    
    function c = truncated_series_kernel(x,r)
      n = size(x, 1);
      g = cubBayesLattice_g.truncated_series(n,r);
      c = g(1 + x*n);
    end
    
    
    function [kernelFunc,constMult] = bernPolyFun(r)
      constMult = -(-1)^(r/2)*((2*pi)^r)/factorial(r);
      % constMult = -(-1)^(b_order/2);
      if r == 2
        bernPoly = @(x)(-x.*(1-x) + 1/6);
      elseif r == 4
        bernPoly = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
      else
        error('Bernoulli order=%d not implemented !', r);
      end
      kernelFunc = @(x) bernPoly(x);
    end
    
    % Computes gradient of the kernel w.r.t. shape parameter
    % dLambda : Eigenvalues of the gradient of the kernel
    function [Lambda,Lambda_ring,dLambda] = dKernel(xun,order,a,avoidCancelError,...
        kernType,gradient,debugEnable)
      
      if kernType==1
        % Bernoulli polynomial order as per the equation
        [kernelFunc,constMult] = cubBayesLattice_g.bernPolyFun(order*2);
      else
        r = order*2;
        constMult = 1;
        kernelFunc = @(x) cubBayesLattice_g.truncated_series_kernel(x,r);
      end
      
      [~, dim] = size(xun);
      if avoidCancelError
        kernFuncValues = kernelFunc(xun);
        % Computes C1m1 = C1 - 1
        % C1_new = 1 + C1m1 indirectly computed in the process
        [C1m1] = cubBayesLattice_g.kernel_t(a*constMult, kernFuncValues);
        % eigenvalues must be real : Symmetric pos definite Kernel
        
        Lambda_ring = real(fft(C1m1));
        if any(find(Lambda_ring < 0))
          Lambda_ring = abs(Lambda_ring);
        end
        
        % Lambda = real(fft(1 + C1m1));
        Lambda = Lambda_ring;
        Lambda(1) = Lambda_ring(1) + length(Lambda_ring);
        
        if gradient==true
          if length(a) > 1
            % different theta per dimension
            try
              partb = 1 - 1./(1+bsxfun(@times,a*constMult,kernFuncValues));
            catch mErr
              error('*** cubBayesLattice_g: Error when evaluating partb.');
            end
            % dC1 = (1./a)*(1 + C1m1).*partb;
            dC1 = bsxfun(@times, (1./a), ...
              bsxfun(@times, (1 + C1m1), partb));
          else
            % common theta
            partb = -(1/dim)*sum(1./(1+a*constMult*kernFuncValues), 2);
            dC1 = (dim/a)*(1 + C1m1 + partb + C1m1.*partb);
          end
          dLambda = real(fft(dC1));
        end
        if debugEnable == true
          % eigenvalues must be real : Symmetric pos definite Kernel
          Lambda_direct = real(fft(C1_alt)); % Note: fft output unnormalized
          if sum(abs(Lambda_direct-Lambda)) > 1
            fprintf('Possible error: check Lambda_ring computation')
          end
        end
      else
        % direct approach to compute first row of the kernel Gram matrix
        temp_ = bsxfun(@times, (a)*constMult, kernelFunc(xun));
        C1 = prod(1 + temp_, 2);
        
        % matlab's builtin fft is much faster and accurate
        % eigenvalues must be real : Symmetric pos definite Kernel
        Lambda = real(fft(C1));
        if gradient==true
          dC1 = (dim/a)*C1.*(1 - (1/dim)*sum(1./(1+temp_), 2));
          dLambda = real(fft(dC1));
        end
      end
      
    end
    
    
    % computes the periodization transform for the given function values
    function f = doPeriodTx(fInput, ptransform)
      
      if strcmp(ptransform,'Baker')
        f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
      elseif strcmp(ptransform,'C0')
        f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
      elseif strcmp(ptransform,'C1')
        % C^1 transform
        f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2);
      elseif strcmp(ptransform,'C1sin')
        % Sidi C^1 transform
        f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(2*sin(pi*x).^2,2);
      elseif strcmp(ptransform,'C2sin')
        % Sidi C^2 transform
        psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16;
        psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
        f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
      elseif strcmp(ptransform,'C3sin')
        % Sidi C^3 transform
        psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi);
        psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
        f=@(x) fInput(psi4(x)).*prod(psi4_1(x),2);
      elseif strcmp(ptransform,'none')
        % do nothing
        f=@(x) fInput(x);
      else
        error('Error: Periodization transform %s not implemented', ptransform);
      end
      
    end
    
    % just returns the generator for rank-1 Lattice point generation
    function [z] = get_lattice_gen_vec(d)
      % 600 dimensional 2^20 points generating vector from Dirks website exod2_base2_m20.txt
      z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599, 151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, 85729, 14597, 94813, 422013, 484367, 355029, 123065, 467905, 41129, 298607, 375981, 256421, 279695, 164795, 256413, 267543, 505211, 225547, 50293, 97031, 86633, 203383, 427981, 221421, 465833, 329843, 212325, 467017, 214065, 98063, 128867, 63891, 426443, 244641, 56441, 357107, 199459, 169327, 407687, 154961, 64579, 436713, 322855, 435589, 220821, 72219, 344125, 315189, 105979, 421183, 212659, 26699, 491987, 310515, 344337, 443019, 174213, 244609, 5979, 85677, 148663, 514069, 172383, 238589, 458305, 460201, 487365, 454835, 452035, 55005, 517221, 85841, 434641, 387469, 24883, 154373, 145103, 416491, 252109, 509385, 296473, 248789, 297219, 119711, 252395, 188293, 23943, 264817, 242005, 26689, 51931, 490263, 155451, 365301, 445277, 311581, 306887, 331445, 208941, 385313, 307593, 359113, 67919, 351803, 335955, 326111, 57853, 52153, 84863, 158013, 272483, 419143, 252581, 372097, 177007, 145815, 350453, 412791, 435559, 387627, 35887, 48461, 389563, 68569, 118715, 250699, 183713, 29615, 168429, 292527, 86465, 450915, 239063, 23051, 347131, 138885, 243505, 201835, 269831, 265457, 496089, 273459, 276803, 225507, 148131, 87909, 115693, 45749, 3233, 194661, 329135, 90215, 104003, 27611, 437589, 422687, 19029, 284433, 348413, 289359, 418785, 293911, 358343, 85919, 501439, 462941, 301185, 292875, 242667, 408165, 137921, 329199, 308125, 48743, 122291, 362643, 90781, 448407, 25389, 78793, 362423, 239423, 280833, 55483, 43757, 138415, 395119, 175965, 253391, 462987, 50655, 67155, 142149, 314277, 452523, 364029, 323001, 105873, 231785, 329547, 517581, 64375, 180745, 30693, 321739, 259327, 523313, 123863, 446629, 112611, 134019, 442879, 516621, 469677, 271077, 83859, 195209, 385581, 3287, 261841, 16525, 243831, 505215, 37669, 275001, 118849, 475943, 56509, 239489, 35893, 31015, 458209, 292255, 94197, 279055, 7573, 233705, 339587, 396313, 310037, 371939, 494279, 261481, 2875, 51129, 204067, 40633, 459101, 226639, 89795, 464665, 439937, 388665, 277539, 370801, 438367, 73733, 166153, 200849, 250477, 148655, 445817, 375723, 373433, 154819, 367247, 462549, 382217, 269073, 15985, 206263, 507895, 335263, 251183, 236851, 285491, 371291, 20143, 471543, 334263, 397501, 52335, 122837, 160981, 332741, 341961, 320455, 144133, 410489, 440261, 274789, 83793, 353867, 310001, 161271, 28267, 400007, 469779, 351385, 158419, 301117, 234521, 260047, 312511, 213851, 332001, 3699, 518163, 119209, 329387, 149889, 485193, 505407, 326067, 149541, 102343, 441707, 499551, 501199, 77817, 355999, 128165, 396261, 247463, 9733, 481107, 411379, 479917, 84085, 380091, 489765, 504237, 47847, 496129, 343905, 496621, 498123, 270835, 459931, 314289, 89077, 505051, 11647, 26765, 349111, 357217, 493937, 179089, 300189, 143621, 205639, 244475, 303281, 180189, 70443, 301471, 17853, 17121, 243179, 377849, 209079, 167565, 357373, 309503, 367039, 136041, 247861, 226573, 63631, 344345, 256401, 138305, 271675, 354845, 420971, 442981, 225321, 342755, 427957, 493767, 488177, 141063, 224621, 9439, 217623, 242451, 508557, 379609, 202291, 266555, 452509, 379789, 89867, 519873, 163115, 237191, 235291, 149683, 187821, 508801, 425951, 239141, 284505, 498919, 493857, 97373, 92147, 492967, 302591, 225277, 16947, 275043, 322807, 377713, 408445, 187103, 185133, 505963, 386109, 96301, 470963, 407939, 6601, 409277, 5031, 128747, 393271, 415197, 114049, 223999, 99373, 482183, 504981, 295837, 34235, 40765, 408397, 216741, 422925, 496079, 300813, 277283, 312489, 368009, 161369, 362997, 6663, 509953, 387903, 97597, 238917, 378851, 190545, 430029, 204931, 466553, 293441, 327939, 183495, 463331, 422655, 428099, 20715, 477503, 465937, 270399, 139589, 129581, 215571, 299645, 125221, 23345, 229345, 138059, 521769, 14731, 318159, 190173, 361381, 485577, 512807, 268009, 185937, 210939, 86965, 113005, 296923, 85753, 381527, 196325, 274565, 182689, 200951, 117371, 489747, 19521, 426587, 168393, 486039, 220941, 392473, 344051, 412275, 501127, 434941, 85569, 406757, 371643, 470783, 466117, 170707, 473019, 494155, 411809, 13371, 202745, 23597, 25621, 64351, 508445, 204947, 38279, 264269, 230499, 405605, 68513, 414481, 301849, 6815, 406425, 62881, 174349, 505503, 329037, 104357, 113815, 137669, 181689, 493057, 296191, 135279, 236891, 82135, 371269, 483993, 394407, 372929, 139823, 114515, 416815, 260309, 489593, 156763, 21523, 189285, 308129, 155369, 213557, 298023, 391439, 379245, 409109, 229765, 28521, 464087, 470911, 435965, 201451, 64371, 370499, 276377, 331635, 196813, 379415, 229547, 430067, 137053, 312839, 390385, 77155, 163911, 514381, 487453];

      z = z(1:d);
    end
    
    % generates rank-1 Lattice points in vander Corput sequence order
    function [xlat,xpts_un,xlat_un,xpts] = simple_lattice_gen(n,d,shift,firstBatch)
      z = cubBayesLattice_g.get_lattice_gen_vec(d);
      
      nmax = n;
      nmin = 1 + n/2;
      if firstBatch==true
        nmin = 1;
      end
      nelem=nmax-nmin+1;
      
      if firstBatch==true
        brIndices=cubBayesLattice_g.vdc(nelem)';
        xpts_un=mod(bsxfun(@times,(0:1/n:1-1/n)',z),1); % unshifted in direct order
      else
        brIndices=cubBayesLattice_g.vdc(nelem)'+1/(2*(nmin-1));
        xpts_un=mod(bsxfun(@times,(1/n:2/n:1-1/n)',z),1); % unshifted in direct order
        
      end
      xpts = mod(bsxfun(@plus,xpts_un,shift),1);  % shifted in direct order
      
      xlat_un = mod(bsxfun(@times,brIndices',z),1);  % unshifted
      xlat = mod(bsxfun(@plus,xlat_un,shift),1);  % shifted in VDC order
    end
    
    % van der Corput sequence in base 2
    function q = vdc(n)
      if n>1
        k=log2(n); % We compute the VDC seq part by part of power 2 size
        q=zeros(2^k,1);
        for l=0:k-1
          nl=2^l;
          kk=2^(k-l-1);
          ptind=repmat([false(nl,1);true(nl,1)],kk,1);
          q(ptind)=q(ptind)+1/2^(l+1);
        end
      else
        q=0;
      end
    end
    
    % fft with decimation in time i.e., input is already in 'bitrevorder'
    function y = fft_DIT(y,nmmin)
      for l=0:nmmin-1
        nl=2^l;
        nmminlm1=2^(nmmin-l-1);
        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
        coefv=repmat(coef,nmminlm1,1);
        evenval=y(ptind);
        oddval=y(~ptind);
        y(ptind)=(evenval+coefv.*oddval);
        y(~ptind)=(evenval-coefv.*oddval);
      end
    end
    
    % using FFT butefly plot technique merges two halves of fft
    function ftildeNew = merge_fft(ftildeNew, ftildeNextNew, mnext)
      ftildeNew=[ftildeNew;ftildeNextNew];
      nl=2^mnext;
      ptind=[true(nl,1); false(nl,1)];
      coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
      coefv=repmat(coef,1,1);
      evenval=ftildeNew(ptind);
      oddval=ftildeNew(~ptind);
      ftildeNew(ptind)=(evenval+coefv.*oddval);
      ftildeNew(~ptind)=(evenval-coefv.*oddval);
    end
    
    % inserts newly generated points with the old set by interleaving them
    % xun - unshifted points
    function [xun, x] = merge_pts(xun, xunnew, x, xnew, n, d)
      temp = zeros(n,d);
      temp(1:2:n-1,:) = xun;
      temp(2:2:n,:) = xunnew;
      xun = temp;
      temp(1:2:n-1,:) = x;
      temp(2:2:n,:) = xnew;
      x = temp;
    end
    
    % to enable GPU for computations
    function y = gpuArray_(x)
      % y = gpuArray(x);  % use GPU
      y = x;  % use CPU
    end
    
    function y = gather_(x)
      % y = gather(x);  % use GPU
      y = x;  % use CPU
    end
    
    function warn_fd()
      warning('GAIL:cubBayesLattice_g:fdnotgiven',...
        'At least, function f and dimension need to be specified');
    end
  end
end

