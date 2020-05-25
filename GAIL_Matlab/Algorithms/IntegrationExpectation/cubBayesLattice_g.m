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
%   [0,1]^d using rank-1 Lattice sampling to within a specified generalized
%   error tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <=
%   tolfun with confidence of at least 99%, where I is the true integral
%   value, Q is the estimated integral value, abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation; however, if | I | is
%   rather small, then abstol determines the accuracy of the estimation.
%   Given the construction of our Lattices, d must be a positive integer
%   with 1 <= d <= 600.
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
%   hyperbox [0,1]^d using rank-1 Lattice sampling. All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.
%
%   [OBJ,Q] = CUBBAYESLATTICE_G(f,dim,absTol,relTol); estimates the integral
%   of f over hyperbox [0,1]^d using rank-1 Lattice sampling. All parameters
%   should be input in the order specified above. The answer is given within
%   the generalized error tolerance tolfun. All other input parameters
%   are initialized with default values as given below.
%
%   [OBJ,Q] = CUBBAYESLATTICE_G(f,dim,inParms); estimates the integral
%   of f over hyperbox [0,1]^d using rank-1 Lattice sampling.
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
%               Default is 2
%     ptransform --- periodization variable transform to use: 'Baker',
%                    'C0', 'C1', 'C1sin', or 'C2sin'. Default is 'C1sin'
%     arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean. Default is 'true'
%     alpha --- confidence level for a credible interval of Q. Default is 0.01
%     mmin --- min number of samples to start with: 2^mmin. Default is 10
%     mmax --- max number of samples allowed: 2^mmax. Default is 22
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
      z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
        151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
        85729, 14597, 94813, 422013, 484367]; %generator
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

