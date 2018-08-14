%CUBBAYESLATTICE_G Bayesian cubature method to estimate the integral
% of a random variable
%
%   OBJ = CUBBAYESLATTICE_G('f',f,'dim',dim,'absTol',absTol,'relTol',relTol,...
%         'order',order, 'ptransform',ptransform, 'arbMean',arbMean);
%   Initializes the object with the given parameters.
%   Q = COMPINTEG(OBJ); estimates the integral of f over hyperbox [0,1]^d
%   using Rank-1 Lattice sampling to within a specified generalized error
%   tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <= tolfun
%   with cofidence of at least 99%, where I is the true integral value,
%   abstol is the absolute error tolerance, and reltol is the relative
%   error tolerance. Usually the reltol determines the accuracy of the
%   estimation, however, if | I | is rather small, then abstol determines
%   the accuracy of the estimation. Input f is a function handle that
%   accepts an n x d matrix input, where d is the dimension of the hyperbox,
%   and n is the number of points being evaluated simultaneously.
%
%   Input Arguments
%
%     f --- the integrand.
%     dim --- number of dimensions of the integrand.
%     absTol --- the absolute error tolerance | I - Q | <= absTol.
%     relTol --- the relative error tolerance | I - Q | <= I*relTol.

%   Optional Input Arguments
%
%     order --- order of the bernoulli polynomial of the kernel.
%     ptransform --- periodization transform to use
%     arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean
%
%  Guarantee
% This algorithm attempts to calculate the integral of function f over the
% hyperbox [0,1]^d to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level 99%. If the algorithm terminates
% without showing any warning messages and provides an answer Q, then the
% following inequality would be satisfied:
%
% Pr(| Q - I | <= tolfun) >= 99%
%
% Please refer to our paper for detailed arguments and proofs.
%
%  Examples
%
% Example 1:
%
% If no parameters are parsed, help text will show up as follows:
% >> cubBayesLattice_g
% ***Bayesian cubature method to estimate the integral ***
%
%
% Example 2: Quadratic
%
% Estimate the integral with integrand f(x) = x.^2 over the interval
% [0,1] with parameters: order=2, ptransform=Baker, abstol=0.01, relTol=0
%
% >> obj = cubBayesLattice_g;
% >> exactInteg = 1.0/3;
% >> muhat=compInteg(obj);
% >> check = double(abs(exactInteg-muhat) < 0.01)
% check = 1
%
% Example 3: ExpCos
%
% Estimate the integral with integrand f(x) = exp(sum(cos(2*pi*x)) over the
% interval [0,1] with parameters: order=4, ptransform=C1sin, abstol=0.01
%
% >> fun = @(x) exp(sum(cos(2*pi*x), 2));
% >> dim=2; absTol=1e-3; relTol=1e-2; fName = 'ExpCos';
% >> exactInteg = besseli(0,1)^dim;
% >> obj=cubBayesLattice_g('f',fun, 'dim',dim, 'absTol',absTol, 'relTol',relTol,...
% >>    'order',4, 'ptransform','C1sin');
% >> muhat=compInteg(obj);
% >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
% check = 1
%
% Example 3: Keister
%
% >> dim=2; absTol=1e-3; relTol=1e-2; fName = 'Keister';
% >> normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
% >> domain = repmat([0;1],[1,dim]);
% >> replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
% >> yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
% >> f1 = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
% >> fKeister = @(x) f1(x,dim); exactInteg = Keistertrue(dim);
% >> obj=cubBayesLattice_g('f',fKeister, 'dim',dim, 'absTol',absTol, 'relTol',relTol,...
% >>    'order',4, 'ptransform','C1','arbMean',false);
% >> muhat=compInteg(obj);
% >> check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
% check = 1
%
% Example 3: MVN
% >> dim=2; absTol=1e-3; relTol=1e-2; fName = 'MVN';
% >> C = [4 1 1; 0 1 0.5; 0 0 0.25]; Cov = C'*C;
% >> a = [-6 -2 -2]; b = [5 2 1];
% >> muBest = 0.676337324357787;
% >> MVNProbMLELattice = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvecMLE, ...
% >>   'errMeth','g','cubMeth','MLELattice','intMeth','Genz', ...
% >>   'BernPolyOrder',2,'ptransform','C1sin', ...
% >>   'fName',fName,'arbMean',true,'absTol',absTol,'relTol',relTol);
% >> muhat = compProb(MVNProbMLELattice);
% >> check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))
% check = 1
%


classdef cubBayesLattice_g < handle
  
  properties
    f = @(x) x.^2; %function to integrate
    dim = 1; %dimension of the integrand
    mmin = 10; %min number of points to start with = 2^mmin
    mmax = 22; %max number of points allowed = 2^mmax
    absTol = 0.01; %absolute tolerance
    relTol = 0; %relative tolerance
    order = 2; %order of the kernel
    alpha = 0.001; % the uncertainty, the default value is 0.1%.
    ptransform = 'Baker'; %periodization transform
    stopAtTol = true; %automatice mode: stop after meeting the error tolerance
    arbMean = true; %by default use zero mean algorithm
    fName = 'None'; %name of the integrand
    figSavePath = ''; %path where to save he figures
    visiblePlot = false; %make plots visible
    debugEnable = false; %enable debug prints
    gaussianCheckEnable = false; %enable plot to check Guassian pdf
    avoidCancelError = true;
    % so the posterior error is a Student-t distribution
    stopCriterion = 'MLE'; % Available options {'MLE', 'GCV', 'full'}
  end
 
  properties (SetAccess = private)
    ff = []; %integrand after the periodization transform
    mvec = []; % n = 2^m
    uncert = 0;  % quantile value for the error bound
    gen_vec = []; % generator for the Lattice points
    
    fullBayes = false; % Full Bayes - assumes m and s^2 as hyperparameters,
    GCV = false; % Generalized cross validation
    vdc_order = false; % use Lattice points generated in vdc order 

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
    function obj = cubBayesLattice_g(varargin)  %Constructor
      
      if nargin > 0
        iStart = 1;
        if isa(varargin{1},'cubBayesLattice_g')
          obj = copy(varargin{1});
          iStart = 2;
        end
        % parse each input argument passed
        if nargin >= iStart
          wh = find(strcmp(varargin(iStart:end),'f'));
          if ~isempty(wh), obj.f = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'dim'));
          if ~isempty(wh), obj.dim = varargin{wh+iStart}; end
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
        end
      end
      
      if strcmp(obj.stopCriterion, 'full')
        % use full Bayes
        obj.fullBayes = true;
      elseif strcmp(obj.stopCriterion, 'GCV')
        % use Generalized cross validation
        obj.GCV = true;
      else
        % use Maximum likelihood estimation
        obj.fullBayes = false;
        obj.GCV = false;
      end
      
      % uncertainity : two sided confidence
      % i.e, alpha percent quantile
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
      obj.errorBdAll = zeros(length_mvec,1);
      obj.muhatAll = zeros(length_mvec,1);
      obj.aMLEAll = zeros(length_mvec,1);
      obj.lossMLEAll = zeros(length_mvec,1);
      obj.timeAll = zeros(length_mvec,1);
      obj.dscAll = zeros(length_mvec,1);
      obj.s_All = zeros(length_mvec,1);
      
    end
    
    % computes the integral
    function [muhat,out] = compInteg(obj)
      
      % comment this line of code to use GPU for computations
      gpuArray = @(x) x;   gather = @(x) x;
      
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
        [stop_flag, muhat] = stopping_criterion(obj, xun_, ftilde_, iter, m);
        
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
      out.ErrBdAll = obj.errorBdAll;
      out.muhatAll = obj.muhatAll;
      out.mvec = obj.mvec;
      out.aMLEAll = obj.aMLEAll;
      out.timeAll = obj.timeAll;
      out.s_All = obj.s_All;
      out.dscAll = obj.dscAll;
      out.absTol = obj.absTol;
      out.relTol = obj.relTol;
      out.shift = shift;
      out.stopAtTol = obj.stopAtTol; 
      
      % convert from gpu memory to local
      muhat=gather(muhat);
      out=gather(out);
      
    end
    
    %% Efficient FFT computation algorithm
    function [ftilde_,xun_,xpts_] = iter_fft(obj,iter,shift,xun,xpts,ftildePrev)
      m = obj.mvec(iter);
      n = 2^m;
      gpuArray = @(x) x;   
      
      % In every iteration, "n" number_of_points is doubled, but FFT is only
      % computed for the newly added points. Previously computed FFT is
      % reused.
      if iter == 1
        % In the first iteration compute full FFT
        xun_ = mod(bsxfun(@times,(0:1/n:1-1/n)',obj.gen_vec),1);
        xpts_ = mod(bsxfun(@plus,xun_,shift),1);  % shifted
        
        % Compute initial FFT
        ftilde_ = fft(gpuArray(obj.ff(xpts_))); %evaluate integrand's fft
      else
        xunnew = mod(bsxfun(@times,(1/n:2/n:1-1/n)',obj.gen_vec),1);
        xnew = mod(bsxfun(@plus,xunnew,shift),1);
        [xun_, xpts_] = obj.merge_pts(xun, xunnew, xpts, xnew, n, obj.dim);
        mnext=m-1;
        
        % Compute FFT on next set of new points
        ftildeNextNew = fft(gpuArray(obj.ff(xnew)));
        if obj.debugEnable
          cubBayesLattice_g.alertMsg(ftildeNextNew, 'Nan', 'Inf');
        end
        
        % combine the previous batch and new batch to get FFT on all points
        ftilde_ = obj.merge_fft(ftildePrev, ftildeNextNew, mnext);
      end
    end
    
    function [ftilde_,xun_,xpts_] = iter_fft_vdc(obj,iter,shift,xun,xpts,ftildePrev)
      m = obj.mvec(iter);
      n = 2^m;
      gpuArray = @(x) x;
      
      function [xun, x] = merge_pts_vdc(xun, xunnew, x, xnew)
        xun = [xun; xunnew];
        x = [x; xnew];
      end
      
      % In every iteration, "n" number_of_points is doubled, but FFT is only
      % computed for the newly added points. Previously computed FFT is
      % reused.
      if iter == 1
        % in the first iteration compute the full FFT
        xun_ = obj.simple_lattice_gen(n,obj.dim,true);
        xpts_ = mod(bsxfun(@plus,xun_,shift),1);  % shifted
        
        % Compute initial FFT
        ftilde_ = obj.fft_DIT(gpuArray(obj.ff(xpts_)),m); %evaluate integrand's fft
      else
        xunnew = obj.simple_lattice_gen(n,obj.dim,false);
        xnew = mod(bsxfun(@plus,xunnew,shift),1);
        mnext=m-1;
        
        % Compute FFT on next set of new points
        ftildeNextNew = obj.fft_DIT(gpuArray(obj.ff(xnew)),mnext);
        if obj.debugEnable
          cubBayesLattice_g.alertMsg(ftildeNextNew, 'Nan', 'Inf');
        end
        
        [xun_, xpts_] = merge_pts_vdc(xun, xunnew, xpts, xnew);
        % combine the previous batch and new batch to get FFT on all points
        ftilde_ = obj.merge_fft(ftildePrev, ftildeNextNew, mnext);
      end
    end
    
    % decides if the user defined error threshold is met
    function [success,muhat] = stopping_criterion(obj, xpts, ftilde, iter, m)
      
      n=2^m;
      success = false;
      %Compute MLE parameter
      lnaMLE = fminbnd(@(lna) ...
        ObjectiveFunction(obj, exp(lna),xpts,ftilde,m), -5,5,optimset('TolX',1e-2));
      aMLE = exp(lnaMLE);
      [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj, aMLE,xpts,ftilde,m);
      
      %Check error criterion
      % compute DSC :
      if obj.fullBayes==true
        % full bayes
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/n);
        else
          DSC = abs((Lambda(1)/n) - 1);
        end
        % 1-alpha two sided confidence interval
        out.ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/(n-1));
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
        out.ErrBd = obj.uncert*sqrt(DSC * (RKHSnorm) /C_Inv_trace);
      else
        % empirical bayes
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/(n + Lambda_ring(1)));
        else
          DSC = abs(1 - (n/Lambda(1)));
        end
        out.ErrBd = obj.uncert*sqrt(DSC * RKHSnorm/n);
      end
      
      if obj.arbMean==true % zero mean case
        muhat = ftilde(1)/n;
      else % non zero mean case
        muhat = ftilde(1)/Lambda(1);
      end
      muminus = muhat - out.ErrBd;
      muplus = muhat + out.ErrBd;
      
      % store intermediate values for post analysis
      % store the debug information
      obj.dscAll(iter) = sqrt(DSC);
      obj.s_All(iter) = sqrt(RKHSnorm/n);
      obj.muhatAll(iter) = muhat;
      obj.errorBdAll(iter) = out.ErrBd;
      obj.aMLEAll(iter) = aMLE;
      obj.lossMLEAll(iter) = loss;
      
      if obj.gaussianCheckEnable == true
        % plots the transformed and scaled integrand values as normal plot
        % Useful to verify the assumption, integrand was an instance of a Gaussian process
        CheckGaussianDensity(obj, ftilde, Lambda)
      end
      
      if 2*out.ErrBd <= ...
          max(obj.absTol,obj.relTol*abs(muminus)) + max(obj.absTol,obj.relTol*abs(muplus))
        if obj.errorBdAll(iter)==0
          obj.errorBdAll(iter) = eps;
        end
        % stopping criterion achieved
        success = true;
      end
      
    end
    
    % objective function to estimate parmaeter theta
    % MLE : Maximum likelihook estimation
    % GCV : Generalized cross validation
    function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,a,xun,ftilde,m)
      
      n = length(ftilde);
      if obj.order==4 || obj.order==2
        [Lambda, Lambda_ring] = obj.kernel(xun,obj.order,a,obj.avoidCancelError,...
          m,obj.vdc_order,obj.debugEnable);
      else
        error('Unsupported Bernoulli polyn order !');
      end
      ftilde = abs(ftilde);  % remove any negative values
      
      % compute RKHSnorm 
      temp = (abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0))) ;
      
      % compute loss 
      if obj.GCV==true
        % GCV 
        temp_gcv = (real(ftilde(Lambda~=0))./(Lambda(Lambda~=0))).^2 ;
        loss1 = log(sum(1./Lambda(Lambda~=0)));
        loss2 = log(sum(temp_gcv(2:end)));
        % ignore all zero val eigenvalues
        loss = log(sum(temp_gcv(2:end))) - 2*log(sum(1./Lambda(Lambda~=0)));
        
        if obj.arbMean==true
          RKHSnorm = sum(temp_gcv(2:end))/n;
        else
          RKHSnorm = sum(temp_gcv)/n;
        end
      else
        % MLE
        if obj.arbMean==true
          RKHSnorm = sum(temp(2:end))/n;
          temp_1 = sum(temp(2:end));
        else
          RKHSnorm = sum(temp)/n;
          temp_1 = sum(temp);
        end
        
        loss1 = sum(log(Lambda(Lambda~=0)));
        loss2 = n*log(temp_1);
        % ignore all zero eigenvalues
        loss = sum(log(Lambda(Lambda~=0))) + n*log(temp_1);        
      end
      
      if obj.debugEnable
        cubMLESobol.alertMsg(RKHSnorm, 'Imag');
        cubMLESobol.alertMsg(loss1, 'Inf');
        cubMLESobol.alertMsg(loss2, 'Inf');
        cubMLESobol.alertMsg(loss, 'Inf', 'Imag', 'Nan');
        cubMLESobol.alertMsg(Lambda, 'Imag');
      end
    end
    
    % Plots the transformed and scaled integrand values as normal plots.
    % This is to verify the assumption, integrand was an instance of
    % gaussian process.
    % Normally distributed :
    %    https://www.itl.nist.gov/div898/handbook/eda/section3/normprp1.htm
    % Short Tails :
    %   https://www.itl.nist.gov/div898/handbook/eda/section3/normprp2.htm
    % Long Tails :
    %    https://www.itl.nist.gov/div898/handbook/eda/section3/normprp3.htm
    function CheckGaussianDensity(obj, ftilde, lambda)
      n = length(ftilde);
      ftilde(1) = 0;  % substract by (m_\MLE \vone) to get zero mean
      w_ftilde = (1/sqrt(n))*real(ftilde)./sqrt(real(lambda));
      if obj.visiblePlot==false
        hFigNormplot = figure('visible','off');
      else
        hFigNormplot = figure();
      end
      set(hFigNormplot,'defaultaxesfontsize',16, ...
        'defaulttextfontsize',16, ... %make font larger
        'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
      %normplot(w_ftilde)
      qqplot(w_ftilde);
      set(hFigNormplot, 'units', 'inches', 'Position', [1 1 8 6])
      
      title(sprintf('Normplot %s n=%d Tx=%s', obj.fName, n, obj.ptransform))
      
      % build filename with path to store the plot
      plotFileName = sprintf('%s%s Normplot d_%d bernoulli_%d Period_%s n_%d.png',...
        obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform, n);
      % plotFileName
      saveas(hFigNormplot, plotFileName)
    end
    
    % plots the objective for the MLE of theta, useful for diagnozing any
    % issues with the code
    function minTheta = plotObjectiveFunc(obj)
      
      numM = length(obj.mvec);
      n = 2.^obj.mvec(end);
      xptsun = cubBayesLattice_g.simple_lattice_gen(n,obj.dim,true);
      fx = obj.ff(xptsun);  % Note: periodization transform already applied
        
      %% plot ObjectiveFunction
      lnTheta = -5:0.2:5;
      % build filename with path to store the plot
      plotFileName = sprintf('%s%s_Cost_d%d_r%d_%s_%s.png',...
        obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform, obj.stopCriterion);
      
      costMLE = zeros(numM,numel(lnTheta));
      tstart = tic;
      
      % loop over all the m values
      for iter = 1:numM
        nii = 2^obj.mvec(iter);
        
        eigvalK = zeros(numel(lnTheta),nii);
        ftilde_iter = fft(bitrevorder(fx(1:nii))); %/nii;
        br_xun = bitrevorder(xptsun(1:nii,:));
        
        tic
        %par
        parfor k=1:numel(lnTheta)
          [costMLE(iter,k),eigvalK(k,:)] = ObjectiveFunction(obj, exp(lnTheta(k)),...
            br_xun,ftilde_iter);
        end
        toc
      end
      
      toc(tstart)
      
      if obj.visiblePlot==false
        hFigCost = figure('visible','off');
      else
        hFigCost = figure();
      end
      
      % semilogx
      semilogx(exp(lnTheta),real(costMLE));
      set(hFigCost, 'units', 'inches', 'Position', [1 1 13.5 11.5])
      xlabel('Shape param, \(\theta\)')
      ylabel('MLE Cost, \( \log \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
      axis tight;
      if obj.arbMean
        mType = '\(m \neq 0\)'; % arb mean
      else
        mType = '\(m = 0\)'; % zero mean
      end
      title(sprintf('%s d=%d r=%d Tx=%s %s', obj.fName, obj.dim, obj.order, obj.ptransform, mType));
      [minVal,Index] = min(real(costMLE),[],2);
      
      % mark the min theta values found using fminbnd
      minTheta = exp(lnTheta(Index));
      hold on;
      semilogx(minTheta, minVal, '.');
      if exist('aMLEAll', 'var')
        semilogx(obj.aMLEAll, obj.lossMLEAll, '+');
      end
      temp = string(obj.mvec);
      temp = strcat('\(2^{',temp,'}\)');
      temp(end+1) = '\(\theta_{min_{true}}\)';
      if exist('aMLEAll', 'var')
        temp(end+1) = '\(\theta_{min_{est}}\)';
      end
      legend(temp,'location','best'); axis tight
      saveas(hFigCost, plotFileName)
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
                fprintf('%s has NaN values', inpvarname);
              end
            case 'Inf'
              if any(isinf(varTocheck))
                fprintf('%s has Inf values', inpvarname);
              end
            case 'Imag'
              if ~all(isreal(varTocheck))
                fprintf('%s has complex values', inpvarname)
              end
            otherwise
              fprintf('%s : unknown type check requested !', inpvarname)
          end
        end
      end
      
    end
    
    % Computes modified kernel Km1 = K - 1
    % Useful to avoid cancellation error in the computation of (1 - n/\lambda_1)
    function [Km1, K] = kernel_t(a, const, Bern)
      theta = a*const;
      d = size(Bern, 2);
      
      Kjm1 = theta*Bern(:,1);  % Kernel at j-dim minus One
      Kj = 1 + Kjm1;  % Kernel at j-dim
      
      for j=2:d
        Kjm1_prev = Kjm1; Kj_prev = Kj;  % save the Kernel at the prev dim
        
        Kjm1 = theta*Bern(:,j).*Kj_prev + Kjm1_prev;
        Kj = 1 + Kjm1;
      end
      
      Km1 = Kjm1; K = Kj;
    end
    
    % bernoulli polynomial based kernel
    % C1 : first row of the kernel
    % Lambda : eigen values of the kernel
    % Lambdahat = fft(C1 - 1)
    function [Lambda, Lambda_ring] = kernel(xun,order,a,avoidCancelError,...
        m,vdc_order,debugEnable)
      
      constMult = -(-1)^(order/2)*((2*pi)^order)/factorial(order);
      % constMult = -(-1)^(order/2);
      if order == 2
        bernPoly = @(x)(-x.*(1-x) + 1/6);
      elseif order == 4
        bernPoly = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
      else
        error('Bernoulli order not implemented !');
      end
      
      if avoidCancelError
        % Computes C1m1 = C1 - 1
        % C1_new = 1 + C1m1 indirectly computed in the process
        [C1m1, C1_alt] = cubBayesLattice_g.kernel_t(a, constMult, bernPoly(xun));
        if vdc_order
          Lambda_ring = real(cubBayesLattice_g.fft_DIT(C1m1,m));
        else
          Lambda_ring = real(fft(C1m1));
        end
        
        Lambda = Lambda_ring;
        Lambda(1) = Lambda_ring(1) + length(Lambda_ring);
        % C1 = prod(1 + (a)*constMult*bernPoly(xun),2);  %
        % Lambda = real(fft(C1));
      else
        % direct appraoch to compute first row of the Kernel Gram matrix
        C1 = prod(1 + (a)*constMult*bernPoly(xun),2);  %
        % matlab's builtin fft is much faster and accurate
        Lambda = real(fft(C1));  % remove any negative values
        Lambda_ring = 0;
      end
      
      if debugEnable == true
        Lambda_direct = abs(fft(C1_alt)); % Note: fft output not normalized
        if sum(abs(Lambda_direct-Lambda)) > 1
          fprintf('Possible error: check Lambda_ring computation')
        end
        
        if sum(C1)==length(C1) || Lambda(1)==length(C1)
          %fprintf('debug');
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
        psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16; psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
        f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
      elseif strcmp(ptransform,'C3sin')
        % Sidi C^3 transform
        psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi); psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
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
    
    % generates rank-1 Lattice points in Vander corput sequence order
    function [xlat, z] = simple_lattice_gen(n,d,firstBatch)
      z = cubBayesLattice_g.get_lattice_gen_vec(d);
      
      nmax = n;
      nmin = 1 + n/2;
      if firstBatch==true
        nmin = 1;
      end
      nelem=nmax-nmin+1;
      
      if firstBatch==true
        brIndices=cubBayesLattice_g.vdc(nelem)';
      else
        brIndices=cubBayesLattice_g.vdc(nelem)'+1/(2*(nmin-1));
      end
      xlat = mod(bsxfun(@times,brIndices',z),1);  % unshifted
    end
    
    % Van der Corput sequence in base 2
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
    
    % fft with deimation in time i.e, input is already in 'bitrevorder'
    function y = fft_DIT(y,nmmin)
      %if nmmin~=log2(length(y))
      %  warning('mismatch')
      %end
      %nmmin = log2(length(y));
      %y = bitrevorder(y);
      for l=0:nmmin-1
        nl=2^l;
        nmminlm1=2^(nmmin-l-1);
        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
        coefv=repmat(coef,nmminlm1,1);
        evenval=y(ptind);
        oddval=y(~ptind);
        y(ptind)=(evenval+coefv.*oddval); %/2;
        y(~ptind)=(evenval-coefv.*oddval); %/2;
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
  end
end

