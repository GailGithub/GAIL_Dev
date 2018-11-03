%CUBBAYESNET_G Bayesian cubature method to estimate the mean of a random
% variable with with Sobol points
%
%   obj = cubBayesNet_g('f',f,'dim',dim,'absTol',absTol,'relTol',relTol,...
%     'order',order, 'arbMean',arbMean);
%   tmu = comIteg(obj) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples used are Sobol points.  The default values are n=2^10 and
%   dim = 1. Input f is a function handle that accepts an n x dim matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%   order : order of the kernel to use
%   ptransform : periodization transform to use
%   arbMean : If True assume 'arbitrary Mean' for the Integrand space
%             else assume 'zero mean'

classdef cubBayesNet_g < handle
  
  properties
    f = @(x) x.^2; %function to integrate
    dim = 1; %dimension of the integrand
    mmin = 4; %min number of points to start with = 2^mmin
    mmax = 20; %max number of points allowed = 2^mmax
    absTol = 0.01; %absolute tolerance
    relTol = 0; %relative tolerance
    order = 2; %order of the kernel
    
    ptransform = 'none'; %periodization transform
    stopAtTol = true; %automatice mode: stop after meeting the error tolerance
    arbMean = false; %by default use zero mean algorithm
    fName = 'None'; %name of the integrand
    figSavePath = ''; %path where to save he figures
    visiblePlot = true; %make plots visible
    debugEnable = true; %enable debug prints
    gaussianCheckEnable = false; %enable plot to check Guassian pdf
    verify_ftilde = false; % enable to verrify iterative computation of ftilde
    net_type = 'Sobol'; % type of nets. Options: NX or Sobol
  end
  
  properties (SetAccess = private)
    avoidCancelError = true;

    gen_h = 0;
    gen_h_un = 0;
    
    %CovProp; %square root and determinant of covariance matrix
    mvec = [];
    
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
    function obj = cubBayesNet_g(varargin)  %Constructor
      if nargin > 0
        iStart = 1;
        if isa(varargin{1},'cubBayesNet_g')
          obj = copy(varargin{1});
          iStart = 2;
        end
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
          % wh = find(strcmp(varargin(iStart:end),'ptransform'));
          % if ~isempty(wh), obj.ptransform = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'arbMean'));
          if ~isempty(wh), obj.arbMean = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'stopAtTol'));
          if ~isempty(wh), obj.stopAtTol = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'figSavePath'));
          if ~isempty(wh), obj.figSavePath = varargin{wh+iStart}; end
          wh = find(strcmp(varargin(iStart:end),'fName'));
          if ~isempty(wh), obj.fName = varargin{wh+iStart}; end
        end
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
      
      if obj.verify_ftilde==true
        warning('Caution: debugEnable is set, this will increase computation time !')
      end
    end
    
    
    % computes the integral
    function [muhat,out] = compInteg(obj)
      % comment this line of code to use GPU for computations
      gpuArray = @(x) x;   gather = @(x) x;
      obj = gpuArray(obj);
      
      tstart = tic; %start the clock
      
      numM = length(obj.mvec);
      
%       %wKernel = @(x)12*( (1/4) - 2.^(floor(log2(x))-1) );
%       sobstr = scramble(sobolset(obj.dim),'MatousekAffineOwen'); %generate a Sobol' sequence
%       % gen_name = 'nx_b2_m30_s4_Cs.col';
%       % gen_name = 'nx_b2_m30_s12_Cs.col';
%       % gen_name = 'nx_s5_alpha2_m32.col';
%       % gen_name = 'nx_s5_alpha3_m32.col';
%       gen_name = 'nxmats\nx_s5_alpha6_m32.col';
%       gen_name = 'sobolmats\sobol_alpha3_Bs53.col';
%       genmat = load(['dirk_nuyens\qmc-generators\DIGSEQ\' gen_name]);
%       digitalseq_b2g('init0', genmat)
      
      % Initialize the points generator
      gen_digital_nets(obj,true);
      
      % no periodization transfrom required for this algorithm
      % obj.ptransform = 'C1sin';  % Hack : remove it soon
      % ff = cubBayesLattice_g.doPeriodTx(obj.f, obj.ptransform) ;
      ff = obj.f;
      
      %% Iteratively find the number of points required for the cubature to meet
      % the error threshold
      for iter = 1:numM
        tstart_iter = tic;
        m = obj.mvec(iter);
        n = 2^m;
        
        if iter == 1
          n0 = 1;
          nnext = n;
          % xpts = sobstr(n0:nnext,1:obj.dim);  % grab Sobol' points
          % xpts = digitalseq_b2g(obj.dim, nnext-n0+1)';
          [xpts,xpts_un] = gen_digital_nets(obj,false,n);
          
          fx = gpuArray(ff(xpts)); %evaluate integrand
          ftilde = cubBayesNet_g.fwht_hs(fx);
        else
          n0=1+n/2;
          nnext = n;
          % xptsnext = sobstr(n0:nnext,1:obj.dim);
          % xptsnext = digitalseq_b2g(obj.dim, nnext-n0+1)';
          [xptsnext,xptsnext_un] = gen_digital_nets(obj,false,n/2);

          xpts = [xpts;xptsnext];
          
          fx = gpuArray(ff(xptsnext));  % initialize for inplace computation
          ftildeNext = cubBayesNet_g.fwht_hs(fx);
          % ftilde = [(ftilde+ftildeNext)/2; (ftilde-ftildeNext)/2];
          ftilde = [(ftilde+ftildeNext); (ftilde-ftildeNext)];
          
          if obj.verify_ftilde == true
            fx = gpuArray(ff(xpts)); %evaluate integrand
            ftilde_direct = cubBayesNet_g.fwht_hs(fx);
            if std(abs(ftilde-ftilde_direct)) > 0.1
              fprintf('fwht_hs computation is wrong');
            end
          end
        end
        
        
        %br_xpts = bitrevorder(gpuArray(xpts));
        br_xpts = gpuArray(xpts);
        
        %Compute MLE parameter
        if 0
          lnaMLE = fminbnd(@(lna) ...
            ObjectiveFunction(obj, exp(lna),br_xpts,ftilde), ...
            -5,5,optimset('TolX',1e-2));
          thetaOpt = exp(lnaMLE);
        else
          lnaMLE = fminsearch(@(lna) ...
            ObjectiveFunction(obj, exp(lna),br_xpts,ftilde), ...
            0,optimset('TolX',1e-2));
          thetaOpt = exp(lnaMLE);
        end
        [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj, thetaOpt,br_xpts,ftilde);
        
        if obj.gaussianCheckEnable==true
          PlotToCheckGaussianDensity(obj, ftilde, Lambda)
        end
        
        %Check error criterion
        % empirical bayes
        if obj.avoidCancelError
          DSC = abs(Lambda_ring(1)/(n + Lambda_ring(1)));
        else
          DSC = abs(1 - (n/Lambda(1)));
        end
        
        % store the debug information
        obj.dscAll(iter) = sqrt(DSC);
        obj.s_All(iter) = sqrt(RKHSnorm/n);
        
        out.ErrBd = 2.58*sqrt(DSC*RKHSnorm/n);
        if obj.arbMean == true % zero mean case
          % muhat is just the sample mean
          muhat = ftilde(1)/n;
          
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
          % if obj.errorBdAll(iter)==0
          %   obj.errorBdAll(iter) = eps; % zero values cause problem in log plots
          % end
          
          % if stopAtTol true, exit the loop
          % else, run for for all 'n' values.
          % Used to compute errorValues for 'n' vs error plotting
          if obj.stopAtTol==true
            break
          end
        end
        
        obj.timeAll(iter) = toc(tstart_iter);
        
        fprintf('thetaOpt=%1.3f, n=%d, ErrBd=%1.2e, time=%1.2f, absTol=%1.1e\n', ...
          thetaOpt,n,out.ErrBd,obj.timeAll(iter),obj.absTol);
        
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
      
      % convert from gpu memory to local
      muhat = gather(muhat);
      out = gather(out);
      muhat;   % let it to print
      out;
    end  %function
    
    
    % generates higher order interlaced digital nets or sobol points
    function [xpts,xpts_un] = gen_digital_nets(obj,init,n)
      if init==true
        if obj.order==1
          obj.gen_h = scramble(sobolset(obj.dim),'MatousekAffineOwen'); %generate a Sobol' sequence
          obj.gen_h_un = sobolset(obj.dim); %generate a Sobol' sequence
        else
          if strcmp(obj.net_type,'NX')
            switch obj.order
              case 2
                gen_name = 'nxmats\nx_s5_alpha2_m32.col';
              case 3
                gen_name = 'nxmats\nx_s5_alpha3_m32.col';
              otherwise
                error('Higher order=%d NX net not supported', obj.order)
            end
          else  % Sobol
            switch obj.order
              case 2
                gen_name = 'sobolmats\sobol_alpha2_Bs53.col';
              case 3
                gen_name = 'sobolmats\sobol_alpha3_Bs53.col';
              otherwise
                error('Higher order=%d sobol net not supported', obj.order)
            end
          end
          genmat = load(['dirk_nuyens\qmc-generators\DIGSEQ\' gen_name]);
          digitalseq_b2g('init0', genmat)
        end
      else
        if obj.order==1
          xpts = obj.gen_h(1:n,1:obj.dim); %grab Sobol' points
          xpts_un = obj.gen_h_un(1:n,1:obj.dim); %grab unscrambled Sobol' points
        else
          xpts_un = digitalseq_b2g(obj.dim, n)';
          xpts = xpts_un;  % Scrambling not yet implemented
        end
      end
    end
    
    
    % plots the objective function for the optimal shape parameter
    function [minTheta, hFigCost] = plotObjectiveFunc(obj)
      
      n = 2^obj.mmax;
      %sobstr = sobolset(obj.dim); %generate a Sobol' sequence
      %xpts = sobstr(1:n,1:obj.dim); %grab Sobol' points
      gen_digital_nets(obj,true);
      [xpts] = gen_digital_nets(obj,false,n);
      fx = obj.f(xpts);  % No periodization transform required
      numM = length(obj.mvec);
      
      %% plot range for the objective function
      lnTheta = -5:0.2:5;
      
      %fullPath = strcat(obj.figSavePath,'/',obj.fName,'/',obj.ptransform,'/');
%       plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s.png', ...
%         obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform);
%       plotFileName  % just to display it
      
      costMLE = zeros(numM,numel(lnTheta));
      tstart=tic;
      
      for iter = 1:numM
        nii = 2^obj.mvec(iter);
        nii
        
        eigvalK = zeros(numel(lnTheta),nii);
        %br_xpts = bitrevorder(xpts(1:nii, :));
        %ftilde = cubBayesNet_g.fwht_hs(obj.f(br_xpts)); %/nii;
        ftilde = cubBayesNet_g.fwht_hs(fx(1:nii)); %, nii, 'hadamard'
        br_xpts = xpts(1:nii,:);
        
        tic
        %par
        parfor k=1:numel(lnTheta)
          [costMLE(iter,k),eigvalK(k,:)] = ObjectiveFunction(obj, exp(lnTheta(k)),...
            br_xpts,ftilde);
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
      % saveas(hFigCost, plotFileName)
    end
    
    
    % Objective function to find the optimal shape parmaeter
    function [loss,Lambda,Lambda_ring,RKHSnorm] = ObjectiveFunction(obj,theta,xpts,ftilde)
      
      n = length(ftilde);
      % Compute the eigenvalues of the covariance matrix
      [Lambda,Lambda_ring] = cubBayesNet_g.kernel(xpts,obj.order,theta,obj.avoidCancelError);
      
      % Not required ftilde = abs(ftilde);  % remove any negative/imaginary values
      temp = ((ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)));
      
      if obj.arbMean==true
        RKHSnorm = sum(temp(2:end))/n;
        % RKHSnorm = sum(temp(2:end)); % already divided by 'n'
        temp_1 = sum(temp(2:end));
      else
        RKHSnorm = sum(temp)/n;
        % RKHSnorm = sum(temp); % already divided by 'n'
        temp_1 = sum(temp);
      end
      cubBayesNet_g.alertMsg(temp_1, 'Imag');

      % loss = mean(log(Lambda)) + log(RKHSnorm);
      
      loss1 = sum(log(Lambda(Lambda~=0)))/n;
      loss2 = log(temp_1);
      % ignore all zero val eigenvalues
      loss = loss1 + loss2;
      %fprintf('Shap %03.3f loss1 %03.3f\t loss2 %03.3f\t Loss %03.3f\t Nzero eigvals %d\n',...
      %    a, loss1, loss2, loss, sum(Lambda~=0)  )
      
      cubBayesNet_g.alertMsg(loss1, 'Inf');
      cubBayesNet_g.alertMsg(loss2, 'Inf');
      cubBayesNet_g.alertMsg(loss, 'Inf', 'Imag', 'Nan');
      cubBayesNet_g.alertMsg(Lambda, 'Imag');
    end
    
    % Plots the transformed and scaled integrand values as normal plots.
    % This is to verify the assumption, integrand was an instance of
    % gaussian process
    function PlotToCheckGaussianDensity(obj ,ftilde, lambda)
      n = length(ftilde);
      w_ftilde = (1/n)*abs(ftilde)./sqrt(abs(lambda));
      figure();
      normplot(w_ftilde)
      title(sprintf('Hist. %s n=%d order=%d', obj.fName, n, obj.order))
    end
  end % end of methods
  
  % static methods are the ones that do not need to access obj
  methods(Static)
    % compute fast walsh transform in 'hadamard' ordering
    function t = fwht_hs(fx)
      [n, ~] = size(fx);
      t = fwht(fx,n,'hadamard')*n;
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
      theta = aconst;
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
    
    % Interface for all the kernels and to call cancellation error fix
    function [Lambda,Lambda_ring] = kernel(xpts,order,theta,avoidCancelError)
      [n,dim] = size(xpts);
      kernelFunc = cubBayesNet_g.BuildKernelFunc(order);
      if avoidCancelError
        [C1m1] = cubBayesLattice_g.kernel_t(theta, kernelFunc(xpts));
        % eigenvalues must be real : Symmetric pos definite Kernel
        Lambda_ring = real(cubBayesNet_g.fwht_hs(C1m1));
        
        Lambda = Lambda_ring;
        if any(Lambda_ring < 0)
          warning('eigen values are negative!\n')
        end
        Lambda(1) = Lambda_ring(1) + n;
      else
        C1 = prod(1 + theta*kernelFunc(xpts),2);
        Lambda = real(cubBayesNet_g.fwht_hs(C1));
        Lambda_ring = Lambda;
        Lambda_ring(1) = Lambda_ring(1) - n;
      end
      
      if 0 % Create the full kernel matrix to compare
        dm = cubBayesNet_g.diffMatrix(xpts);
        K = 1 + theta*kernelFunc(dm);
        if dim > 1
          K = prod(K, 3);  % for d>1 reduce the 3rd dim by multiplication
        else
        end
        lambda = eig(K);
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
        kernFunc = @(x)(12*( (1/6) - 2.^(floor(log2(x))-1) ));
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
    
    %{
    % walsh kernel
    function [K, Lambda] = kernel_low(xpts,order,theta)
      
      if order==1
        wKernel_1D = @(x, shp)(1 + shp*12*( (1/6) - 2.^(floor(log2(x))-1) ));
        [~, dim] = size(xpts);
        if dim > 1
          wKernel2 = @(x, shp)prod(wKernel_1D(x,shp), ndims(x));
        else
          wKernel2 = wKernel_1D;
        end
        K = wKernel2(xpts, theta);
      else
        error('kernel order not yet supported');
      end
      Lambda = abs(cubBayesNet_g.fwht_hs(K));
      cubMLESobol.alertMsg(Lambda, 'Imag');
    end
    
    % High order walsh kernel
    function [K, Lambda] = kernel_high_(xpts,order,theta)
      
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
      
      if order==2
        %omega2 = @(x, a)prod(1 + theta*(s1(x) + ts2(x) - 1), ndims(x));
        % to avoid subtracting "1", s1 is used directly
        %omega2 = @(x, a)prod(1 + theta*(-2*x + ts2(x)), ndims(x));
        omega2_1D = @(x, shp)(1 + shp*(1-2*x + ts2(x)));
        [n, dim] = size(xpts);
        if dim > 1
          omega2 = @(x, shp)prod(omega2_1D(x,shp), ndims(x));
        else
          omega2 = omega2_1D;
        end
        % figure(2); x=[0:0.001:1]; plot(x, omega2(x), '.'); grid on; axis([0 1 -1 2])
        
        if theta==1 % debug
          %fprintf('Shape a==1\n')
        end
        K = omega2(xpts, theta);
        
      elseif order==3
        omega3_1D = @(x, shp)(1 + shp*(1-2*x + s2(x) + ts3(x)));
        [~, dim] = size(xpts);
        if dim > 1
          omega3 = @(x, shp)prod(omega3_1D(x,shp), ndims(x));
        else
          omega3 = omega3_1D;
        end
        K = omega3(xpts, theta);
      else
        error('kernel order not yet supported');
      end
      Lambda = abs(cubBayesNet_g.fwht_hs(K));
      
      if 0 % Create the full kernel matrix to compare
        dm = cubBayesNet_g.diffMatrix(xpts);
        lambda = eig(omega2(dm,a))/length(xpts);
        temp = lambda./sort(Lambda);
        if max(temp) > 2
          fprintf('debug')
          figure; plot(temp)
        end
      end
      
      cubMLESobol.alertMsg(Lambda, 'Imag');
    end
    %}
       
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
        figure(); plot(xpts, C1, '.', 'MarkerSize', 10); grid on; axis([0 1 -1 3])
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
  end % end of static functions
end