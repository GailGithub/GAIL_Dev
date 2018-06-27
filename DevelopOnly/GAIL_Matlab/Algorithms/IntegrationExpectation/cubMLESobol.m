%CUBMLESobol Bayesian cubature method to estimate the mean of a random
% variable with with Sobol points
%
%   obj = cubMLESobol('f',f,'dim',dim,'absTol',absTol,'relTol',relTol,...
%     'order',order, 'ptransform',ptransform, ...
%     'arbMean',arbMean);
%   tmu = comIteg(obj) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples used are Sobol points.  The default values are n=2^10 and
%   dim = 1. Input f is a function handle that accepts an n x dim matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%   order : order of the kernel to use
%   ptransform : periodization transform to use
%   arbMean : If True assume 'arbitrary Mean' for the Integrand space
%             else assume 'zero mean'

classdef cubMLESobol < handle
  
  properties
    f = @(x) x.^2; %function to integrate
    dim = 1; %dimension of the integrand
    mmin = 10; %min number of points to start with = 2^mmin
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
    debugEnable = false; %enable debug prints
    gaussianCheckEnable = false; %enable plot to check Guassian pdf
    verify_ftilde = false; % enable to verrify iterative computation of ftilde
  end
  
  properties (SetAccess = private)
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
    function obj = cubMLESobol(varargin)  %Constructor
      if nargin > 0
        iStart = 1;
        if isa(varargin{1},'cubMLESobol')
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
      
      %wKernel = @(x)12*( (1/4) - 2.^(floor(log2(x))-1) );
      sobstr = scramble(sobolset(obj.dim),'MatousekAffineOwen'); %generate a Sobol' sequence
      
      % no periodization transfrom required for this algorithm
      obj.ptransform = 'Baker';  % Hack : remove it soon
      ff = cubMLELattice.doPeriodTx(obj.f, obj.ptransform) ;
      % ff = obj.f;
      
      % pick a random value to apply as shift
      %shift = rand(1,obj.dim);
      
      %% Iteratively find the number of points required for the cubature to meet
      % the error threshold
      for iter = 1:numM
        tstart_iter = tic;
        m = obj.mvec(iter);
        n = 2^m;
        
        if iter == 1
          n0 = 1;
          nnext = n;
          xpts = sobstr(n0:nnext,1:obj.dim);  % grab Sobol' points
          fx = gpuArray(ff(xpts)); %evaluate integrand
          ftilde = cubMLESobol.fwht_hs(fx);
        else
          n0=1+n/2;
          nnext = n;
          xptsnext = sobstr(n0:nnext,1:obj.dim);
          xpts = [xpts;xptsnext];
          
          fx = gpuArray(ff(xptsnext));  % initialize for inplace computation
          ftildeNext = cubMLESobol.fwht_hs(fx);
          % ftilde = [(ftilde+ftildeNext)/2; (ftilde-ftildeNext)/2];
          ftilde = [(ftilde+ftildeNext); (ftilde-ftildeNext)];
          
          if obj.verify_ftilde == true
            fx = gpuArray(ff(xpts)); %evaluate integrand
            ftilde_direct = cubMLESobol.fwht_hs(fx);
            if sum(abs(ftilde-ftilde_direct)) > 1
              fprintf('fwht_hs computation is wrong');
            end
          end
        end
        
        
        %br_xpts = bitrevorder(gpuArray(xpts));
        br_xpts = gpuArray(xpts);
        
        %Compute MLE parameter
        lnaMLE = fminbnd(@(lna) ...
          MLEKernel(obj, exp(lna),br_xpts,ftilde), ...
          -5,5,optimset('TolX',1e-2));
        aMLE = exp(lnaMLE);
        [loss,Lambda,RKHSnorm] = MLEKernel(obj, aMLE,br_xpts,ftilde);
        
        if obj.gaussianCheckEnable==true
          PlotToCheckGaussianDensity(obj, ftilde, Lambda)
        end
        
        %Check error criterion
        % DSC = abs(1 - (n/Lambda(1)));
        % Lambda already divided by 'n'
        % DSC = abs(1 - (1/Lambda(1)));
        DSC = abs(1 - (n/Lambda(1)));
        %DSC = abs(Kthat_new(1)/(n + Kthat_new(1)));
        
        % store the debug information
        obj.dscAll(iter) = sqrt(DSC);
        obj.s_All(iter) = sqrt(RKHSnorm/n);
        
        out.ErrBd = 2.58*sqrt(DSC)*sqrt(RKHSnorm/n);
        if obj.arbMean == true % zero mean case
          % compute sum of 'f' divided by 'n'
          muhat = ftilde(1)/n;
          % muhat = ftilde(1);  % ftilde already divided by 'n'
          
        else % non zero mean case
          muhat = ftilde(1)/Lambda(1);
        end
        muminus = muhat - out.ErrBd;
        muplus = muhat + out.ErrBd;
        obj.muhatAll(iter) = muhat;
        obj.errorBdAll(iter) = out.ErrBd;
        obj.aMLEAll(iter) = aMLE;
        obj.lossMLEAll(iter) = loss;
        
        if 2*out.ErrBd <= ...
            max(obj.absTol,obj.relTol*abs(muminus)) + max(obj.absTol,obj.relTol*abs(muplus))
          if obj.errorBdAll(iter)==0
            obj.errorBdAll(iter) = eps; % zero values cause problem in log plots
          end
          
          % if stopAtTol true, exit the loop
          % else, run for for all 'n' values.
          % Used to compute errorValues for 'n' vs error plotting
          if obj.stopAtTol==true
            break
          end
        end
        
        obj.timeAll(iter) = toc(tstart_iter);
        
      end %iteration loop
      
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
    
    
    % plots the objective function for the MLE of theta
    function minTheta = plotMLE_Loss(obj)
      
      n = 2^obj.mmax;
      sobstr = sobolset(obj.dim); %generate a Sobol' sequence
      xpts = sobstr(1:n,1:obj.dim); %grab Sobol' points
      fx = obj.f(xpts);  % No periodization transform required
      numM = length(obj.mvec);
      
      %% plot MLEKernel cost function
      lnTheta = -5:0.2:5;
      
      %fullPath = strcat(obj.figSavePath,'/',obj.fName,'/',obj.ptransform,'/');
      plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s.png', ...
        obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform);
      plotFileName  % just to display it
      
      costMLE = zeros(numM,numel(lnTheta));
      tstart=tic;
      
      for iter = 1:numM
        nii = 2^obj.mvec(iter);
        nii
        
        eigvalK = zeros(numel(lnTheta),nii);
        %br_xpts = bitrevorder(xpts(1:nii, :));
        %ftilde = cubMLESobol.fwht_hs(obj.f(br_xpts)); %/nii;
        ftilde = cubMLESobol.fwht_hs(fx(1:nii)); %, nii, 'hadamard'
        br_xpts = xpts(1:nii,:);
        
        tic
        %par
        parfor k=1:numel(lnTheta)
          [costMLE(iter,k),eigvalK(k,:)] = MLEKernel(obj, exp(lnTheta(k)),...
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
      saveas(hFigCost, plotFileName)
    end
    
    
    % MLE objective function to find the optimal shape parmaeter
    function [loss,Lambda,RKHSnorm,K] = MLEKernel(obj, a, xpts, ftilde)
      
      n = length(ftilde);
      if obj.order==4
        [K, Lambda] = cubMLESobol.kernel(xpts,obj.order,a);
      elseif obj.order==2
        [K, Lambda] = cubMLESobol.kernel(xpts,obj.order,a);
      else
        error('Unsupported Kernel order !');
      end
      
      ftilde = abs(ftilde);  % remove any negative values
      temp = ((ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0))) ;
      
      if obj.arbMean==true
        RKHSnorm = sum(temp(2:end))/n;
        % RKHSnorm = sum(temp(2:end)); % already divided by 'n'
        temp_1 = sum(temp(2:end));
      else
        RKHSnorm = sum(temp)/n;
        % RKHSnorm = sum(temp); % already divided by 'n'
        temp_1 = sum(temp);
      end
      cubMLESobol.alertMsg(temp_1, 'Imag');
      RKHSnormSq = sqrt(RKHSnorm);
      
      cubMLESobol.alertMsg(RKHSnormSq, 'Nan');
      % loss = mean(2*log(LambdaSq)) + log(RKHSnorm);
      % loss = mean(log(Lambda)) + log(RKHSnorm);
      
      loss1 = sum(log(Lambda(Lambda~=0)));
      cubMLESobol.alertMsg(loss1, 'Inf');
      loss2 = n*log(temp_1);
      cubMLESobol.alertMsg(loss2, 'Inf');
      % ignore all zero val eigenvalues
      loss = sum(log(Lambda(Lambda~=0))) + n*log(temp_1);
      %fprintf('Shap %03.3f loss1 %03.3f\t loss2 %03.3f\t Loss %03.3f\t Nzero eigvals %d\n',...
      %    a, loss1, loss2, loss, sum(Lambda~=0)  )
      
      cubMLESobol.alertMsg(loss, 'Inf', 'Imag', 'Nan');
      cubMLESobol.alertMsg(Lambda, 'Imag');
    end
    
    % Plots the transformed and scaled integrand values as normal plots.
    % This is to verify the assumption, integrand was an instance of
    % gaussian process
    function PlotToCheckGaussianDensity(obj ,ftilde, lambda)
      n = length(ftilde);
      w_ftilde = (1/n)*abs(ftilde)./sqrt(abs(lambda));
      figure();
      normplot(w_ftilde)
      title(sprintf('Hist. %s n=%d Tx=%s', obj.fName, n, obj.ptransform))
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
                fprintf('%s has NaN values', inpvarname);
              end
            case 'Inf'
              if isinf(varTocheck)
                fprintf('%s has Inf values', inpvarname);
              end
            case 'Imag'
              if ~isreal(varTocheck)
                fprintf('%s has complex values', inpvarname)
              end
            otherwise
              fprintf('%s : unknown type check requested !', inpvarname)
          end
        end
      end
      
    end
    
    % Builds the difference matrix to compute kernel matrix, useful for
    % debugging or verification purpose
    % x : input points of size [n,d]
    function dm = diffMatrix(x)
      [n, d] = size(x);
      A = reshape(x, n,1,d);
      dm = repmat(A*n, [1, n, 1]);
      A1 = reshape(x, 1,n,d);
      dm1 = repmat(A1*n, [n, 1, 1]);
      dm = bitxor(dm, dm1)/n;  % bitwise subtraction
    end
    
    % walsh kernel
    function [K, Lambda] = kernel(xpts,order,theta)
      
      if order==2
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
      Lambda = abs(cubMLESobol.fwht_hs(K));
      cubMLESobol.alertMsg(Lambda, 'Imag');
    end
    
    % High order walsh kernel
    function [K, Lambda] = kernel_high(xpts,order,theta)
      
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
        omega2_1D = @(x, shp)(1 + shp*(-2*x + ts2(x)));
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
        omega3_1D = @(x, shp)(1 + shp*(-2*x + s2(x) + ts3(x)));
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
      Lambda = abs(cubMLESobol.fwht_hs(K));
      
      if 0 % Create the full kernel matrix to compare
        dm = cubMLESobol.diffMatrix(xun);
        lambda = eig(omega2(dm,a))/length(xun);
        temp = lambda./sort(Lambda);
        if max(temp) > 2
          fprintf('debug')
          figure; plot(temp)
        end
      end
      
      cubMLESobol.alertMsg(Lambda, 'Imag');
    end
    
    % plots the kernel with different params
    % useful for visualization and debugging
    % For Ex:
    %   obj = cubMLESobol(); obj.demoKernel(1024,2,2,1)
    function demoKernel(npts,ndims,order,theta)
      
      sobstr=sobolset(ndims); %generate a Sobol' sequence
      xpts = sobstr(1:npts,1:ndims); %grab Sobol' points
      K1 = cubMLESobol.kernel(xpts,order,theta);
      if ndims==1
        figure(); plot(xpts, K1, '.'); grid on; axis([0 1 -1 2])
      elseif ndims==2
        ns = sqrt(npts);
        if floor(ns) ~= ns
          error('number of points should be n = 2^m');
        end
        x = sobstr(1:ns,1:1); %grab Sobol' points
        x = sort(x);
        [X,Y] = meshgrid(x,x);
        xpts = [X(:) Y(:)];
        K1 = cubMLESobol.kernel(xpts,order,theta);
        Z = reshape(K1, ns,ns);
        figure(); surf(X,Y, Z); axis tight
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$\omega_2$')
      elseif ndims==2 && false
        figure(); plot3(xpts(:,1), xpts(:,2), K1,'.'); axis tight; grid on
        xlabel('$x_1$')
        ylabel('$x_2$')
        zlabel('$\omega_2$')
        
      elseif ndims==2 && false  % does not work
        ns = sqrt(npts);
        %xpts = sortrows(xpts);
        X = reshape(xpts(:,1), ns,ns); Y= reshape(xpts(:,2), ns,ns);
        Z = reshape(K1, ns,ns);
        figure(); surf(X,Y, Z ,'FaceAlpha',0.5); axis tight
      else
        error('demoKernel: ndims > 2 not implemented');
      end
    end
  end % end of static functions
end