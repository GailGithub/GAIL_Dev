%CUBMLESobol Monte Carlo method to estimate the mean of a random variable
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
        f = @(x) x.^2; %function
        dim = 1; %dimension
        mmin = 5; %min number of points to start with = 2^mmin
        mmax = 20; %max number of points allowed = 2^mmax
        absTol = 0.01; %absolute tolerance
        relTol = 0; %relative tolerance
        order = 2; %order of the kernel
        
        ptransform = 'none'; %periodization transform
        stopAtTol = false; %stop after meeting the error tolerance
        arbMean = false; %by default use zero mean algorithm
        fName = 'None'; %name of the integrand
        figSavePath = ''; %path where to save he figures
        visiblePlot = true; %make plots visible
        debugEnable = false; %enable debug prints
    end
    
    properties (SetAccess = private)
        %CovProp; %square root and determinant of covariance matrix
        
        % variables to save debug info in each iteration
        errorBdAll = zeros(length(mvec),1);
        muhatAll = zeros(length(mvec),1);
        aMLEAll = zeros(length(mvec),1);
        lossMLEAll = zeros(length(mvec),1);
        timeAll = zeros(length(mvec),1);
        dscAll = zeros(length(mvec),1);
        s_All = zeros(length(mvec),1);
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
        end
        
        % computes the integral
        function [muhat,out] = compInteg(obj)
            
            % comment this line of code to use GPU for computations
            gpuArray = @(x) x;   gather = @(x) x;
            
            tstart = tic; %start the clock
            % define min and max number of points allowed in the automatic cubature
            
            mvec = obj.mmin:obj.mmax;
            numM = length(mvec);
            
            %wKernel = @(x)12*( (1/4) - 2.^(floor(log2(x))-1) );
            sobstr = sobolset(obj.dim); %generate a Sobol' sequence
            
            % no periodization transfrom required for this algorithm
            ff = obj.f;
            
            for iter = 1:numM
                tstart_iter = tic;
                m = mvec(iter);
                n = 2^m;
                
                if iter==1
                    n0=1;
                    nnext=n;
                    xpts=sobstr(n0:nnext,1:obj.dim); %grab Sobol' points
                    fx=gpuArray(ff(xpts)); %evaluate integrand
                    ftilde = fwht(fx);
                else
                    n0=1+n/2;
                    nnext=n;
                    xptsnext=sobstr(n0:nnext,1:out_param.d); %grab Sobol' points
                    xpts=[xpts;xptsnext];
                    fx=gpuArray(ff(xptsnext));  % initialize for inplace computation
                    ftildeNext=fwht(fx);
                    ftilde=[ftilde;ftildeNext];
                end
                
                %br_xpts = bitrevorder(gpuArray(xpts));
                br_xpts = gpuArray(xpts);
                
                %Compute MLE parameter
                lnaMLE = fminbnd(@(lna) ...
                    MLEKernel(exp(lna),br_xpts,ftilde), ...
                    -5,5,optimset('TolX',1e-2));
                aMLE = exp(lnaMLE);
                [loss,Ktilde,RKHSnorm,K] = MLEKernel(aMLE,br_pts,ftilde);
                
                %Check error criterion
                DSC = abs(1 - (n/Ktilde(1)));
                %DSC = abs(Kthat_new(1)/(n + Kthat_new(1)));
                
                % store the debug information
                obj.dscAll(iter) = sqrt(DSC);
                obj.s_All(iter) = sqrt(RKHSnorm/n);
                
                out.ErrBd = 2.58*sqrt(DSC)*sqrt(RKHSnorm/n);
                if obj.arbMean==true % zero mean case
                    muhat = ftilde(1)/n;
                else % non zero mean case
                    muhat = ftilde(1)/Ktilde(1);
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
            out.time = toc(tstart);   % let it to print
            out.ErrBdAll = obj.errorBdAll;
            out.muhatAll = obj.muhatAll;
            out.mvec = mvec;
            out.aMLEAll = obj.aMLEAll;
            out.timeAll = obj.timeAll;
            out.s_All = obj.s_All;
            out.dscAll = obj.dscAll;

            % convert from gpu memory to local
            muhat=gather(muhat);
            out=gather(out);
            muhat
            out
            
        end  %function

        
        % plots the objective function for the MLE of theta
        function minTheta = plotMLE_Loss(obj)
            
            useArbMean = obj.arbMean;
            n = 2^obj.mmax;
            sobstr=sobolset(obj.dim); %generate a Sobol' sequence
            xpts = sobstr(1:n,1:obj.dim); %grab Sobol' points
            fx = obj.f(xpts);  % No periodization transform required
            mvec = obj.mmin:obj.mmax;
            numM = length(mvec);
            
            %% plot MLEKernel cost function
            lnTheta = -5:0.2:5;
            
            %fullPath = strcat(obj.figSavePath,'/',obj.fName,'/',obj.ptransform,'/');
            plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s.png', ...
                obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform);
            plotFileName
            
            costMLE = zeros(numM,numel(lnTheta));
            tstart=tic;
            
            for iter = 1:numM
                nii = 2^mvec(iter);
                nii
                
                eigvalK = zeros(numel(lnTheta),nii);
                %br_xpts = bitrevorder(xpts(1:nii, :));
                %ftilde = fwht(obj.f(br_xpts)); %/nii;
                ftilde = fwht(fx(1:nii)); %, nii, 'hadamard'
                br_xpts = xpts(1:nii,:);
                
                tic
                %par
                for k=1:numel(lnTheta)
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
            
            % semilogx
            loglog(exp(lnTheta),real(costMLE));
            set(hFigCost, 'units', 'inches', 'Position', [0 0 13.5 11.5])
            xlabel('Shape param, \(\theta\)')
            ylabel('MLE Cost, \( \log \log \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
            % ylabel('Log MLE Obj. fun.')
            axis tight;
            if useArbMean
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
            loglog(minTheta,minVal, '.');
            if exist('aMLEAll', 'var')
                loglog(obj.aMLEAll,obj.lossMLEAll, '+');
            end
            temp = string(mvec); temp=strcat('2^',temp);
            temp(end+1) = '\(\theta_{min_{true}}\)';
            if exist('aMLEAll', 'var')
                temp(end+1) = '\(\theta_{min_{est}}\)';
            end
            
            legend(temp,'location','best'); axis tight
            saveas(hFigCost, plotFileName)
            
        end
        
        
        % MLE objective function to find the optimal shape parmaeter
        function [loss,Ktilde,RKHSnorm,K] = MLEKernel(obj, a, xpts, ftilde)
            
            useArbMean = obj.arbMean;
            n = length(ftilde);
            if obj.order==4
                [K, Ktilde] = cubMLESobol.kernel(xpts,obj.order,a);
            elseif obj.order==2
                [K, Ktilde] = cubMLESobol.kernel(xpts,obj.order,a);
            else
                error('Unsupported Bernoulli polyn order !');
            end
            
            ftilde = abs(ftilde);  % remove any negative values
            temp = (abs(ftilde(Ktilde~=0).^2)./(Ktilde(Ktilde~=0))) ;
            
            if useArbMean==true
                RKHSnorm = sum(temp(2:end))/n;
                temp_1 = sum(temp(2:end));
            else
                RKHSnorm = sum(temp)/n;
                temp_1 = sum(temp);
            end
            RKHSnormSq = sqrt(RKHSnorm);
            
            if isnan(RKHSnormSq)
                fprintf('RKHSnormSq NaN \n');
            end
            % loss = mean(2*log(KtildeSq)) + log(RKHSnorm);
            % loss = mean(log(Ktilde)) + log(RKHSnorm);
            
            loss1 = sum(log(Ktilde(Ktilde~=0)));
            if isinf(loss1)
                fprintf('loss1 is infinity \n');
            end
            loss2 = n*log(temp_1);
            if isinf(loss2)
                fprintf('loss2 is infinity \n');
            end
            % ignore all zero val eigenvalues
            loss = sum(log(Ktilde(Ktilde~=0))) + n*log(temp_1);
            fprintf('Shap %03.3f loss1 %03.3f\t loss2 %03.3f\t Loss %03.3f\t Nzero eigvals %d\n',...
                a, loss1, loss2, loss, sum(Ktilde~=0)  )
            
            if isinf(loss)
                fprintf('loss is infinity \n');
            end
            
            if ~isreal(loss)
                fprintf('Total loss has complex vals \n');
            end
            
            if ~isreal(Ktilde)
                fprintf('Ktilde has complex vals \n');
            end
            
            if ~isreal(temp_1)
                fprintf('temp_1 has complex vals \n');
            end
            
            if isnan(loss)
                fprintf('loss NaN \n');
            end
        end
   
    end
    
    methods(Static)
        % Builds the difference matrix to compute kernel matrix
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
        function [K, Ktilde] = kernel(xpts,order,a)
            
            order; % not used for now
            
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
            s1 = @(x)(1-2*x);
            ts2 = @(x)((1-5*t1(x))/2 - (a1(x)-2).*x);
            
            %omega2 = @(x, a)prod(1 + a*(s1(x) + ts2(x) - 1), ndims(x));
            % to avoid subtracting "1", s1 is used directly
            %omega2 = @(x, a)prod(1 + a*(-2*x + ts2(x)), ndims(x));
            omega2_1D = @(x, a)(1 + a*(-2*x + ts2(x)));
            [n, dim] = size(xpts);
            if dim > 1
                omega2 = @(x, a)prod(omega2_1D(x,a), ndims(x));
            else
                omega2 = omega2_1D;
            end
            % figure(2); x=[0:0.001:1]; plot(x, omega2(x), '.'); grid on; axis([0 1 -1 2])
            
            if a==1
                fprintf('Shape a==1\n')
            end
            
            K = omega2(xpts, a);
            Ktilde = abs(fwht(K)); %, n, 'hadamard'

            if 0 % Create the full kernel matrix to compare
                dm = cubMLESobol.diffMatrix(xun);
                lambda = eig(omega2(dm,a))/length(xun);
                temp = lambda./sort(Ktilde);
                if max(temp) > 2
                    fprintf('debug')
                    figure; plot(temp)
                end
            end
            
            if ~isreal(Ktilde)
                fprintf('Ktilde has complex vals \n');
            end
            
        end
        
    end
end