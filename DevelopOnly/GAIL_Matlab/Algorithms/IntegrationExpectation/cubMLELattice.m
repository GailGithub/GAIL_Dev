%CUBMLE Monte Carlo method to estimate the mean of a random variable
%
%   tmu = cubMLELattice(f,absTol,relTol,alpha,nSig,inflate) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples may be of one of several kinds.  The default values are n=2^10 and
%   d = 1 Input f is a function handle that accepts an n x d matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%
%
% This is a heuristic algorithm based on a Central Limit Theorem
% approximation
classdef cubMLELattice < handle
    
    properties
        f = @(x) x.^2; %function to integrate
        dim = 1; %dimension of the integrand
        mmin = 8; %min number of points to start with = 2^mmin
        mmax = 20; %max number of points allowed = 2^mmax
        absTol = 0.01; %absolute tolerance
        relTol = 0; %relative tolerance
        order = 2; %order of the kernel
        ptransform = 'Baker'; %periodization transform
        stopAtTol = false; %stop after meeting the error tolerance
        arbMean = false; %by default use zero mean algorithm
        fName = 'None'; %name of the integrand
        figSavePath = ''; %path where to save he figures
        visiblePlot = true; %make plots visible
        debugEnable = false; %enable debug prints
    end
    
    properties (SetAccess = private)
        ff = []; %integrand after the periodization transform
        mvec = [];
        
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
        function obj = cubMLELattice(varargin)  %Constructor
    
            if nargin > 0
                iStart = 1;
                if isa(varargin{1},'cubMLELattice')
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
            
            % apply periodization transformation to the function
            obj.ff = cubMLELattice.doPeriodTx(obj.f, obj.ptransform);
            
            obj.mvec = obj.mmin:obj.mmax;
            length_mvec = length(obj.mvec);
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
            
            %% Iteratively find the number of points required for the cubature to meet
            % the error threshold
            for iter = 1:numM
                tstart_iter = tic;
                m = obj.mvec(iter);
                n = 2^m;
                
                %Update function values
                %% Efficient FFT computation algorithm
                if iter == 1
                    % in the first iteration compute the full FFT
                    xun = cubMLELattice.simple_lattice_gen(n,obj.dim,true);
                    x = mod(bsxfun(@plus,xun,shift),1);  % shifted
                    
                    ftildeNew=gpuArray(obj.ff(x)); %evaluate integrand
                    
                    % Compute initial FFT
                    for l=0:obj.mmin-1
                        nl=2^l;
                        nmminlm1=2^(obj.mmin-l-1);
                        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
                        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
                        coefv=repmat(coef,nmminlm1,1);
                        evenval=ftildeNew(ptind);
                        oddval=ftildeNew(~ptind);
                        ftildeNew(ptind)=(evenval+coefv.*oddval);
                        ftildeNew(~ptind)=(evenval-coefv.*oddval);
                    end
                else
                    xunnew = cubMLELattice.simple_lattice_gen(n,obj.dim,false);
                    xnew = mod(bsxfun(@plus,xunnew,shift),1);
                    xun = [xun;xunnew];
                    x = [x;xnew];
                    
                    mnext=m-1;
                    ftildeNextNew=gpuArray(obj.ff(xnew));  % initialize for inplace computation
                    
                    if obj.debugEnable==true
                        if any(isnan(ftildeNextNew)) || any(isinf(ftildeNextNew))
                            fprintf('ftildeNextNew NaN \n');
                        end
                    end
                    
                    % Compute initial FFT on next set of new points
                    for l=0:mnext-1
                        nl=2^l;
                        nmminlm1=2^(mnext-l-1);
                        ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
                        coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
                        coefv=repmat(coef,nmminlm1,1);
                        evenval=ftildeNextNew(ptind);
                        oddval=ftildeNextNew(~ptind);
                        ftildeNextNew(ptind)=(evenval+coefv.*oddval);
                        ftildeNextNew(~ptind)=(evenval-coefv.*oddval);
                        
                        if obj.debugEnable==true
                            if any(isnan(ftildeNextNew)) || any(isinf(ftildeNextNew))
                                fprintf('ftildeNextNew NaN \n');
                            end
                        end
                    end
                    
                    % combine the previous batch and new batch to get FFT on all points
                    ftildeNew=[ftildeNew;ftildeNextNew];
                    nl=2^mnext;
                    ptind=[true(nl,1); false(nl,1)];
                    coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
                    coefv=repmat(coef,nmminlm1,1);
                    evenval=ftildeNew(ptind);
                    oddval=ftildeNew(~ptind);
                    ftildeNew(ptind)=(evenval+coefv.*oddval);
                    ftildeNew(~ptind)=(evenval-coefv.*oddval);
                end
                
                ftilde = ftildeNew;
                % ftilde(1) = sum(ff(x)); % correction to avoid round off error
                br_xun = bitrevorder(gpuArray(xun));
                
                %Compute MLE parameter
                lnaMLE = fminbnd(@(lna) ...
                    MLEKernel(obj, exp(lna),br_xun,ftilde), ...
                    -5,5,optimset('TolX',1e-2));
                aMLE = exp(lnaMLE);
                [loss,Ktilde,Kthat_new,RKHSnorm,K] = MLEKernel(obj, aMLE,br_xun,ftilde);
                
                %Check error criterion
                % DSC = abs(1 - (n/Ktilde(1)));
                DSC = abs(Kthat_new(1)/(n + Kthat_new(1)));
                
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
                        obj.errorBdAll(iter) = eps;
                    end
                    
                    % if stopAtTol true, exit the loop
                    % else, run for for all 'n' values.
                    % Used to compute error values for 'n' vs error plotting
                    if obj.stopAtTol==true
                        break
                    end
                end
                
                obj.timeAll(iter) = toc(tstart_iter);
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
            
            % convert from gpu memory to local
            muhat=gather(muhat);
            out=gather(out);
            muhat   % let it to print
            out
            
        end
        
        
        % plots the objective for the MLE of theta
        function minTheta = plotMLE_Loss(obj)
            
            numM = length(obj.mvec);
            n = 2.^obj.mvec(end);
            xun = cubMLELattice.simple_lattice_gen(n,d,true);
            fx = obj.ff(xun);  % Note: periodization transform already applied
            
            %% plot MLEKernel cost function
            lnTheta = -5:0.2:5;
            % fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
            plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s.png',...
                obj.figSavePath, obj.fName, obj.dim, obj.order, obj.ptransform);
            plotFileName
            
            costMLE = zeros(numM,numel(lnTheta));
            tstart = tic;
            
            for iter = 1:numM
                nii = 2^obj.mvec(iter);
                nii
                
                eigvalK = zeros(numel(lnTheta),nii);
                ftilde = fft(bitrevorder(fx(1:nii))); %/nii;
                br_xun = bitrevorder(xun(1:nii,:));
                
                tic
                %par
                parfor k=1:numel(lnTheta)
                    [costMLE(iter,k),eigvalK(k,:)] = MLEKernel(obj, exp(lnTheta(k)),...
                        br_xun,ftilde,obj.order,obj.arbMean);
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
            set(hFigCost, 'units', 'inches', 'Position', [4 4 13.5 11.5])
            %title(lgd,'Sample Size, \(n\)'); legend boxoff
            xlabel('Shape param, \(\theta\)')
            % ylabel('MLE Cost, \( \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
            ylabel('Log MLE Obj. fun.')
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
            semilogx(obj.aMLEAll, obj.lossMLEAll, '+');
            temp = string(obj.mvec);
            temp(end+1) = '\(\theta_{min_{true}}\)';
            temp(end+1) = '\(\theta_{min_{est}}\)';
            legend(temp,'location','best'); axis tight
            saveas(hFigCost, plotFileName)
            
        end        

       
        % MLE objective function to find the optimal shape parmaeter
        function [loss,Ktilde,Kthat_new,RKHSnorm,K] = MLEKernel(obj, a,xun,ftilde)
            
            n = length(ftilde);
            if obj.order==4
                [K, Ktilde, Kthat_new] = cubMLELattice.kernel(xun,obj.order,a);
            elseif obj.order==2
                [K, Ktilde, Kthat_new] = cubMLELattice.kernel(xun,obj.order,a);
            else
                error('Unsupported Bernoulli polyn order !');
            end
            
            ftilde = abs(ftilde);  % remove any negative values
            
            Ktilde(1) = Kthat_new(1) + n;
            Ktilde(2:end) = Kthat_new(2:end);
            
            %RKHSnorm = mean(abs(ftilde).^2./Ktilde);
            
            % temp = (abs(ftilde(KtildeSq~=0))./(KtildeSq(KtildeSq~=0))).^2 ;
            temp = (abs(ftilde(Ktilde~=0).^2)./(Ktilde(Ktilde~=0))) ;
            
            if obj.arbMean==true
                RKHSnorm = sum(temp(2:end))/n;
                temp_1 = sum(temp(2:end));
            else
                RKHSnorm = sum(temp)/n;
                temp_1 = sum(temp);
            end
            cubMLESobol.alertMsg(temp_1, 'Imag');
            RKHSnormSq = sqrt(RKHSnorm);
            
            cubMLESobol.alertMsg(RKHSnormSq, 'Nan');
            % loss = mean(2*log(KtildeSq)) + log(RKHSnorm);
            % loss = mean(log(Ktilde)) + log(RKHSnorm);
            
            loss1 = sum(log(Ktilde(Ktilde~=0)));
            cubMLESobol.alertMsg(loss1, 'Inf');
            loss2 = n*log(temp_1);
            cubMLESobol.alertMsg(loss2, 'Inf');
            % ignore all zero val eigenvalues
            loss = sum(log(Ktilde(Ktilde~=0))) + n*log(temp_1);
            
            cubMLESobol.alertMsg(loss, 'Inf', 'Imag', 'Nan');
            cubMLESobol.alertMsg(Ktilde, 'Imag');
        end
    end
    
    methods(Static)  
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

        % computes Ktilde = K - 1
        % to avoid cancellation error in the computation of (1 - n/\lambda_1)
        function [Kt, K] = kernel_t(a, const, Bern)
            theta = a*const;
            d = size(Bern, 2);
            
            Kjt = theta*Bern(:,1);
            Kj = 1 + Kjt;
            
            for j=2:d
                Kjt_1 = Kjt; Kj_1 = Kj;
                
                Kjt = theta*Bern(:,j).*Kj_1 + Kjt_1;
                Kj = 1 + Kjt;
            end
            
            Kt = Kjt; K = Kt;
        end
        
        
        % bernoulli polynomial based kernel
        function [K, Ktilde, Kthat_new] = kernel(xun,order,a)
            
            constMult = -(-1)^(order/2)*((2*pi)^order)/factorial(order);
            if order == 2
                bernPloy = @(x)(-x.*(1-x) + 1/6);
            elseif order == 4
                bernPloy = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
            else
                error('Bernoulli order not implemented !');
            end
            K = prod(1 + (a)*constMult*bernPloy(xun),2);
            
            [Kt_new, K_new] = cubMLELattice.kernel_t(a, constMult, bernPloy(xun));
            
            Kthat_new = abs(fft(Kt_new));
            Khat_new = abs(fft(K_new));
            
            % matlab's builtin fft is much faster and accurate
            Ktilde = abs(fft(K));  % remove any negative values
            
            if sum(K)==length(K) || Ktilde(1)==length(K)
                %fprintf('debug');
            end
        end
        
        % computes the periodization transform for the given function values
        function f = doPeriodTx(fInput, ptransform)
            
            if strcmp(ptransform,'Baker')
                f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
            elseif strcmp(ptransform,'C0')
                f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
            elseif strcmp(ptransform,'C1')
                f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
            elseif strcmp(ptransform,'C1sin')
                f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(2*sin(pi*x).^2,2); % Sidi C^1 transform
            elseif strcmp(ptransform,'C2sin')
                psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16; psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
                f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
            elseif strcmp(ptransform,'C3sin')
                psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi); psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
                f=@(x) fInput(psi4(x)).*prod(psi4_1(x),2);
            elseif strcmp(ptransform,'none')
                % do nothing
                f=@(x) fInput(x);
            else
                error('Error: Periodization transform %s not implemented', ptransform);
            end
            
        end
        
        function xlat = simple_lattice_gen(n,d,firstBatch)
            if d<=10
                % this gives best accuracy
                z = [1, 364981, 245389, 97823, 488939, 62609, 400749, 385317, 21281, 223487]; % generator from Hickernell's paper
                %z = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599]; %generator
            else
                z = [1 182667 302247 433461 160317 94461 481331 252345 358305 221771 48157 489023 438503 399693 200585 169833 308325 247437 281713 424209 244841 205461 336811 359375 86263 370621 422443 284811 231547 360239 505287 355195 52937 344561 286935 312429 513879 171905 50603 441451 164379 139609 371213 152351 138607 441127 157037 510073 281681 380297 208143 497641 482925 233389 238553 121499 137783 463115 168681 70699];
            end
            
            z = z(1:d);
            
            if false
                % this is very slow
                if firstBatch==true
                    brIndices = bitrevorder((0:1/n:1-1/n));
                else
                    brIndices = bitrevorder((1/n:2/n:1-1/n));
                end
            else
                nmax = n;
                nmin = 1 + n/2;
                if firstBatch==true
                    nmin = 1;
                end
                nelem=nmax-nmin+1;
                
                if firstBatch==true
                    brIndices=cubMLELattice.vdc(nelem)';
                else
                    brIndices=cubMLELattice.vdc(nelem)'+1/(2*(nmin-1));
                end
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
        function y = fft_DIT( y )
            nmmin = log2(length(y));
            %y = bitrevorder(y);
            for l=0:nmmin-1
                nl=2^l;
                nmminlm1=2^(nmmin-l-1);
                ptind=repmat([true(nl,1); false(nl,1)],nmminlm1,1);
                coef=exp(-2*pi()*sqrt(-1)*(0:nl-1)'/(2*nl));
                coefv=repmat(coef,nmminlm1,1);
                evenval=y(ptind);
                oddval=y(~ptind);
                y(ptind)=(evenval+coefv.*oddval)/2;
                y(~ptind)=(evenval-coefv.*oddval)/2;
            end
        end
        
        
    end
end