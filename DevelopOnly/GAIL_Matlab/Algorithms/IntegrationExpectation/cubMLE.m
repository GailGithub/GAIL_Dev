function [muhat,out]=cubMLE(fInput,nvec,domain,whSample,whKer,powerFuncMethod,BernPolynOrder,ptransform,fName,figSavePath)
%CUBMLE Monte Carlo method to estimate the mean of a random variable
%
%   tmu = cubML(Yrand,absTol,relTol,alpha,nSig,inflate) estimates the mean,
%   mu, of a f(X) using nvec samples of a random variable X in [0,1]^d.
%   The samples may be of one of several kinds.  The default values are n=2^10 and
%   d = 1 Input f is a function handle that accepts an n x d matrix of
%   n points in [0,1]^d and returns an n x 1 vector of f values.
%   powerFuncMethod:
%      'cauchy' : using the Cauchy interlacing theorem
%      'thompson' : using the technique form R.C.Thompson paper
%
% This is a heuristic algorithm based on a Central Limit Theorem
% approximation
if nargin < 7
    if ~exist('powerFuncMethod','var') || isempty(powerFuncMethod)
        powerFuncMethod='Cauchy';  %technique to compute power function%
    end
    if nargin < 6
        if nargin < 5
            whKer = 'Mat1'; %type of kernel
            if nargin < 4
                whSample = 'Sobol'; %type of sampling, scrambled Sobol
                if nargin < 3
                    domain = [0;1]; %dimension
                    if nargin < 2
                        nvec = 2^10; %number of samples
                        if nargin < 1
                            f = @(x) x.^2; %function
                        end
                    end
                end
            end
        end
    end
end

nmax = max(nvec);
nn = numel(nvec);
d = size(domain,2);
if strcmp(whSample,'Sobol')
    x = net(scramble(sobolset(d),'MatousekAffineOwen'),nmax);
end


%% Generate Rank-1 Lattice points
if strcmp(whSample,'Lattice1')
    x = gail.lattice_gen(1,nmax,d);
    %latticeseq_b2('init0');
    %x = latticeseq_b2(d,nmax)';
end


%% Periodization transform
%ptransform = 'C1sin'; %default option

if strcmp(ptransform,'Baker')
    f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(ptransform,'C0')
    f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(ptransform,'C1')
    f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2); % C^1 transform
elseif strcmp(ptransform,'C1sin')
    %f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(1-cos(2*pi*x),2); % Sidi C^1 transform
    % 1-cos(2x) = 2 sin^2(x)
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


yInput = f(x);
%% precompute the Bernoulli polynominal values to speedup the computation
out_aMLE(nn,1) = 0;
out_costMLE(nn,1) = 0;
muhat(nn,1) = 0;
out_disc2(nn,1) = 0;
out_ErrBd(nn,1) = 0;


if strcmp(whKer,'Fourier')
    xdiff = abs(x - x(1,:));
    if BernPolynOrder==2
        BernPolynX = (6*xdiff.*(xdiff - 1) + 1)/6;
    elseif BernPolynOrder==4
        %BernPolynX = (30*(xdiff.^2.*(xdiff.^2 - 2*xdiff + 1)) - 1)/30;
        %BernPolynX = xdiff.^4 - 2*xdiff.^3 + xdiff.^2 -1/30;
        BernPolynX = (xdiff.^2).*(xdiff.^2 - 2*xdiff + 1) -(1.0/30);
    elseif BernPolynOrder==6
        BernPolynX = x.^6 - 3*x.^5 + (5/2)*x.^4 - (1/2)*x.^2 + (1/42);
    else
        error('Error: Bernoulli polynomial order %d not implemented, exiting!!', BernPolynOrder);
    end
end

if  true %strcmp(whKer,'Fourier') % disabled in release code
    enableMovie = false;
    %% plot MLEKernel cost function
    lnTheta = -3:0.05:0;
    plotFileName = sprintf('%s%s Cost d_%d bernoulli_%d Period_%s', figSavePath, fName, d, BernPolynOrder, ptransform);
    if enableMovie
        figH = figure(21);
        set(figH, 'position', [10 10 1920 1080]);
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        
        loops = nn*numel(lnTheta);
        Frames(loops) = struct('cdata',[],'colormap',[]);
        warning('off','MATLAB:handle_graphics:MException:SceneNode');
        warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
        
        %set(figH, 'Visible', 'off');
        %set(gcf,'Visible','off')              % turns current figure "off"
        %set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"
        FramesT = reshape(Frames, nn, numel(lnTheta)); %zeros(nn, numel(lnTheta));
    end
    
    costMLE = zeros(nn,numel(lnTheta));
    tic
    
    for ii = 1:nn
        nii = nvec(ii)
        xReord = x(1:nii,:); yReord = yInput(1:nii); BernPolynXr = BernPolynX(1:nii,:);
        if true
            %[~,w2] = sort(xReord(:,1));
            w2 = bitrevorder(1:nii);
            xReord = xReord(w2,:); yReord = yReord(w2); BernPolynXr = BernPolynXr(w2,:); %reorder to get symmetric matrix
        end
        eigvalK = zeros(numel(lnTheta),nii); ffy = zeros(numel(lnTheta),nii);
        tic
        %par
        for k=1:numel(lnTheta)
            [costMLE(ii,k),eigvalK(k,:),ffy(k,:)] = MLEKernel(exp(lnTheta(k)),xReord,yReord,whKer,domain,BernPolynXr,BernPolynOrder);
            %costMLE(ii,k) = MLEKernel(exp(lnTheta(k)),xReord,yReord,whKer,domain,BernPolynXr,BernPolynOrder);
        end
        
        if enableMovie
            for k=1:numel(lnTheta)
                
                warning('off','MATLAB:handle_graphics:MException:SceneNode');
                warning('off','MATLAB:handle_graphics:exceptions:SceneNode')
                %drawnow
                create_cost_plots(eigvalK(k),ffy(k),exp(lnTheta),k,costMLE(ii,:));
                set(figH, 'position', [10 10 1920 1080]);
                %Frames((ii-1)*numel(lnTheta) + k) = getframe(figH);
                FramesT(ii,k) = getframe(figH);
            end
        end
        
        toc
    end
    
    toc
    if enableMovie
        Frames1 = reshape(FramesT',1,nn*numel(lnTheta));
        myVideo = VideoWriter(strcat(plotFileName, '_movie3.avi'), 'Motion JPEG AVI');
        open(myVideo);
        writeVideo(myVideo,Frames1)
        close(myVideo)
    end
    
    hFigCost = figure; semilogx(exp(lnTheta),costMLE); %lgd = legend(string(nvec),'location','north'); axis tight
    set(hFigCost, 'units', 'inches', 'Position', [4 4 10 7])
    %title(lgd,'Sample Size, \(n\)'); legend boxoff
    xlabel('Shape param, \(\theta\)')
    ylabel('MLE Cost, \( \frac{y^T K_\theta^{-1}y}{[\det(K_\theta^{-1})]^{1/n}} \)')
    axis tight;
    title(sprintf('%s Cost d=%d Bernoulli=%d PeriodTx=%s', fName, d, BernPolynOrder, ptransform));
    [minVal,Index] = min(costMLE,[],2);
    hold on; plot(exp(lnTheta(Index)),minVal, '.');
    
    minTheta = exp(lnTheta(Index));
end


tstart = tic; %start the clock
parfor ii = 1:nn
    nii = nvec(ii);
    
    xReord = x(1:nii,:); yReord = yInput(1:nii); BernPolynXr = BernPolynX(1:nii,:); % default Initial values
    if true % disable %strcmp(whKer,'Fourier')
        %% Reorder the x, yInput to get symmetric circulant (Toeplitz) matrix
        w2 = bitrevorder(1:nii);
        %[~,w2] = sort(xReord(:,d-1));
        xReord = xReord(w2,:);
        yReord = yReord(w2);
        BernPolynXr = BernPolynX(w2,:);
    end
    if strcmp(whKer,'Fourier')
        ax = -3; bx = 0; % limit theta within (0,1]
    else
        ax = -5; bx = 5;
    end
    
    %% Estimate optimal \theta
    [lnaMLE, costMin] = fminbnd(@(lna) ...
        MLEKernel(exp(lna),xReord,yReord,whKer,domain,BernPolynXr,BernPolynOrder), ...
        ax,bx,optimset('TolX',1e-4));
    if 0
        costFunc = @(lna)MLEKernel(exp(lna),xReord,yReord,whKer,domain,BernPolynXr,BernPolynOrder);
        [xmin, costMin] = fminsearch(costFunc,0.5*ones(1,d) );
        aMLE = exp(xmin);
    end
    aMLE = exp(lnaMLE);
    
    %% Use the optimal \theta
    out_aMLE(ii,:) = aMLE;
    out_costMLE(ii) = costMin;
    [K,kvec,k0] = kernelFun(xReord,whKer,aMLE,domain,BernPolynXr,BernPolynOrder);
    
    if strcmp(whKer,'Fourier')
        cn = K;
        
        ffy = abs(fft(yReord-1))/nii;
        ffy(1)=1;
        eigval = fft(cn'-1)/nii; % subtract 1 and then after fft computation, put it back in fft
        eigval(1) = 1;  % this helps to avoid zero values in fft
        if any(eigval==0)
            fprintf('Zero eigval in estimating Kinvy \n')
        end
        eigvaln = eigval + (eigval==0)*eps;
        Kinvy = ifft(ffy./eigvaln);
        
        %% compute the approximate mu
        muhat(ii) = sum(yReord)/sum(cn);
        muhat_old = real(kvec'*Kinvy);
        
        % [muhat_old muhat(ii) muhat_old-muhat(ii)]
    else
        Kinv = pinv(K);
        %w = Kinv*kvec;
        Kinvy = Kinv*yInput(1:nii);
        
        %% compute the approximate mu
        muhat(ii) = real(kvec'*Kinvy);
    end
    
    %out_Kcond(ii) = cond(K);
    
    %% compute the discriminant
    if strcmp(whKer,'Fourier') == true
        %disc2_old = k0 - kvec'*ifft(fft_DIT(kvec)./fft_DIT(cn'));
        disc2 = (k0*sum(cn) - nii)/sum(cn);
        %ffy = abs(fft(yReord))/nii;
        val2 = sum((ffy(eigval~=0).^2)./abs(eigval(eigval~=0)));
        out_ErrBd(ii) = 2.58*sqrt(disc2*val2/nii);
    else
        if strcmp(powerFuncMethod, 'Cauchy')
            eigK = eig(K);
            eigKaug = eig([k0 kvec'; kvec K]);
            disc2 = exp(sum(log(eigKaug(1:nii)) - log(eigK)))*eigKaug(end);
        else
            % from R.C.Thompson paper
            [eigVecKaug, eigValKaug] = eig([k0 kvec'; kvec K]);
            if isvector(eigValKaug) == false
                eigValKaug = diag(eigValKaug);
            end
            uii = eigVecKaug(:,1);
            disc2 = 1.0/sum(uii.^2 ./ eigValKaug);
        end
        out_ErrBd(ii) = 2.58*sqrt(disc2*(yReord'*Kinvy))/nii;
    end
    out_disc2(ii) = disc2;
    
    
end

% add the actual min theta and from fminbnd to the plot
if true %strcmp(whKer,'Fourier') == true
    if isvalid(hFigCost)
        lgdText = string(nvec);
        for i=1:length(nvec)
            lgdText(i) = sprintf('%-7s %.3f', lgdText(i), out_aMLE(i));
        end
        legend(lgdText,'location','north'); legend boxoff
        plot(out_aMLE,out_costMLE, 'c+');
        saveas(hFigCost, strcat(plotFileName, '.jpg'))
    end
end

out.aMLE = out_aMLE;
out.disc2 = out_disc2;
out.ErrBd = abs(out_ErrBd);
out.time = toc(tstart)
out.BernPolynOrder = BernPolynOrder;
out.ptransform = ptransform;

end


% debug function to create plots
function create_cost_plots(eigval,ffy,theta,k,costMLE)
shape = theta(k);
%figure(21);
[~, I] = sort(eigval,'descend');
zi2 = ffy.^2;
temp = zi2(I)./eigval(I);
%thresh = quantile(temp, 0.95);
%temp(temp > thresh) = 0;
cumz = cumsum(temp);
subplot(2,3,1); loglog(cumz)
xlabel('count')
ylabel('cumsum of \(\frac{z_i^2}{\gamma_i}\)'); axis tight manual
title(sprintf('shape %0.3f, n %d', shape, length(ffy)))

subplot(2,3,2); loglog(eigval(eigval~=0), (ffy(eigval~=0).^2)./abs(eigval(eigval~=0)), '.')
xlabel('Eigvals of K, \(\gamma\)')
ylabel('DFT(y)**2/eigvals, \(\frac{z_i^2}{\gamma_i}\)'); axis tight manual

subplot(2,3,3); loglog(eigval(eigval~=0), ffy(eigval~=0).^2, '.');  xlabel('Eigvals of K, \(\gamma\)'); ylabel('DFT(y)**2, \(z_i^2\)'); axis tight manual
title('');

%subplot(2,3,4); loglog(eigval(eigval~=0)); xlabel('count'); ylabel('Eigvals of K, \(\gamma\)'); axis tight manual
subplot(2,3,4); loglog(eigval(I)); xlabel('count'); ylabel('Eigvals of K, \(\gamma\)'); axis tight manual

%subplot(2,3,5); loglog(ffy(eigval~=0)); xlabel('count'); ylabel('DFT of y, \(z_i\)'); axis tight manual
subplot(2,3,5); loglog(ffy(I)); xlabel('count'); ylabel('DFT of y, \(z_i\)'); axis tight manual
subplot(2,3,6); semilogx(theta,costMLE); xlabel('theta'); ylabel('cost'); axis tight manual
end


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
