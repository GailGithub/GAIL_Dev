function TractAve()
kMax = 10; N = 2.^(2:kMax); % number of data points
M = 7; % number of t points
nTest = 1000;
nFunc = 10;
dMax = 10;
AverageError = zeros(kMax-1,dMax);
SqrtNTest = sqrt(nTest);

%% kernel
alpha0 = 1; gamma0 = 1;
rAlpha = 2; rGamma = 2;
alpha = alpha0*(1:dMax)'.^(-rAlpha);
gamma = gamma0*(1:dMax)'.^(-rGamma);
kernel = @(a,g,x,t) prod(1-a+a.*exp(-g.*(x-t).^2));

%% main program
for k = 1:kMax-1
    for d = 1:dMax
        P = sobolset(d);
        tsites = net(scramble(P,'MatousekAffineOwen'),M);
        kernelD = @(x,t) kernel(alpha(1:d),gamma(1:d),x,t);
        KtMat = symKernelMat(kernelD,tsites);
        [V,sig,~] = svd(KtMat);
        dsig = diag(sig);
        lastBigSig = find(dsig>1e-10,1,'last');
        B = bsxfun(@times,V(:,1:lastBigSig),1./sqrt(dsig(1:lastBigSig)'));
        dsites = net(scramble(P,'MatousekAffineOwen'),N(k));
        ctrs = dsites;
        epoints = net(scramble(P,'MatousekAffineOwen'),nTest);
        IM = symKernelMat(kernelD,dsites);
        EM = kernelMat(kernelD,epoints,ctrs);
        for j = 1:nFunc
            coef = B*randn(lastBigSig,1);
            frandJ = @(x) frand(coef,kernelD,tsites,x);
            rhs = fctVec(frandJ,dsites);
            Pf = EM * (IM\rhs);
            exact = fctVec(frandJ,epoints);
            AverageError(k,d) = AverageError(k,d)+norm(Pf-exact);
        end
        AverageError(k,d) = AverageError(k,d)/SqrtNTest/nFunc;
    end
end
%% plot
fontSize = 16;
loglog(N,AverageError);
hold on;
plot(N,4*N.^(-sqrt(rAlpha*rGamma)),'m--o');
xlabel('Number of Data Points','FontSize',fontSize);
ylabel('RMSE','FontSize',fontSize);
%title('Interpolation With General Gaussian (Average Case)','FontSize',14);
myLegend = cell(dMax+1,1);
for d = 1:dMax
    myLegend(d) = cellstr(['d = ',num2str(d)]);
end
myLegend(dMax+1) = cellstr('Theoretical Convergence Rate');
legend(myLegend,'Location','SouthWest');
hFig = figure(1);
set(hFig, 'Position', [0 0 1080 720]);
end

function fval = frand(coef,kernel,tsites,x)
n = length(coef);
fval = 0;
for j = 1:n
    fval = fval+coef(j)*kernel(x,tsites(j,:)');
end
end