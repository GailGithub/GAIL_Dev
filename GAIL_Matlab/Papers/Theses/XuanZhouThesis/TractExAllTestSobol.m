% DistanceMatrixFit
% Script that uses Euclidean distance matrices to perform
% scattered data interpolation for arbitrary space dimensions
% Calls on: DistanceMatrixB, CreatePoints, testfunctionsD,
%           PlotSurf, PlotError2D, PlotSlices, PlotErrorSlices
% Uses:     various routines called by CreatePoints
clear all; close all;
kmax = 12;
d = 10;
nTestFun = 5;
ep = 1; r = 2;
N = 2.^(1:kmax);
M = 1000;
RMS_err = zeros(kmax,d,nTestFun);
rGenz = zeros(nTestFun,1);
rGenz(1) = 2;
rGenz(2) = .3;
rGenz(3) = 2;
rGenz(4) = 2;
rGenz(5) = 2;
%rGenz(6) = 2;
alpha = zeros(nTestFun,d);
alpha(1,:) = (1:d).^(-rGenz(1));
alpha(2,:) = (1:d).^rGenz(2);
alpha(3,:) = (1:d).^(-rGenz(3));
alpha(4,:) = (1:d).^(-rGenz(4));
alpha(5,:) = (1:d).^(-rGenz(5));
%alpha(6,:) = (1:d).^(-rGenz(6));
beta = zeros(nTestFun,d);
beta(2,:) = linspace(0,1,d);
beta(4,:) = linspace(0,1,d);
beta(5,:) = linspace(0,1,d);
%beta(6,:) = rand(1,d);
sqrtM = sqrt(M);
for k = 1:kmax
    for s = 1:d
        p = sobolset(s);
        dsites = net(p,N(k));
        ctrs = dsites;
        epoints = CreatePoints(M,s,'u');
        IM = GaussHermite(ep,r,dsites,ctrs);
        EM = GaussHermite(ep,r,epoints,ctrs);
        for index = 1:nTestFun
            rhs = genz_test_fun(dsites,index,s,alpha(index,:),beta(index,:),rGenz(index));
            Pf = EM * (IM\rhs);
            exact = genz_test_fun(epoints,index,s,alpha(index,:),beta(index,:),rGenz(index));
            RMS_err(k,s,index) = norm(Pf-exact)/sqrtM;
        end
    end
end
%save('TractExAllTestSobol.mat','RMS_err','ep','r','rGenz','alpha','beta');
figure;
for index = 1:nTestFun
    subplot(2,3,index);
    loglog(N,RMS_err(:,:,index));
    xlabel('n','FontSize',12);
    ylabel('RMSE','FontSize',12);
    title(['Interpolation with Anisotropic Gaussian (Genz ',num2str(index),')'],'FontSize',14);
end
myLegend = cell(d,1);
for s = 1:d
    myLegend(s) = cellstr(['d = ',num2str(s)]);
end
legend(myLegend);
