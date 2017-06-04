% DistanceMatrixFit
% Script that uses Euclidean distance matrices to perform
% scattered data interpolation for arbitrary space dimensions
% Calls on: DistanceMatrixB, CreatePoints, testfunctionsD,
%           PlotSurf, PlotError2D, PlotSlices, PlotErrorSlices
% Uses:     various routines called by CreatePoints
clear all; close all;
kmax = 6;
ep = 1; r = 1;
N = 2.^(1:kmax);
M = 1000;
d = 10;
RMS_err = zeros(kmax,d);
alpha = (1:d).^(-r);
beta = (1:d).^(-r);
index = 1;
sqrtM = sqrt(M);
for k = 1:kmax
    for s = 1:d
        dsites = CreatePoints(N(k),s,'h');
        ctrs = dsites;
        epoints = CreatePoints(M,s,'u');
        rhs = genz_test_fun(dsites,index,s,alpha,beta,r);
        IM = GaussHermite(ep,r,dsites,ctrs);
        EM = GaussHermite(ep,r,epoints,ctrs);
        Pf = EM * (IM\rhs);
        exact = genz_test_fun(epoints,index,s,alpha,beta,r);
        RMS_err(k,s) = norm(Pf-exact)/sqrtM;
    end
end
loglog(N,RMS_err);
xlabel('n','FontSize',12);
ylabel('RMSE','FontSize',12);
title(['Interpolation with Anisotropic Gaussian (Genz ',num2str(index),')'],'FontSize',14);
myLegend = cell(d,1);
for s = 1:d
    myLegend(s) = cellstr(['d = ',num2str(s)]);
end
legend(myLegend);
print(gcf,'-depsc',['TractExGenz',num2str(index),'.eps']);