function [tmu,out_param]=meanMC_CV(YX,meanX,absTol,relTol,alpha)
% YOPTPRICE_CV creates the control variate output for option pricing using
% the |optPayoff| object

beta =@(YX) (bsxfun(@minus,YX(:,2:end),mean(YX(:,2:end),1))\YX(:,1)); %optimal beta
YCV =@(YX) (YX(:,1) - bsxfun(@minus,YX(:,2:end),meanX)*beta(YX)); %control variate random variablea
[tmu, out_param] = meanMC_g(@(n) (YCV(YX(n))), absTol, relTol,alpha);