function [tmu,out_param]=meanMCCV_g(YXrand,meanX,absTol,relTol)
% meanMCCV_g is a function which can estimate the mean of a random
% variable Y, while using another random variable X to generate ycv to 
% replace the mean of Y, which has the same mean as Y but with lower
% variance. Then we will need fewer samples, and save a lot of time.
%
% YXrand is the function generating a n*(p+1) random variable matrix, the
% first column is the sample of Y, and the other p columns represent X
%
% meanX is the known true means of X, with size of 1*p
% absTol and relTol are just the same as in meanMC_g
[tmu,out_param]=meanMC_g(@(n)controlvariate(YXrand,n,meanX),absTol,relTol);
end

function ycv = controlvariate(YXrand,n,meanX)
% this is the control variate part, follows the formula of control variate
YX = YXrand(n);% generate the matrix YX
Y =YX(:,1); % get the first column
X =YX(:,2:end); % get the rest columns
beta = bsxfun(@minus,X,mean(X,1))\(Y-mean(Y)); % the optimal beta estimated
ycv = Y - bsxfun(@minus,X,meanX)*beta; % ycv is the control variate random variable
end