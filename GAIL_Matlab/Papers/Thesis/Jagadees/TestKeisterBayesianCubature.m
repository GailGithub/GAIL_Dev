%% Keister's Example of Multidimensional Integration
%
% B. D. Keister, Multidimensional quadrature algorithms, _Computers in
% Physics_, *10*, pp. 119-122, 1996, presents the following
% multidimensional integral, inspired by a physics application:
%
% \[ I = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, \mathrm{d} \boldsymbol{x},
% \qquad d = 1, 2, \ldots. \]

function [muhat,aMLE,err,out] = TestKeisterBayesianCubature(dim,BernPolyOrder,...
  ptransform,figSavePath,visiblePlot,arbMean,stopAtTol)

normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t

domain = repmat([0;1],[1,dim]);

replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));  %using erfcinv is more accurate than erfinv with -1
fKeister = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;


%% Bayesian Cubature
fName='Keister';

f1 = @(x) fKeister(x,dim);

fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

tic
    absTol = 1E-15;
    relTol = 0;
    order = BernPolyOrder;
    if ~exist('stopAtTol','var')
      stopAtTol = false;  % to plot the error, run for for all n values
    end
    tic
    [muhatFinal,out]=cubMLELattice(f1,dim,absTol,relTol,order,ptransform,stopAtTol,fullPath,fName,arbMean);
    toc
    nvec = 2.^out.mvec;
    muhat = out.muhatAll;
    ErrBd = out.ErrBdAll;
toc


%% plot error
exactInteg = Keistertrue(dim);
errCubMLE = abs(exactInteg - muhat);

plotCubatureError(dim, nvec, errCubMLE, ErrBd, fName, BernPolyOrder, ptransform, ...
    fullPath,visiblePlot,arbMean, out.s_All, out.dscAll)

figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
        fullPath, fName, dim, BernPolyOrder, ptransform);
plot_nvec_vs_computeTime(nvec, out.timeAll, visiblePlot, figSavePathName)

fprintf('Done\n')
