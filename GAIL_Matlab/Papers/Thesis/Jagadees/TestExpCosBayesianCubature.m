%% integrat Exp(cos) using Bayesian cubature
function [nvec,muhat,aMLE,errCubMLE,out] = TestExpCosBayesianCubature(dim,BernPolyOrder,...
  ptransform,figSavePath,visiblePlot,arbMean)

% define the integrand function
f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
fName = 'Exp(cos)';

% set the output dir to save the plots
fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

tic
    absTol = 1E-15;
    relTol = 0;
    order = BernPolyOrder;
    stopAtTol = false;  % to plot the error, run for for all n values
    [muhatFinal,out]=cubMLELattice(f,dim,absTol,relTol,order,ptransform,stopAtTol,fullPath,fName,arbMean);
    nvec = 2.^out.mvec;
    muhat = out.muhatAll;
    ErrBd = out.ErrBdAll;
toc

%% plot error
errCubMLE = abs(exactInteg - muhat);
plotCubatureError(dim, nvec, errCubMLE, ErrBd, fName, BernPolyOrder, ptransform, ...
    fullPath,visiblePlot,arbMean, out.s_All, out.dscAll)

% plot computation time vs number of points 
figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
        fullPath, fName, dim, BernPolyOrder, ptransform);
plot_nvec_vs_computeTime(nvec, out.timeAll, visiblePlot, figSavePathName)

fprintf('done');


end
