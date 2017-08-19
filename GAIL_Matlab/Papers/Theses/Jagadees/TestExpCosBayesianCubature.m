%% Test out MLE
function [nvec,muhat,aMLE,errCubMLE,out] = TestExpCosBayesianCubature(dim,BernPolyOrder,ptransform,figSavePath,visiblePlot,arbMean)

f = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
fName = 'Exp(cos)';
%figSavePath = '/home/jagadees/MyWriteup/Apr1stweek/';


fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

tic
    absTol = 1E-15;
    relTol = 0;
    order = BernPolyOrder;
    testAll = true;  % to plot the error, run for for all n values
    [muhatFinal,out]=cubMLELattice(f,dim,absTol,relTol,order,ptransform,testAll,fullPath,fName,arbMean);
    % function [muhat,out]=cubBayesLatticeRef(f,d,absTol,relTol,order,arbMean)
    % arbMean = true;
    % [muhatFinal,out]=cubBayesLatticeRef(f,dim,absTol,relTol,order,arbMean,ptransform,testAll,figSavePath,fName);
    nvec = 2.^out.mvec;
    muhat = out.muhatAll;
    ErrBd = out.ErrBdAll;
toc

%% plot error
errCubMLE = abs(exactInteg - muhat);
plotCubatureError(dim, nvec, errCubMLE, ErrBd, fName, BernPolyOrder, ptransform, ...
    fullPath,visiblePlot,arbMean, out.s_All, out.dscAll)

figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
        fullPath, fName, dim, BernPolyOrder, ptransform);
plot_nvec_vs_computeTime(nvec, out.timeAll, visiblePlot, figSavePathName)

fprintf('done');


end
