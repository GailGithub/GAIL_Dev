%% Test out MLE
function [nvec,muhat,aMLE,errCubMLE,out] = TestExpCos(dim,BernPolyOrder,ptransform,figSavePath,visiblePlot)

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
    [muhatFinal,out]=cubMLELattice(f,dim,absTol,relTol,order,ptransform,testAll,figSavePath,fName);
    nvec = 2.^out.mvec;
    muhat = out.muhatAll;
    ErrBd = out.ErrBdAll;
toc

%% plot error
errCubMLE = abs(exactInteg - muhat);
plotCubatureError(dim, nvec, errCubMLE, ErrBd, fName, BernPolyOrder, ptransform, fullPath,visiblePlot)

fprintf('done');


end
