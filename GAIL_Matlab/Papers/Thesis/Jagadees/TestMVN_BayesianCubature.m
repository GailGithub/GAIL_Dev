%% Test Multivariate Normal Probabilities

function [muhat,aMLE,err,outAll] = TestMVN_BayesianCubature(dim,BernPolyOrder,...
          ptransform,figSavePath,visiblePlot,arbMean)

%gail.InitializeWorkspaceDisplay %clean up
%format long

nvec = 2.^(7:20)';
nlarge = nvec(end)*2;
nn = numel(nvec);
mu = 0;

if dim==2
    C = [4 1 1; 0 1 0.5; 0 0 0.25];
    Cov = C'*C;
    a = [-6 -2 -2];
    b = [5 2 1];
elseif dim==3
    C = [4 1 1 1; 0 1 0.5 0.5; 0 0 0.25 0.25; 0 0 0 0.25];
    Cov = C'*C;
    a = [-6 -2 -2 -2];
    b = [5 2 1 2];
else
    error('wrong dimension')
end
nRep = 100;
alpha = 0.1;

if exist('MVNProbExampleAllData.mat','file')
    load MVNProbExampleAllData
    MVNProbBestArch = MVNProbBest;
    nRepGoldArch = nRepGold;
    nRepArch = nRep;
    MVNProbIIDGnArch = MVNProbIIDGn;
    MVNProbSobolGnArch = MVNProbSobolGn;
    MVNProbuSobolGnArch = MVNProbuSobolGn;
    MVNProbMLESobolGnArch = MVNProbMLESobolGn;
    if exist('MVNProbMLELatticeGn', 'var')
        MVNProbMLELatticeGnArch = MVNProbMLELatticeGn;
    end
end

fName = 'MVN' ;
%figSavePath = '/home/jagadees/MyWriteup/Apr1stweek/';
fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
    mkdir(fullPath);
end

%% First compute a high accuracy answer
nGold = 2^27;
nRepGold = nRep;
MVNProbBest = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nGold, ...
    'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compGold = true;
if exist('MVNProbExampleAllData.mat','file')
    if sameProblem(MVNProbBest,MVNProbBestArch) && ...
            nRepGoldArch == nRepGold
        disp('Already have gold standard answer')
        compGold = false;
    end
end
if compGold
    disp('(Re-)computing gold standard answer')
    muBestvec = zeros(1,nRepGold);
    tic
    parfor i = 1:nRepGold
        i
        tic
        muBestvec(1,i) = compProb(MVNProbBest);
        toc
    end
    toc
    muBest = mean(muBestvec);
end
disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])

%% Try MLE Bayseian cubature with Fourier kernel and Rank1 Lattice points
nvecMLE = 2.^(10:23)';
nnMLE = numel(nvecMLE);
MVNProbMLELatticeGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvecMLE, ...
    'errMeth','n','cubMeth','MLELattice','intMeth','Genz', ...
    'BernPolyOrder',BernPolyOrder,'ptransform',ptransform, ...
    'fName',fName,'figSavePath',fullPath,'arbMean',arbMean);
compMLELattice = true;

if compMLELattice
    datetime
    tic
    nRep = 2; % reduced it
    outAll = {'a', 'b'};
    muMVNProbMLELatticeGn = zeros(nnMLE,nRep);
    aMLE = zeros(nnMLE,nRep);
    errbdvecMBVProbMLELatticeGn(nnMLE,nRep) = 0;
    for i = 1:nRep
        if i/1 == floor(i/1), i, end
        tic
        [~, out] = compProb(MVNProbMLELatticeGn);
        muMVNProbMLELatticeGn(:,i) = out.muhatAll;
        toc
        errbdvecMBVProbMLELatticeGn(:,i) = out.ErrBdAll;
        aMLE(:,i) = out.aMLEAll;
        outAll{i} = out;
    end
    muhat = muMVNProbMLELatticeGn;
    
    % loglog(2.^(out.mvec) , (abs(muBest - muMVNProbMLELatticeGn(:,1:i))), 2.^(out.mvec) , (abs(errbdvecMBVProbMLELatticeGn(:,1:i)))); axis tight
    
    errvecMVNProbMLELatticeGn = abs(muBest - muMVNProbMLELatticeGn);
    errCubMLE = median(errvecMVNProbMLELatticeGn,2);
    err = errCubMLE;
    errtopMVNProbMLELatticeGn = quantile(errvecMVNProbMLELatticeGn,1-alpha,2);
    ErrBd = quantile(errbdvecMBVProbMLELatticeGn,1-alpha,2);
    toc
    datetime

   
    plotCubatureError(dim, nvecMLE, errCubMLE, ErrBd, fName, BernPolyOrder, ...
      ptransform, fullPath,visiblePlot,arbMean,outAll{1}.s_All, outAll{1}.dscAll)

    figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
        fullPath, fName, dim, BernPolyOrder, ptransform);
    plot_nvec_vs_computeTime(nvecMLE, outAll{1}.timeAll, visiblePlot, figSavePathName)
    fprintf('done\n');

end


%% Save output
%save MVNProbExampleAllData.mat








