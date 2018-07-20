%% Test Multivariate Normal Probabilities
%function [muhat,err,time,outVec] = TestMVN_BayesianCubature(dim,BernPolyOrder,...
%  ptransform,figSavePath,visiblePlot,arbMean,stopAtTol,samplingMethod,absTol,relTol)
function [muhat,err,time,outVec] = TestMVN_BayesianCubature(varargin)


dim = get_arg('dim', varargin);
ptransform = get_arg('ptransform', varargin);
stopAtTol = get_arg('stopAtTol', varargin);
figSavePath = get_arg('figSavePath', varargin);
visiblePlot = get_arg('visiblePlot', varargin);
samplingMethod = get_arg('samplingMethod', varargin);
BernPolyOrder = get_arg('order', varargin);
arbMean = get_arg('arbMean', varargin);

if dim==2
  % d-3 problem reduced to d-2 using Genz method
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
    %disp('Already have gold standard answer')
    compGold = false;
  end
end
if compGold
  disp('(Re-)computing gold standard answer')
  % limit to 2 parallel threads when memory is limited
  % parpool('local',2)  
  muBestvec = zeros(1,nRepGold);
  tic
  %par
  parfor i = 1:nRepGold
    i
    tic
    muBestvec(1,i) = compProb(MVNProbBest);
    toc
  end
  toc
  muBest = mean(muBestvec);
end
%disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])
% muBest = '0.676337324357787483819492990733124315738677978515625';

%% Try MLE Bayseian cubature with Fourier kernel and Rank1 Lattice points
nvecMLE = 2.^(10:20)';
nnMLE = numel(nvecMLE);

if exist('absTol','var')==false
  absTolVal = 1e-2;
else
  absTolVal = absTol;
end

if exist('relTol','var')==false
  relTolVal = 1e-2;
else
  relTolVal = relTol;
end

if strcmp(samplingMethod,'Sobol')
  cubMethod='MLESobol';
else
  cubMethod='MLELattice';
end

if exist('stopAtTol','var')==false
  stopAtTol = false;
end
MVNProbMLELatticeGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvecMLE, ...
  'errMeth','g','cubMeth',cubMethod,'intMeth','Genz', ...
  'BernPolyOrder',BernPolyOrder,'ptransform',ptransform, ...
  'fName',fName,'figSavePath',fullPath,'arbMean',arbMean,...
  'absTol',absTolVal, 'relTol',relTolVal, 'stopAtTol',stopAtTol);
compMLELattice = true;

integrand = @(t) Genz(t,MVNProbMLELatticeGn);

inputArgs = varargin;
inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;

% initialise the object based on the sampling method
if exist('samplingMethod','var') && ...
    strcmp(samplingMethod,'Sobol') % use Sobol points
  objCubMLE=cubMLESobol(inputArgs{:});
else % use Lattice points
  objCubMLE=cubMLELattice(inputArgs{:});
end

if compMLELattice
  tic
  
  nRep = 30; %100; % increase it for gail plots
  %outVec = {'a', 'b'};
  muMVNProbMLELatticeGn = zeros(nnMLE,nRep);
  aMLE = zeros(nnMLE,nRep);
  errbdvecMBVProbMLELatticeGn(nnMLE,nRep) = 0;
  muhatVec(nRep) = 0;
  %timeVec(nRep) = 0;
  nPointsVec(nRep) = 0;
  for i = 1:nRep
    %if i/1 == floor(i/1), i, end
    %[muhat, out] = compProb(MVNProbMLELatticeGn);
    [muhatVec(i),outVec(i)]=compInteg(objCubMLE);
    muMVNProbMLELatticeGn(:,i) = outVec(i).muhatAll;
    errbdvecMBVProbMLELatticeGn(:,i) = outVec(i).ErrBdAll;
    aMLE(:,i) = outVec(i).aMLEAll;
    nPointsVec(i) = outVec(i).n;
  end
  timeVec = [outVec(:).time];
  % loglog(2.^(out.mvec) , (abs(muBest - muMVNProbMLELatticeGn(:,1:i))), 2.^(out.mvec) , (abs(errbdvecMBVProbMLELatticeGn(:,1:i)))); axis tight
  
  errvecMVNProbMLELatticeGn = abs(muBest - muMVNProbMLELatticeGn);
  errCubMLE = median(errvecMVNProbMLELatticeGn,2);
  errVec = abs(muBest - muhatVec);
  errtopMVNProbMLELatticeGn = quantile(errvecMVNProbMLELatticeGn,1-alpha,2);
  ErrBd = quantile(errbdvecMBVProbMLELatticeGn,1-alpha,2);
  fprintf('\nMedian error_n %1.2e, worst_n %d, worst_time %1.3f, absTol %1.2e, relTol %1.2e ', ...
    median(errVec), quantile(nPointsVec,1-alpha), quantile(timeVec,1-alpha), ...
    outVec(1).absTol, outVec(1).relTol);
  
  %time = quantile(timeVec,1-alpha);
  time = median([outVec(:).time]);
  err = median(errVec);
  muhat = median(muhatVec);
  if err/outVec(1).absTol > 1
    error('wait')
  end
  toc
  
  if exist('stopAtTol','var')==false || stopAtTol==false
    nvecMLE = 2.^outVec(1).mvec;
    plotCubatureError(dim, nvecMLE, errCubMLE, ErrBd, fName, BernPolyOrder, ...
      ptransform, fullPath,visiblePlot,arbMean,outVec{1}.s_All, outVec{1}.dscAll)
    figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
      fullPath, fName, dim, BernPolyOrder, ptransform);
    plot_nvec_vs_computeTime(nvecMLE, outVec{1}.timeAll, visiblePlot, figSavePathName, samplingMethod)
  end
  
  fprintf('done\n');
  
end

end
%% Save output
%save MVNProbExampleAllData.mat


% picks the input argument from the varargin cell array
function output = get_arg(argName, inputArgs, defaultVal)
iStart = 1;
wh = find(strcmp(inputArgs(iStart:end),argName));
if ~isempty(wh), output=inputArgs{wh+iStart}; else, output=defaultVal; end
end


