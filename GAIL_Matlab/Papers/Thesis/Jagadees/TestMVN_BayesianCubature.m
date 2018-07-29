%% Test Multivariate Normal Probabilities
%
function [muhat,errVec,timeVec,outVec] = TestMVN_BayesianCubature(varargin)


dim = get_arg('dim', varargin);
ptransform = get_arg('ptransform', varargin);
figSavePath = get_arg('figSavePath', varargin);
visiblePlot = get_arg('visiblePlot', varargin);
samplingMethod = get_arg('samplingMethod', varargin);
BernPolyOrder = get_arg('order', varargin);
arbMean = get_arg('arbMean', varargin);

if dim==2
  % d-3 problem reduced to d-2 using Genz method
  MVNParams.C = [4 1 1; 0 1 0.5; 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2];
  MVNParams.b = [5 2 1];
  MVNParams.mu = 0;
elseif dim==3
  MVNParams.C = [4 1 1 1; 0 1 0.5 0.5; 0 0 0.25 0.25; 0 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2 -2];
  MVNParams.b = [5 2 1 2];
  MVNParams.mu = 0;
else
  error('wrong dimension')
end
MVNParams.CovProp.C = chol(MVNParams.Cov)';

fName = 'MVN' ;
fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
  mkdir(fullPath);
end

alpha = 0.1;

muBest = computeGoldenMuhat(MVNParams);

%% Try MLE Bayseian cubature with Fourier kernel and Rank1 Lattice points

compMLELattice = true;

integrand = @(t) GenzFunc(t,MVNParams);

inputArgs = varargin;
inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;

% initialise the object based on the sampling method
if exist('samplingMethod','var') && ...
    strcmp(samplingMethod,'Sobol') % use Sobol points
  objCubBayes=cubMLESobol(inputArgs{:});
else % use Lattice points
  objCubBayes=cubBayesLattice_g(inputArgs{:});
end

nvecMLE = 2.^(objCubBayes.mvec');
nnMLE = numel(nvecMLE);

if compMLELattice
  tic
  
  nRep = 100; % increase it for gail plots
  muMVNProbMLELatticeGn = zeros(nnMLE,nRep);
  aMLE = zeros(nnMLE,nRep);
  errbdvecMBVProbMLELatticeGn(nnMLE,nRep) = 0;
  muhatVec(nRep,1) = 0;
  nPointsVec(nRep,1) = 0;
  for i = 1:nRep
    %if i/1 == floor(i/1), i, end
    [muhatVec(i),outVec(i)]=compInteg(objCubBayes);
    muMVNProbMLELatticeGn(:,i) = outVec(i).muhatAll;
    errbdvecMBVProbMLELatticeGn(:,i) = outVec(i).ErrBdAll;
    aMLE(:,i) = outVec(i).aMLEAll;
    nPointsVec(i) = outVec(i).n;
  end
  timeVec = [outVec(:).time]';
  
  errvecMVNProbMLELatticeGn = abs(muBest - muMVNProbMLELatticeGn);
  errCubMLE = median(errvecMVNProbMLELatticeGn,2);
  errVec = abs(muBest - muhatVec);
  %errtopMVNProbMLELatticeGn = quantile(errvecMVNProbMLELatticeGn,1-alpha,2);
  ErrBd = quantile(errbdvecMBVProbMLELatticeGn,1-alpha,2);
  fprintf('\nMedian error_n %1.2e, worst_n %d, worst_time %1.3f, absTol %1.2e, relTol %1.2e ', ...
    median(errVec), quantile(nPointsVec,1-alpha), quantile(timeVec,1-alpha), ...
    outVec(1).absTol, outVec(1).relTol);
  
  %time = quantile(timeVec,1-alpha);
  %time = median(timeVec);
  err = median(errVec);
  muhat = median(muhatVec);
  if err/outVec(1).absTol > 1
    error('wait')
  end
  toc
  
  if outVec(1).stopAtTol==false
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

function fval = GenzFunc(w,obj)
dim = numel(obj.a);
nn = size(w,1);
am = obj.a - obj.mu;
bm = obj.b - obj.mu;
a1 = am(1)/obj.CovProp.C(1,1);
b1 = bm(1)/obj.CovProp.C(1,1);
d = gail.stdnormcdf(a1);
e = gail.stdnormcdf(b1);
fval = (e-d)*ones(nn,1);
y = zeros(nn,dim-1);
for i = 2:dim
  y(:,i-1) = gail.stdnorminv(d+w(:,i-1).*(e-d));
  aux = sum(bsxfun(@times,obj.CovProp.C(i,1:i-1),y(:,1:i-1)),2);
  a1 = (am(i)-aux)/obj.CovProp.C(i,i);
  b1 = (bm(i)-aux)/obj.CovProp.C(i,i);
  d = gail.stdnormcdf(a1);
  e = gail.stdnormcdf(b1);
  fval = fval .* (e-d);
end
end

function muBest = computeGoldenMuhat(MVNParams)

if 0
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
end

if 0
%% First compute a high accuracy answer
nRep = 100;

nGold = 2^27;
nRepGold = nRep;
MVNProbBest = multivarGauss('a',MVNParams.a,'b',MVNParams.b,'Cov',MVNParams.Cov, ...
  'n',nGold,'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compGold = true;
if 0  %exist('MVNProbExampleAllData.mat','file')
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
  for i = 1:nRepGold
    i
    tic
    muBestvec(1,i) = compProb(MVNProbBest);
    toc
  end
  toc
  muBest = mean(muBestvec);
  
  save('MVNProbExampleDim4.mat','muBest', 'muBestvec', 'nRepGold', 'MVNProbBest');
end
%disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])
end

dim = length(MVNParams.a)-1;
if dim==2
  muBest = '0.676337324357787483819492990733124315738677978515625';
elseif dim==3
  muBest = '0.67451648307312195296248091835877858102321624755859375';
else
  error('Unsupported dim')
end
muBest =str2double(muBest);

end
    
% picks the input argument from the varargin cell array
function output = get_arg(argName, inputArgs, defaultVal)
iStart = 1;
wh = find(strcmp(inputArgs(iStart:end),argName));
if ~isempty(wh), output=inputArgs{wh+iStart}; else, output=defaultVal; end
end


