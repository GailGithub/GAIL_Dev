%% Interface to test Multivariate Normal Probabilities
%
function [muhat,errVec,timeVec,outVec] = TestMVN_BayesianCubature(varargin)

dim = get_arg('dim', varargin);
samplingMethod = get_arg('samplingMethod', varargin);

if dim==2
  % d=3 problem reduced to d=2 using Genz method
  MVNParams.C = [4 1 1; 0 1 0.5; 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2];
  MVNParams.b = [5 2 1];
  MVNParams.mu = 0;
elseif dim==3
  % d=4 problem reduced to d=3 using Genz method
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

muBest = computeGoldenMuhat(MVNParams);

%% Bayseian cubature with Shif invariant kernel and rank-1 Lattice points
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

tic

nRep = 100; % increased it for gail plots

muhatVec(nRep,1) = 0;
for i = 1:nRep
  [muhatVec(i),outVec(i)]=compInteg(objCubBayes);
end
timeVec = [outVec(:).time]';
errVec = abs(muBest - muhatVec);

toc

muhat = median(muhatVec);
fprintf('done\n');

end

function muBest = computeGoldenMuhat(MVNParams)

dim = length(MVNParams.a)-1;
if 0 % if we need to recompute

if dim==2
  dataFileName = 'MVNProbExampleDim3.mat';
elsif dim==3
  dataFileName = 'MVNProbExampleDim4.mat';
end

if exist(dataFileName,'file')
  S = load(dataFileName);
  MVNProbBestArch = S.MVNProbBest;
  nRepGoldArch = S.nRepGold;
end

%% First compute a high accuracy answer
nRep = get_arg('nRepAuto', varargin, 100);

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
  
  save(dataFileName,'muBest', 'muBestvec', 'nRepGold', 'MVNProbBest');
end
%disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])
end

% use precomputed values
if dim==2
  muBest = '0.67633732435778748381';
elseif dim==3
  muBest = '0.674516483073121952962';
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


