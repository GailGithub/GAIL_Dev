%% Generate Examples of Asian Arithmetic Mean Option Pricing
%function [muhat,aMLE,errVec,outAll] = TestAsianArithmeticMeanOptionAutoExample(dim,BernPolyOrder_1,...
%  ptransform,figSavePath,visiblePlot_1,arbMean_1,stopAtTol,samplingMethod,absTol_arg)
function [muhat,err,timeVec,outVec] = TestAsianArithmeticMeanOptionAutoExample(varargin)

% sampling = 'Lattice';
% figSavePath = 'D:/Mega/MyWriteupBackup/';
% newPath = strcat(figSavePath, sampling, '/', 'arbMean/');
% 
% inputArgs = {'dim',12, 'absTol',1e-2, 'order',2, ...
%                 'ptransform','Baker', 'stopAtTol',true, ...
%                 'figSavePath',newPath, 'arbMean',true, ...
%                 'samplingMethod',sampling, 'visiblePlot',true};
% varargin = inputArgs;

whichExample = 'Pierre';
dataFileName = [whichExample 'AsianCallExampleAllData.mat'];

if exist(dataFileName,'file')
  load(dataFileName)
  ArchEuroCall = AsianCall;
  ArchAsianCall = AsianCall;
  ArchAsianCallCV = AsianCall;
  Archnvec = nvec;
  ArchabsTolGold = absTolGold;
  ArchnGoldRep = nGoldRep;
  ArchAbsTol = absTol;
  ArchRelTol = relTol;
  ArchnRepAuto = nRepAuto;
  ArchnRep = nRep;
end

%% Parameters for the Asian option, Fred's Original
if strcmp(whichExample,'Fred')
  absTol = 1e-4;
  relTol = 0;
  inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for one quarter
  inp.assetParam.initPrice = 100; %initial stock price
  inp.assetParam.interest = 0.01; %risk-free interest rate
  inp.assetParam.volatility = 0.5; %volatility
  inp.payoffParam.strike = 100; %strike price
  inp.priceParam.absTol = absTol; %absolute tolerance of a penny
  inp.priceParam.relTol = relTol; %zero relative tolerance
  inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
  inp.bmParam.assembleType = 'PCA';
  inp.payoffParam.putCallType = {'call'};
  nvec = 2.^(7:17)';
  absTolGold = 1e-6;
  
  %% Parameters for the Asian option, Pierre's Example
elseif strcmp(whichExample,'Pierre')
  absTol = 1e-2; %1e-4;
  relTol = 0;
  inp.timeDim.timeVector = 1/12:1/12:1; %weekly monitoring for one quarter
  inp.assetParam.initPrice = 100; %initial stock price
  inp.assetParam.interest = 0.05; %risk-free interest rate
  inp.assetParam.volatility = 0.5; %volatility
  inp.payoffParam.strike = 100; %strike price
  inp.priceParam.absTol = absTol; %absolute tolerance of a penny
  inp.priceParam.relTol = relTol; %zero relative tolerance
  inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
  inp.bmParam.assembleType = 'PCA';
  inp.payoffParam.putCallType = {'call'};
  absTolGold = 1e-5;
  nvec = 2.^(7:17)';
end
nmax = max(nvec);
nRep = 100;
nlarge = nmax*2;
nn = numel(nvec);
alpha = 0.1;

%% Construct some different options
EuroCall = optPrice(inp); %construct a European optPrice object
AsianCall = optPrice(EuroCall); %construct an Asian optPrice object
AsianCall.payoffParam = struct( ...
  'optType',{{'amean'}},...
  'putCallType',{{'call'}});
AsianCallCV = optPrice(EuroCall); %construct an Asian and European optPayoff object for CV
AsianCallCV.payoffParam = struct( ...
  'optType',{{'amean','gmean'}},...
  'putCallType',{{'call','call'}});

%% Construct a very accurate answer
disp('Gold standard')
compGold = true;
nGoldRep = 100;
if exist('callPriceExact','var') && ...
    absTolGold == ArchabsTolGold && nGoldRep == ArchnGoldRep && ...
    all(ArchAsianCall.timeDim.timeVector == AsianCall.timeDim.timeVector) && ...
    ArchAsianCall.assetParam.initPrice == AsianCall.assetParam.initPrice && ... %initial stock price
    ArchAsianCall.payoffParam.strike == ArchAsianCall.payoffParam.strike, %strike price
  compGold = false;
  disp('Already have gold standard Asian Call')
end

if compGold
  fCV.func = @(x) genOptPayoffs(AsianCallCV,x);
  fCV.cv = AsianCallCV.exactPrice(2:end);
  d = AsianCallCV.timeDim.nSteps;
  callPriceGold(nGoldRep,1) = 0;
  tic
  for ii = 1:nGoldRep
    gail.TakeNote(ii,1) %print out every 10th ii
    callPriceGold(ii) = ...
      cubSobol_g(fCV,[zeros(1,d); ones(1,d)],'uniform',absTolGold,relTol);
    %       x = net(scramble(sobolset(AsianCall.timeDim.nSteps), ...
    %          'MatousekAffineOwen'),nGold);
    %       payoffAsianEuro = genOptPayoffs(AsianCall,x);
    %       callPriceGold(ii) = mean(payoffAsianEuro(:,1));
    %       payoffAsianEuro = genOptPayoffs(AsianCallCV,x);
    %       temp = payoffAsianEuro(:,1) - AsianCallCV.exactPrice(2) + payoffAsianEuro(:,2);
    %       callPriceGold(ii) = mean(temp);
  end
  callPriceExact = mean(callPriceGold);
  toc
end
disp(['mu  = ' num2str(callPriceExact,15) ' +/- ' num2str(2*std(callPriceGold),10)])
%return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try Bayesian automatic cubature
%disp('Automatic Bayesian cubature')
nRepAuto = 100;
timeVecAsianCallBayesAutoo(nRepAuto,1) = 0;
nVecAsianCallBayesAutoo(nRepAuto,1) = 0;
compCallAutoBayes = true;
if exist(dataFileName,'file')
  if exist('muAsianCallBayesAuto','var') && ...
      numel(nRepAuto) == numel(ArchnRepAuto) && ...
      ArchAbsTol == absTol && ...
      ArchRelTol == relTol
    compCallAutoBayes = false;
    disp('Already have automatic scrambled Bayesian Asian Call')
  end
end
if true
  integrand = @(x) genOptPayoffs(AsianCall,x);
  dim = AsianCall.timeDim.nSteps;
  muAsianCallBayesAuto(nRepAuto,1) = 0;
  %outCallBayes(nRepAuto,1) = 0;
  
  % input params initializations
  ptransform = get_arg('ptransform', varargin);
  stopAtTol = get_arg('stopAtTol', varargin);
  figSavePath = get_arg('figSavePath', varargin);
  visiblePlot = get_arg('visiblePlot', varargin);
  samplingMethod = get_arg('samplingMethod', varargin);

  fName = 'optPrice';
  
  inputArgs = varargin;
  inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
  inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;
  %inputArgs{end+1} = 'dim'; inputArgs{end+1} = dim;
  inputArgs = set_arg('dim', inputArgs, dim);

  % initialise the object based on the sampling method
  if exist('samplingMethod','var') && ...
      strcmp(samplingMethod,'Sobol') % use Sobol points
    obj=cubMLESobol(inputArgs{:});
  else % use Lattice points
    obj=cubBayesLattice_g(inputArgs{:});
  end

  
  tic
  for i =  1:nRepAuto
    gail.TakeNote(i,10)

    [muAsianCallBayesAuto(i),outCallBayes(i)] = compInteg(obj);
    timeVecAsianCallBayesAutoo(i) = outCallBayes(i).time;
    nVecAsianCallBayesAutoo(i) = outCallBayes(i).n;
  end
  toc
end
errvecAsianCallBayesAutoo = abs(callPriceExact - muAsianCallBayesAuto);
errmedAsianCallBayesAutoo = median(errvecAsianCallBayesAutoo);
errtopAsianCallBayesAutoo = quantile(errvecAsianCallBayesAutoo,1-obj.alpha);
rangeAsianCallBayesAutoo  = range(muAsianCallBayesAuto);
timetopAsianCallBayesAutoo = quantile(timeVecAsianCallBayesAutoo,1-obj.alpha);
ntopAsianCallBayesAutoo    = quantile(nVecAsianCallBayesAutoo,1-obj.alpha);
successAsianCallBayesAutoo = mean(errvecAsianCallBayesAutoo <= obj.absTol);

fprintf('\nError: Median %1.2e, Worst %1.2e, Range %1.2e, \n worstN %d, worstTime %1.3f, SuccessRatio %1.2f, \n absTol %1.2e, relTol %1.2e\n', ...
  errmedAsianCallBayesAutoo, errtopAsianCallBayesAutoo, rangeAsianCallBayesAutoo, ...
  ntopAsianCallBayesAutoo, timetopAsianCallBayesAutoo, ...
  successAsianCallBayesAutoo, outCallBayes(1).absTol, outCallBayes(1).relTol);


muhat = median(muAsianCallBayesAuto);
err = errvecAsianCallBayesAutoo;
outVec = outCallBayes;
timeVec = timeVecAsianCallBayesAutoo;

%% Save output
%% save(dataFileName)
end


% picks the input argument from the varargin cell array
function output = get_arg(argName, inputArgs, defaultVal)
iStart = 1;
wh = find(strcmp(inputArgs(iStart:end),argName));
if ~isempty(wh), output=inputArgs{wh+iStart}; else, output=defaultVal; end
end

function outArgs = set_arg(argName, inputArgs, val)
outArgs=inputArgs;
iStart=1;
wh=find(strcmp(outArgs(iStart:end),argName));
if ~isempty(wh), outArgs{wh+1}=val; else, outArgs{end+1}=argName; outArgs{end+2}=val; end
end

