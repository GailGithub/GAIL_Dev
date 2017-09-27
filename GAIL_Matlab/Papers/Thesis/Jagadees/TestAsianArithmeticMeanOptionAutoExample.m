%% Generate Examples of Asian Arithmetic Mean Option Pricing
function [muhat,aMLE,errVec,outAll] = TestAsianArithmeticMeanOptionAutoExample(dim,BernPolyOrder_1,...
  ptransform_1,figSavePath_1,visiblePlot_1,arbMean_1,testAll_arg,absTol_arg)

whichExample = 'Pierre'
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
disp('Automatic Bayesian cubature')
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
  f = @(x) genOptPayoffs(AsianCall,x);
  dim = AsianCall.timeDim.nSteps;
  muAsianCallBayesAuto(nRepAuto,1) = 0;
  outCallBayes = cell(nRepAuto,1);
  
  if ~exist('figSavePath_1','var')
    figSavePath_1 = '/home/jagadees/MyWriteup/Sep2ndweek_optprice/';
  end
  if ~exist('ptransform_1','var')
    ptransform_1 = 'Baker';
  end
  if ~exist('testAll_1','var')
    testAll_1 = false;
  end
  if ~exist('fName_1','var')
    fName_1 = 'optPrice';
  end
  if ~exist('BernPolyOrder_1','var')
    BernPolyOrder_1 = 4;
  end
  if ~exist('arbMean_1','var')
    arbMean_1 = false;
  end
  if ~exist('absTol_arg','var')
    absTol_arg = absTol;
  end
  
  tic
  for i =  1:nRepAuto
    gail.TakeNote(i,10)
    [muAsianCallBayesAuto(i),outCallBayes{i}] = ...
      cubMLELattice(f,dim,absTol_arg,relTol,...
            BernPolyOrder_1,ptransform_1,testAll_1,figSavePath_1,...
            fName_1,arbMean_1);
    timeVecAsianCallBayesAutoo(i) = outCallBayes{i}.time;
    nVecAsianCallBayesAutoo(i) = outCallBayes{i}.n;
  end
  toc
end
errvecAsianCallBayesAutoo = abs(callPriceExact - muAsianCallBayesAuto);
errmedAsianCallBayesAutoo = median(errvecAsianCallBayesAutoo);
errtopAsianCallBayesAutoo = quantile(errvecAsianCallBayesAutoo,1-alpha);
rangeAsianCallBayesAutoo  = range(muAsianCallBayesAuto);
timetopAsianCallBayesAutoo = quantile(timeVecAsianCallBayesAutoo,1-alpha);
ntopAsianCallBayesAutoo    = quantile(nVecAsianCallBayesAutoo,1-alpha);
successAsianCallBayesAutoo = mean(errvecAsianCallBayesAutoo <= absTol);

muhat = median(muAsianCallBayesAuto);
aMLE = 0;
errVec = errvecAsianCallBayesAutoo;
outAll = outCallBayes;

%% Save output
%% save(dataFileName)


