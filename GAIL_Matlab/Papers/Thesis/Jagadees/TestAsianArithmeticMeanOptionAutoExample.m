%% Generate Examples of Asian Arithmetic Mean Option Pricing
function [muhat,err,timeVec,outVec] = TestAsianArithmeticMeanOptionAutoExample(varargin)

whichExample = 'Pierre';
dataFileName = [whichExample 'AsianCallExampleAllData.mat'];

if exist(dataFileName,'file')
  S = load(dataFileName);
  callPriceExact = S.callPriceExact;
end

%% Parameters for the Asian option, Pierre's Example
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

%% Construct some different options
EuroCall = optPrice(inp); %construct a European optPrice object
AsianCall = optPrice(EuroCall); %construct an Asian optPrice object
AsianCall.payoffParam = struct( ...
  'optType',{{'amean'}},...
  'putCallType',{{'call'}});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bayesian automatic cubature
% read input params
samplingMethod = get_arg('samplingMethod', varargin, 'Lattice');
nRepAuto = get_arg('nRepAuto', varargin, 100);

timeVecAsianCallBayesAutoo(nRepAuto,1) = 0;
nVecAsianCallBayesAutoo(nRepAuto,1) = 0;
integrand = @(x) genOptPayoffs_fixNan(AsianCall,x);

dim = AsianCall.timeDim.nSteps;
muAsianCallBayesAuto(nRepAuto,1) = 0;
%outCallBayes(nRepAuto,1) = 0;

fName = 'optPrice';
inputArgs = varargin;
inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;
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
  tic
  [muAsianCallBayesAuto(i),outCallBayes(i)] = compInteg(obj);
  toc
  timeVecAsianCallBayesAutoo(i) = outCallBayes(i).time;
  nVecAsianCallBayesAutoo(i) = outCallBayes(i).n;
end
toc

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

end

% reset NanN vlaues to zero
function y = genOptPayoffs_fixNan(AsianCall,x)
y = genOptPayoffs(AsianCall,x);
y(isnan(y)) = 0;
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
if ~isempty(wh), outArgs{wh+1}=val; else, outArgs{end+1}=argName; outArgs{end+1}=val; end
end

