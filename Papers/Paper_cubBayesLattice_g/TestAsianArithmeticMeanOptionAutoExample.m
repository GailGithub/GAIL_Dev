%% Interface to test Asian Arithmetic Mean Option Pricing
% using Bayesian cubature algorithm
function [muhat,err,timeVec,outVec,errTolVec] = TestAsianArithmeticMeanOptionAutoExample(varargin)

whichExample = 'Pierre';
dataFileName = [whichExample 'AsianCallExampleAllData.mat'];

if exist(dataFileName,'file')
  S = load(dataFileName);
  callPriceExact = S.callPriceExact;
else
  callPriceExact = str2double('13.121994282389803033');
end

%% Parameters for the Asian option, Pierre's Example
absTol = 1e-2;
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
log10ErrVec = get_arg('log10ErrVec', varargin);

integrand = @(x) genOptPayoffs_fixNan(AsianCall,x);

dim = AsianCall.timeDim.nSteps;
muAsianCallBayesAuto(nRepAuto,1) = 0;
errTolVec(nRepAuto,1) = 0;

log10ErrTol_a = log10ErrVec(1) - 0.3;
log10ErrTol_b = log10ErrVec(end) + 0.3;
randErrTol = @() 10^(log10ErrTol_a + (log10ErrTol_b - log10ErrTol_a)*rand());

fName = 'optPrice';
inputArgs = varargin;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;

tStart=tic;
for i =  1:nRepAuto
  errTolVec(i) = randErrTol();
  inputArgs = set_arg('absTol', inputArgs, errTolVec(i));
  
  % initialise the object based on the sampling method
  if exist('samplingMethod','var') && ...
      strcmp(samplingMethod,'Net') % use Sobol points
    obj=cubBayesNet_g(integrand,dim,inputArgs{:});
  else % use Lattice points
    obj=cubBayesLattice_g(integrand,dim,inputArgs{:});
  end

  [muAsianCallBayesAuto(i),outCallBayes(i)] = compInteg(obj);
end
disp(' ');
toc(tStart)

timeVecAsianCallBayesAuto = [outCallBayes(:).time]';
nVecAsianCallBayesAuto = [outCallBayes(:).n]';

errvecAsianCallBayesAuto = abs(callPriceExact - muAsianCallBayesAuto);
errmedAsianCallBayesAuto = median(errvecAsianCallBayesAuto);
errtopAsianCallBayesAuto = quantile(errvecAsianCallBayesAuto,1-obj.alpha);
rangeAsianCallBayesAuto  = range(muAsianCallBayesAuto);
timetopAsianCallBayesAuto = quantile(timeVecAsianCallBayesAuto,1-obj.alpha);
ntopAsianCallBayesAuto    = quantile(nVecAsianCallBayesAuto,1-obj.alpha);
successAsianCallBayesAuto = mean(errvecAsianCallBayesAuto <= obj.absTol);

fprintf('\n err: Median %1.2e, Worst %1.2e, Range %1.2e, \n worstN %d, worstTime %1.3f, SuccessRatio %1.2f, \n absTol %1.2e, relTol %1.2e\n', ...
  errmedAsianCallBayesAuto, errtopAsianCallBayesAuto, rangeAsianCallBayesAuto, ...
  ntopAsianCallBayesAuto, timetopAsianCallBayesAuto, ...
  successAsianCallBayesAuto, outCallBayes(1).optParams.absTol, ...
  outCallBayes(1).optParams.relTol);

muhat = median(muAsianCallBayesAuto);
err = errvecAsianCallBayesAuto;
outVec = outCallBayes;
timeVec = timeVecAsianCallBayesAuto;

end

% reset NaN values to zero
function y = genOptPayoffs_fixNan(AsianCall,x)
y = genOptPayoffs(AsianCall,x);
y(isnan(y)) = 0;
end

% picks the value of given argument from varargin cell array
function output = get_arg(argName, inputArgs, defaultVal)
iStart = 1;
wh = find(strcmp(inputArgs(iStart:end),argName));
if ~isempty(wh), output=inputArgs{wh+iStart}; else, output=defaultVal; end
end

% sets the argument of the given argument
function outArgs = set_arg(argName, inputArgs, val)
outArgs=inputArgs;
iStart=1;
wh=find(strcmp(outArgs(iStart:end),argName));
if ~isempty(wh), outArgs{wh+1}=val; else, outArgs{end+1}=argName; outArgs{end+1}=val; end
end

