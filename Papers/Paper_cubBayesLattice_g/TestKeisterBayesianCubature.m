%% Keister's Example of Multidimensional Integration using Bayesian cubature
%
% B. D. Keister, Multidimensional quadrature algorithms, _Computers in
% Physics_, *10*, pp. 119-122, 1996, presents the following
% multidimensional integral, inspired by a physics application:
%
% \[ I = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, \mathrm{d} \boldsymbol{x},
% \qquad d = 1, 2, \ldots. \]

function [muhat,errVec,timeVec,outVec,errTolVec] = TestKeisterBayesianCubature(varargin)

% input params initializations
dim = get_arg('dim', varargin);
samplingMethod = get_arg('samplingMethod', varargin);
log10ErrVec = get_arg('log10ErrVec', varargin);

% define the integrand 
fName='Keister';
integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
exactInteg = Keistertrue(dim);

inputArgs = varargin;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;

nRep = get_arg('nRepAuto', varargin, 100);
muhatVec = zeros(nRep,1);
errTolVec(nRep,1) = 0;

log10ErrTol_a = log10ErrVec(1) - 0.3;
log10ErrTol_b = log10ErrVec(end) + 0.3;
randErrTol = @() 10^(log10ErrTol_a + (log10ErrTol_b - log10ErrTol_a)*rand());

tic


for i=1:nRep
  errTolVec(i) = randErrTol();
  inputArgs = set_arg('absTol', inputArgs, errTolVec(i));
  
  % initialise the object based on the sampling method
  if exist('samplingMethod','var') && ...
      strcmp(samplingMethod,'Net') % use Sobol points
    obj=cubBayesNet_g(integrand,dim,inputArgs{:});
  else % use Lattice points
    obj=cubBayesLattice_g(integrand,dim,inputArgs{:});
  end
    [muhatVec(i),outVec(i)]=compInteg(obj);
  end
muhat = median(muhatVec);

toc

errVec = abs(exactInteg - muhatVec);
timeVec = [outVec(:).time]';

end

% picks the input argument from the varargin cell array
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
