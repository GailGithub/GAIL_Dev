%% Keister's Example of Multidimensional Integration using Bayesian cubature
%
% B. D. Keister, Multidimensional quadrature algorithms, _Computers in
% Physics_, *10*, pp. 119-122, 1996, presents the following
% multidimensional integral, inspired by a physics application:
%
% \[ I = \int_{\mathbb{R}^d} \cos(\lVert \boldsymbol{x} \rVert)
% \exp(-\lVert \boldsymbol{x} \rVert^2) \, \mathrm{d} \boldsymbol{x},
% \qquad d = 1, 2, \ldots. \]

function [muhat,errVec,timeVec,outVec] = TestKeisterBayesianCubature(varargin)

% input params initializations
dim = get_arg('dim', varargin);
samplingMethod = get_arg('samplingMethod', varargin);

% define the integrand 
fName='Keister';
integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
exactInteg = Keistertrue(dim);

tic

inputArgs = varargin;
inputArgs{end+1} = 'f'; inputArgs{end+1} = integrand;
inputArgs{end+1} = 'fName'; inputArgs{end+1} = fName;

% initialise the object based on the sampling method
if exist('samplingMethod','var') && ...
    strcmp(samplingMethod,'Sobol') % use Sobol points
  obj=cubMLESobol(inputArgs{:});
else % use Lattice points
  obj=cubBayesLattice_g(inputArgs{:});
end
%plotObjectiveFunc(obj);

nRep = get_arg('nRepAuto', varargin, 100);
muhatVec = zeros(nRep,1);

%plotObjectiveFunc(obj)
for i=1:nRep
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


