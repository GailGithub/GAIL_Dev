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
ptransform = get_arg('ptransform', varargin);
stopAtTol = get_arg('stopAtTol', varargin);
figSavePath = get_arg('figSavePath', varargin);
visiblePlot = get_arg('visiblePlot', varargin);
samplingMethod = get_arg('samplingMethod', varargin);

% define the integrand 
fName='Keister';
integrand = @(x) keisterFunc(x,dim,0.8); % a=0.8
exactInteg = Keistertrue(dim);

% set the output dir to save the plots        
fullPath = strcat(figSavePath,'/',fName,'/',ptransform,'/');
if exist(fullPath,'dir')==false
  mkdir(fullPath);
end

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

nRep = 100;
muhatVec = zeros(nRep,1);

%plotObjectiveFunc(obj)
for i=1:nRep
  [muhatVec(i),outVec(i)]=compInteg(obj);
end
muhat = median(muhatVec);

%plotMLE_Loss(obj);
toc

errVec = abs(exactInteg - muhatVec);
timeVec = [outVec(:).time]';
out = outVec;


%% plot error
if stopAtTol==false
  BernPolyOrder = outVec{end}.order;
  nvec = 2.^outVec{end}.mvec;
  temp = [outVec(:).muhatAll];
  muhatVec = median(temp,2);
  temp = [outVec(:).ErrBdAll];
  ErrBdVec = median(temp,2);
  
  errCubatureVec = abs(exactInteg - muhatVec);
  plotCubatureError(dim, nvec, errCubatureVec, ErrBdVec, fName, BernPolyOrder, ...
    ptransform, fullPath,visiblePlot,arbMean, out.s_All, out.dscAll)

  % plot computation time vs number of points
  figSavePathName = sprintf('%s%s computeTime d_%d bernoulli_%d Period_%s.png', ...
    fullPath, fName, dim, BernPolyOrder, ptransform);
  plot_nvec_vs_computeTime(nvec, out.timeAll, visiblePlot, figSavePathName, samplingMethod)
end
fprintf('Done\n');

end

% picks the input argument from the varargin cell array
function output = get_arg(argName, inputArgs, defaultVal)
iStart = 1;
wh = find(strcmp(inputArgs(iStart:end),argName));
if ~isempty(wh), output=inputArgs{wh+iStart}; else, output=defaultVal; end
end


