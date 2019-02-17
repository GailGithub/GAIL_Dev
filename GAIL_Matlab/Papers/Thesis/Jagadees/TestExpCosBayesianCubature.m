%% Integrate Exp(cos) using Bayesian cubature
%
%
function [muhat,err,time,out] = TestExpCosBayesianCubature(varargin)

% input params initializations
dim = get_arg('dim', varargin);
ptransform = get_arg('ptransform', varargin);
stopAtTol = get_arg('stopAtTol', varargin);
figSavePath = get_arg('figSavePath', varargin);
visiblePlot = get_arg('visiblePlot', varargin);
samplingMethod = get_arg('samplingMethod', varargin);

% define the integrand 
integrand = @(x) exp(sum(cos(2*pi*x), 2));
exactInteg = besseli(0,1)^dim;
fName = 'Exp(cos)';

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
  obj=cubBayesNet_g(inputArgs{:});
else % use Lattice points
  obj=cubBayesLattice_g(inputArgs{:});
end

nRep = 100;
muhatVec = zeros(nRep,1);

for i=1:nRep
  [muhatVec(i),outVec(i)]=compInteg(obj);
end
muhat = median(muhatVec);

%plotMLE_Loss(obj);
toc

err = abs(exactInteg - muhat);
time = median([outVec(:).time]);
out = outVec;

if 0
  % Example to plot cubature error
  trueError = abs(outVec(1).optParams.muhatAll - exactInteg);
  estError = outVec(1).optParams.ErrBdAll;
  visiblePlot = true;
  scale = outVec(1).optParams.s_All;
  dsc = outVec(1).optParams.dscAll;
  plotCubatureError(obj.dim, 2.^obj.mvec, trueError, estError, obj.fName, obj.order, ...
    obj.ptransform, obj.figSavePath, visiblePlot, obj.arbMean, scale, dsc)
end

if err/obj.absTol > 1
  error('wait')
end
%% plot error
if stopAtTol==false
  
  nvec = 2.^outVec(1).mvec;

  temp = [outVec(:).muhatAll];
  muhatVec = median(temp,2);
  temp = [outVec(:).ErrBdAll];
  ErrBdVec = median(temp,2);
  
  errCubatureVec = abs(exactInteg - muhatVec);
  plotCubatureError(dim, nvec, errCubatureVec, ErrBdVec, fName, BernPolyOrder, ptransform, ...
    fullPath,visiblePlot,arbMean, out.s_All, out.dscAll)

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

