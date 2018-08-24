%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isunix
  !synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
end
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

if isunix
  figSavePath = '/home/jagadees/MyWriteup/';
else
  figSavePath = 'D:/Mega/MyWriteupBackup/';
end
figSavePath = strcat(figSavePath, 'Jul_2ndweek2018/');

if exist(figSavePath,'dir')==false
  mkdir(figSavePath);
end

% log the results
completereport = strcat(figSavePath,...
  '_tests-logs-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

visiblePlot=true;

%
% https://www.mathworks.com/matlabcentral/answers/
% 98969-how-can-i-temporarily-avoid-figures-to-be-displayed-in-matlab
%
if visiblePlot==false
  set(0,'DefaultFigureVisible','off')
else
  set(0,'DefaultFigureVisible','on')
end

rng(202326) % control random number generation

stopAtTol = true;
alpha = 0.01;

dim = 2;
errTol = 0.001;
relTol = 0;
stopAtTol = true;

log10ErrVec = -2:-1:-5;  % 1E-6 cannot be computed
errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
errTolVec = 10.^log10ErrVec;

% d-3 problem reduced to d-2 using Genz method
MVNParams.C = [4 1 1; 0 1 0.5; 0 0 0.25];
MVNParams.Cov = MVNParams.C'*MVNParams.C;
MVNParams.a = [-6 -2 -2];
MVNParams.b = [5 2 1];
MVNParams.mu = 0;
MVNParams.CovProp.C = chol(MVNParams.Cov)';
muBest = '0.676337324357787483819492990733124315738677978515625';
muBest = str2double(muBest);

integrand = @(t) GenzFunc(t,MVNParams);

nRep = 100;
muhatVec(nRep,1) = 0;
indx = 1;

for errTol=errTolVec(1:end)
  for i=1:nRep
    fprintf('%1.1e \n', errTol)
    inputArgs = {'fName','MVN','dim',dim, 'absTol',errTol,'relTol',relTol, ....
      'stopAtTol',stopAtTol,'f',integrand };
    tic
    [muhatVec(indx),outVec(indx)] = cubBayesMLE_Matern_g(inputArgs{:});
    toc
    indx = indx+1;
  end
end

save('matern_guranteed.mat', 'muhatVec', 'outVec','inputArgs',...
  'MVNParams','muBest')

fprintf('done')

