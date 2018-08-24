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
fName = 'MVN';
nRep = 100;
muhatVec(nRep,1) = 0;
indx = 1;

for errTol=errTolVec(1:end)
  for i=1:nRep
    fprintf('%1.1e \n', errTol)
    inputArgs = {'fName',fName,'dim',dim, 'absTol',errTol,'relTol',relTol, ....
      'stopAtTol',stopAtTol,'f',integrand };
    tic
    [muhatVec(indx),outVec(indx)] = cubBayesMLE_Matern_g(inputArgs{:});
    toc
    indx = indx+1;
  end
end

timeStamp = datetime('now','Format','y-MMM-d');

save('matern_guranteed.mat', 'muhatVec', 'outVec','inputArgs',...
  'MVNParams','muBest')

fprintf('done')

timeVec = [outVec.time];
errVec = abs(muhatVec-muBest);
S.timeVec = reshape(timeVec, [], length(errTolVec));
S.errVec = reshape(errVec, [], length(errTolVec));
S.log10ErrVec = log10ErrVec;
S.tolVec = repmat(errTolVec, 100,1);
S.testFunArg = struct(inputArgs{:});
S.timeStamp = timeStamp;
pointSize = 30;
pointShapes = {'o','s','d','^','v','<','>','p','h'};
errVecLimits = [1E-7, 1E1];

figH = figure();
set(figH, 'units', 'inches', 'Position', [1 1 9 6])
%timeLimits = [1E-3, 1E-2];
timeAxisLimits(1) = min(S.timeVec(:))/2;
timeAxisLimits(2) = max(S.timeVec(:))*2;
plot([1, 1], timeAxisLimits, 'r', 'LineWidth',1)
hold on

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.timeVec(:,i),pointSize,log10(S.tolVec(:,i)),...
    pointShapes{i},'filled')
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
ylabel('Time (secs)')
[V,I] = sort(S.log10ErrVec);
c = colorbar('Direction','reverse', 'Ticks',V, ...
  'TickLabels',errTolVecText(I), 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';
% axis tight; not required

%assert(max(timeVec(:)) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))
axis([errVecLimits(1) errVecLimits(2) timeAxisLimits(1) timeAxisLimits(2) ])

timeTicksLimits(1) = floor(log10(min(S.timeVec(:))));
timeTicksLimits(2) = floor(log10(max(S.timeVec(:))));
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^(timeTicksLimits(1) :1: timeTicksLimits(2))))

figSavePathName = sprintf('%s%s_guaranteed_time_Matern_d%d_%s.png', ...
  figSavePath, S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
saveas(figH, figSavePathName)


figH1 = figure();
set(figH1, 'units', 'inches', 'Position', [1 1 9 6])

scatter(S.tolVec(:),S.timeVec(:),pointSize,'o','filled')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Error Tolerance, $\varepsilon$')
ylabel('Time (secs)')
errTolLimits = [1E-6 1E-1];
axis([errTolLimits(1) errTolLimits(2) timeAxisLimits(1) timeAxisLimits(2)])

set(gca,'Xtick',(10.^(log10(errTolLimits(1)):2:log10(errTolLimits(2)))), ...
  'YTick',(10.^(timeTicksLimits(1) :1: timeTicksLimits(2))))

figSavePathName = sprintf('%s%s_rapid_time_Matern_d%d_%s.png', ...
  figSavePath, S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
saveas(figH1, figSavePathName)


