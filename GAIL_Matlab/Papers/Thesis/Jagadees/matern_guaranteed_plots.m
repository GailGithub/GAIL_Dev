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

alpha = 0.01;
errTol = 0.001;
relTol = 0;
stopAtTol = true;

log10ErrVec = -2:-1:-5;  % 1E-6 cannot be computed
errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
errTolVec = 10.^log10ErrVec;
orig_dim = 3;  % MVN problem dimension

if orig_dim==3
  % d=3 problem reduced to d=2 using Genz method
  dim = 2;  % Genz reduced dim
  MVNParams.C = [4 1 1; 0 1 0.5; 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2];
  MVNParams.b = [5 2 1];
  MVNParams.mu = 0;
  MVNParams.CovProp.C = chol(MVNParams.Cov)';
  muBest = '0.676337324357787483819492990733124315738677978515625';
  muBest = str2double(muBest);
elseif orig_dim == 4
  dim = 3;
  MVNParams.C = [4 1 1 1; 0 1 0.5 0.5; 0 0 0.25 0.25; 0 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2 -2];
  MVNParams.b = [5 2 1 2];
  MVNParams.mu = 0;
  MVNParams.CovProp.C = chol(MVNParams.Cov)';
  muBest = '0.67451648307312195296248091835877858102321624755859375';
  muBest = str2double(muBest);
end
  
integrand = @(t) GenzFunc(t,MVNParams);
fName = 'MVN';
nRep = 100;
muhatVec(nRep,1) = 0;
indx = 1;

matFilePath = 'D:\Dropbox\fjhickernellGithub\GAIL_Dev-BayesianCubature\GAIL_Matlab\Papers\Thesis\Jagadees\Paper2018\figures\';
filename = 'matern_guranteed_2018-Aug-21';
matFileName = [matFilePath filename '.mat'];
if exist(matFileName, 'file') == 2
  S = load(matFileName);
  
else

  for errTol=errTolVec(end:-1:1)
    tic
    for i=1:nRep
      fprintf('%d  %1.1e', indx, errTol)
      inputArgs = {'fName',fName,'dim',dim, 'absTol',errTol,'relTol',relTol, ....
        'stopAtTol',stopAtTol,'f',integrand };
      [muhatVec(indx),outVec(indx)] = cubBayesMLE_Matern_g(inputArgs{:});
      fprintf(', time %1.3f n %d \n', outVec(indx).time, outVec(indx).n)
      indx = indx+1;
    end
    toc
  end

  timeStamp = datetime('now','Format','y-MMM-d');

  save(sprintf('matern_guranteed_%s.mat', timeStamp), 'muhatVec', 'outVec',...
    'inputArgs','MVNParams','muBest','timeStamp')
  S = load(sprintf('matern_guranteed_%s.mat', timeStamp));
end

fprintf('done')

timeVec = [S.outVec.time];
errVec = [abs(S.muhatVec-S.muBest)];
tolVec = [S.outVec.absTol]';
nVec = [S.outVec.n]';
errVec = errVec./tolVec;
S.timeVec = reshape(timeVec, [], length(errTolVec));
S.errVec = reshape(errVec, [], length(errTolVec));
S.nVec = reshape(nVec, [], length(errTolVec));
S.log10ErrVec = log10ErrVec;
% S.tolVec = repmat(errTolVec, length(timeVec)/length(errTolVec),1);
S.tolVec = reshape(tolVec, [], length(errTolVec));
S.testFunArg = struct(S.inputArgs{:});
timeStamp = S.timeStamp;

pointSize = 30;
pointShapes = {'o','s','d','^','v','<','>','p','h'};
errVecLimits = [1E-6, 1E1];

figH = figure();
set(figH, 'units', 'inches', 'Position', [1 1 9 6])
%timeLimits = [1E-3, 1E-2];
%timeAxisLimits(1) = min(S.timeVec(:))/2;
%timeAxisLimits(2) = max(S.timeVec(:))*2;
%plot([1, 1], timeAxisLimits, 'r', 'LineWidth',1)
timeTicksLimits(1) = floor(log10(min(S.timeVec(:)))); 
timeTicksLimits(2) = ceil(log10(max(S.timeVec(:))));
timeTickInterval = ceil(diff(timeTicksLimits)/5);
plot([1, 1], 10.^timeTicksLimits, 'r', 'LineWidth',1)
hold on

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.timeVec(:,i),pointSize,log10(S.tolVec(:,i)),...
    pointShapes{i},'filled')
end

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\({\vert\mu-\widehat{\mu} \vert}/{\varepsilon}\)')
ylabel('Time (secs)')
[V,I] = sort(S.log10ErrVec);
c = colorbar('Direction','reverse', 'Ticks',V, ...
  'TickLabels',errTolVecText(I), 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Error Tolerance, $\varepsilon$';
% axis tight; not required

%assert(max(timeVec(:)) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))
axis([errVecLimits(1) errVecLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2) ])
timeTicksVec = (timeTicksLimits(1) :timeTickInterval: timeTicksLimits(2));
if timeTicksVec(end) ~= timeTicksLimits(2)
  timeTicksVec(end+1) = timeTicksLimits(2);
end
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^timeTicksVec))

figSavePathName = sprintf('%s%s_guaranteed_time_Matern_d%d_%s.png', ...
  figSavePath, S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
saveas(figH, figSavePathName)


figH1 = figure();
set(figH1, 'units', 'inches', 'Position', [1 1 9 6])

scatter(S.errVec(:),S.timeVec(:),pointSize,'o','filled')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Error Tolerance, $\varepsilon$')
ylabel('Time (secs)')
errTolLimits = [1E-6 1E-1];
axis([errTolLimits(1) errTolLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2)])

set(gca,'Xtick',(10.^(log10(errTolLimits(1)):2:log10(errTolLimits(2)))), ...
   'YTick',(10.^(timeTicksVec)) )
set(gca, 'XDir','reverse')

figSavePathName1 = sprintf('%s%s_rapid_time_Matern_d%d_%s.png', ...
  figSavePath, S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
saveas(figH1, figSavePathName1)
 

%
% Num. points vs time growth
figH2 = figure();
set(figH2, 'units', 'inches', 'Position', [1 1 9 6])

scatter(S.nVec(:),S.timeVec(:),pointSize,'o','filled')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Num. Samples, $n$')
ylabel('Time (secs)')
nvecLimits = [1E1 1E4];
axis([nvecLimits(1) nvecLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2)])

set(gca,'Xtick',(10.^(log10(nvecLimits(1)):1:log10(nvecLimits(2)))), ...
  'YTick',(10.^(timeTicksVec)) )

% hold on
% loglog([min(S.nVec(:)) max(S.nVec(:))],...
%   min(S.timeVec(:))*[1 exp(nvecLimits(end)/nvecLimits(1)) ], 'g--')

figSavePathName2 = sprintf('%s%s_rapid_n_vs_time_Matern_d%d_%s.png', ...
  figSavePath, S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
saveas(figH2, figSavePathName2)

fprintf('done')



