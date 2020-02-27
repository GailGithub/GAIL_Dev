%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Creates the plots for the Matern kernel example
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

% path where the GAIL is installed
GAIL_path = GAILstart(0);
logSavePath=[GAIL_path,'OutputFiles',filesep, 'Paper_cubBayesLattice_g', filesep];

if exist(logSavePath,'dir')==false
  mkdir(logSavePath);
end

% log the results
completereport = strcat(logSavePath, filesep, ...
  '_tests-logs-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

visiblePlot=true;
% temporarily avoid figures to be displayed in matlab
if visiblePlot==false
  set(0,'DefaultFigureVisible','off')
else
  set(0,'DefaultFigureVisible','on')
end

rng(202326) % initialize random number generation to enable reproducability

alpha = 0.01;
errTol = 0.001;
relTol = 0;
stopAtTol = true;

log10ErrVec = -2:-1:-5;  % 1E-6 cannot be computed
errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
log10ErrTol_a = min(log10ErrVec) - 0.3;
log10ErrTol_b = max(log10ErrVec) + 0.3;
randErrTol = @() 10^(log10ErrTol_a + (log10ErrTol_b - log10ErrTol_a)*rand());
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
  muBest = '0.676337324357787483819';
  muBest = str2double(muBest);
elseif orig_dim == 4
  dim = 3;
  MVNParams.C = [4 1 1 1; 0 1 0.5 0.5; 0 0 0.25 0.25; 0 0 0 0.25];
  MVNParams.Cov = MVNParams.C'*MVNParams.C;
  MVNParams.a = [-6 -2 -2 -2];
  MVNParams.b = [5 2 1 2];
  MVNParams.mu = 0;
  MVNParams.CovProp.C = chol(MVNParams.Cov)';
  muBest = '0.674516483073121952962';
  muBest = str2double(muBest);
end
  
integrand = @(t) GenzFunc(t,MVNParams);
fName = 'MVN';
nRep = 100;
muhatVec(nRep,length(log10ErrVec)) = 0;
indx = 1;


% path to save/read the .mat files
figSavePath = logSavePath;
filename = 'matern_guranteed_2019-Jun-29';
matFileName = [figSavePath filesep filename '.mat'];
if exist(matFileName, 'file') == 2
  % data already exists, just load it
  S = load(matFileName);
else
  % recompute the data required to plot
  for k=1:length(log10ErrVec)
    tic
    for i=1:nRep
      errTol = randErrTol();
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
  matfilename = [figSavePath, filesep, ...
    sprintf('matern_guranteed_%s', timeStamp), '.mat'];
  save(matfilename, ...
    'muhatVec', 'outVec','inputArgs','MVNParams','muBest','timeStamp')
  S = load(matfilename);
end

fprintf('done')

timeVec = [S.outVec.time];
errVec = [abs(S.muhatVec-S.muBest)];
tolVec = reshape([S.outVec.absTol], size(errVec));
nVec = [S.outVec.n];
errVec = errVec./tolVec;
S.timeVec = reshape(timeVec, [], length(log10ErrVec));
S.errVec = reshape(errVec, [], length(log10ErrVec));
S.nVec = reshape(nVec, [], length(log10ErrVec));
S.log10ErrVec = log10ErrVec;
% S.tolVec = repmat(errTolVec, length(timeVec)/length(errTolVec),1);
S.tolVec = reshape(tolVec, [], length(log10ErrVec));
S.testFunArg = struct(S.inputArgs{:});
timeStamp = S.timeStamp;

pointSize = 30;
pointShapes = {'o','s','d','^','v','<','>','p','h'};
errVecLimits = [1E-6, 1E1];

% guaranteed plot with colorbar showing \varepsilon, the error tolerance
figH = figure();
set(figH, 'units', 'inches', 'Position', [1 1 9 5])
timeTicksLimits(1) = floor(log10(min(S.timeVec(:)))); 
timeTicksLimits(2) = ceil(log10(max(S.timeVec(:))));
timeTickInterval = ceil(diff(timeTicksLimits)/5);
plot([1, 1], 10.^timeTicksLimits, 'r', 'LineWidth',1)
hold on


% pointshapes represent number of samples
i=1;
for nval=unique(S.nVec)'
  index = find(S.nVec==nval);
  scatter(S.errVec(index),S.timeVec(index),pointSize,log10(S.tolVec(index)),...
    pointShapes{i},'filled')
  i = i+1;
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

temp = arrayfun(@(x){sprintf('${n={%d}}$', x)}, unique(S.nVec));
temp = [{'${\vert\mu-\widehat{\mu} \vert}/{\varepsilon}=1$'}, temp'];

axis([errVecLimits(1) errVecLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2) ])
timeTicksVec = (timeTicksLimits(1) :timeTickInterval: timeTicksLimits(2));
if timeTicksVec(end) ~= timeTicksLimits(2)
  timeTicksVec(end+1) = timeTicksLimits(2);
end
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^timeTicksVec))

figSavePathName = sprintf('%s_guaranteed_time_Matern_d%d_%s', ...
  S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
gail.save_image(figH, 'Paper_cubBayesLattice_g', figSavePathName)


% colorbar shows 'n', the number of samples
figH1 = figure();
set(figH1, 'units', 'inches', 'Position', [1 1 9 5])
timeTicksLimits(1) = floor(log10(min(S.timeVec(:)))); 
timeTicksLimits(2) = ceil(log10(max(S.timeVec(:))));
timeTickInterval = ceil(diff(timeTicksLimits)/5);
plot([1, 1], 10.^timeTicksLimits, 'r', 'LineWidth',1)
hold on

for i=1:size(S.errVec,2)
  scatter(S.errVec(:,i),S.timeVec(:,i),pointSize,log2(S.nVec(:,i)),...
    pointShapes{i},'filled')
end

temp = arrayfun(@(x){sprintf('${\\varepsilon=10^{%d}}$', x)}, log10ErrVec(end:-1:1));
temp = [{'${\vert\mu-\widehat{\mu} \vert}/{\varepsilon}=1$'}, temp];
legend(temp,'location','best','Interpreter','latex'); axis tight

set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('\({\vert\mu-\widehat{\mu} \vert}/{\varepsilon}\)')
ylabel('Time (secs)')
% cticks = log2(sort(unique(S.nVec(:))));
cticks = [4 6 8 10 12];
cticks_text = arrayfun(@(x){sprintf('2^{%d}', x)}, cticks);
c = colorbar('Ticks',cticks, ...
  'TickLabels',cticks_text, 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'Sample Size, $n$';
% axis tight; not required

axis([errVecLimits(1) errVecLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2) ])
timeTicksVec = (timeTicksLimits(1) :timeTickInterval: timeTicksLimits(2));
if timeTicksVec(end) ~= timeTicksLimits(2)
  timeTicksVec(end+1) = timeTicksLimits(2);
end
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):3:log10(errVecLimits(2)))), ...
  'YTick',(10.^timeTicksVec))

figSavePathName = sprintf('%s_guaranteed_time_with_n_Matern_d%d_%s', ...
  S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
gail.save_image(figH1, 'Paper_cubBayesLattice_g', figSavePathName)

%
% Num. points vs time growth
figH2 = figure();
set(figH2, 'units', 'inches', 'Position', [1 1 9 5])

scatter(S.nVec(:),S.timeVec(:),pointSize,'o','filled')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Number of Samples, $n$')
ylabel('Time (secs)')
nvecLimits = [1E1 1E4];
axis([nvecLimits(1) nvecLimits(2) 10^timeTicksLimits(1) 10^timeTicksLimits(2)])

set(gca,'Xtick',(10.^(log10(nvecLimits(1)):1:log10(nvecLimits(2)))), ...
  'YTick',(10.^(timeTicksVec)) )

% hold on
% loglog([min(S.nVec(:)) max(S.nVec(:))],...
%   min(S.timeVec(:))*[1 exp(nvecLimits(end)/nvecLimits(1)) ], 'g--')

figSavePathName2 = sprintf(%s_rapid_n_vs_time_Matern_d%d_%s', ...
  S.testFunArg.fName,...
  S.testFunArg.dim, S.timeStamp );
gail.save_image(figH2, 'Paper_cubBayesLattice_g', figSavePathName2)

fprintf('done')



