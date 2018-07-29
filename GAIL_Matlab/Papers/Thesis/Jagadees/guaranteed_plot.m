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

stopAtTol = true;
alpha = 0.01;

testFunArgs(1)=struct('fName','MVN','dim',2,'order',2','varTx','C0',...
  'sampling','Lattice','arbMean',true,'fullBayes',false,'GCV',false);
testFunArgs(2)=struct('fName','Keister','dim',4,'order',4','varTx','C2sin',...
  'sampling','Lattice','arbMean',true,'fullBayes',false,'GCV',false);
testFunArgs(3)=struct('fName','optPrice','dim',12,'order',2','varTx','Baker',...
  'sampling','Lattice','arbMean',true,'fullBayes',false,'GCV',false);

for i=1:3
  testFunArgs(i+3)=testFunArgs(i);
  testFunArgs(i+3).fullBayes=true;
end
for i=1:3
  testFunArgs(i+6)=testFunArgs(i);
  testFunArgs(i+6).GCV=true;
end

for testFunArg=testFunArgs(1:end)
  
  stopCrit='MLE';
  if testFunArg.fullBayes
    stopCrit='FB';
  end
  if testFunArg.GCV
    stopCrit='GCV';
  end
  
  fName = testFunArg.fName;
%   if strcmp(fName, 'Keister')
%     continue
%   end
%   if strcmp(fName, 'MVN')
%     continue
%   end
  tstart=tic;
  muhatVec = [];
  errVec = [];
  nptsVec = [];
  timeVec = [];
  tolVec = [];
  outStructVec = {};
  indx = 1;
  %pdTx = {'Baker', 'C0',};  % 'C1','C1sin', 'C2sin', };  %, 'none'
  %arbMeanType = [true,false];
  %samplingMethod = {'Lattice',}; %'Sobol',
  %log10ErrVec = -5:1:-2;
  log10ErrVec = -4:1:-1;
  %errTolVecText = {'1e-5','1e-4','1e-3','1e-2',};
  errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
  errTolVec = 10.^log10ErrVec;
  
  sampling = testFunArg.sampling;
  
  for errTol=errTolVec
    errTol;
    
    arbMean=testFunArg.arbMean;
    if arbMean==true
      newPath = strcat(figSavePath, sampling, '/', 'arbMean/');
    else
      newPath = strcat(figSavePath, sampling, '/', 'zeroMean/');
    end
    
    vartx=testFunArg.varTx;
    dim=testFunArg.dim;
    bern=testFunArg.order;
    
    inputArgs = {'dim',dim, 'absTol',errTol, 'order',bern, 'GCV',testFunArg.GCV...
      'ptransform',vartx, 'stopAtTol',stopAtTol, 'fullBayes',testFunArg.fullBayes...
      'figSavePath',newPath, 'arbMean',arbMean, 'alpha',alpha ...
      'samplingMethod',sampling, 'visiblePlot',visiblePlot};
    testFun = '';
    switch fName
      case 'Exp(cos)'
        testFun = @()TestExpCosBayesianCubature(inputArgs{:});
      case 'Keister'
        testFun = @()TestKeisterBayesianCubature(inputArgs{:});
      case 'MVN'
        if dim~=4
          testFun = @()TestMVN_BayesianCubature(inputArgs{:});
        end
      case 'optPrice'
        testFun = @()TestAsianArithmeticMeanOptionAutoExample(inputArgs{:});
      otherwise
        error('Unknown Integrand !');
    end
    
    if ~strcmp(testFun,'')
      [muhat,err,time,out] = testFun();
      
      errVec = [errVec err/errTol];
      if any(err > errTol)
        warning 'Error exceeded given threshold'
      end
      
      nptsVec = [nptsVec [out.n]'];
      timeVec = [timeVec [out.time]'];
      tolVec = [tolVec repmat(errTol,size(err))];
      outStructVec{indx} = out;
      muhatVec(indx) = muhat;
      indx = indx + 1;
    end
  end
  
  toc(tstart)
  timeStamp = datetime('now','Format','d-MMM-y HH-mm-ss');
  
  datFileName=sprintf('%sGuaranteed_plot_data_%s_%s_%s_%s.mat',...
    figSavePath,fName,stopCrit,vartx,timeStamp);
  save(datFileName,...
    'errVec','timeVec','tolVec', 'errTolVec',...
    'outStructVec','testFunArg','log10ErrVec','errTolVecText','fName',...
    'timeStamp','figSavePath','nptsVec','stopCrit');
  
  
  % matFilePath = 'D:\Dropbox\fjhickernellGithub\GAIL_Dev-BayesianCubature\GAIL_Matlab\Papers\Thesis\Jagadees\Paper2018\figures\';
  % load([matFilePath 'Guaranteed_plot_data_Exp(cos)_19-Jul-2018 08-14-22.mat'])
  % timeStamp = '19-Jul-2018 08-14-22';
  %
  % load([matFilePath 'Guaranteed_plot_data_MVN_19-Jul-2018 00-01-41.mat'])
  % timeStamp = '19-Jul-2018 00-01-41';
  %
  % load([matFilePath 'Guaranteed_plot_data_Keister_20-Jul-2018 20-54-43.mat'])
  % timeStamp = '20-Jul-2018 20-54-43';
  
  % for ind=1:length(S.outStructVec)
  %   temp = [S.outStructVec{ind}.n];
  %   nptsVec(ind) = median(temp);
  % end
  
  figHn = figure();
  set(figHn, 'units', 'inches', 'Position', [1 1 9 6])
  errVecLimits = [1E-12, 1E2];
  mvec = outStructVec{1}(1).mvec;
  nptsLimits = [2^(mvec(1)-1), 2^(mvec(end)+1)];
  plot([1, 1], nptsLimits, 'r', 'LineWidth',1)
  hold on
  pointSize=30; %point size
  pointShapes = {'o','s','d','^','v','<','>','p','h'};
  
  for i=1:size(errVec,2)
    scatter(errVec(:,i),nptsVec(:,i),pointSize,log10(tolVec(:,i)),pointShapes{i},'filled')
  end
    
  assert(max(nptsVec(:)) <= nptsLimits(2), ...
    sprintf('nume samples greater than max limit %d', nptsLimits(2)))
  
  set(gca,'xscale','log')
  set(gca,'yscale','log')
  %set(gca,'zscale','log')
  xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
  ylabel('Num. Samples')
  c = colorbar('Direction','reverse', 'Ticks',log10ErrVec, ...
    'TickLabels',errTolVecText, 'TickLabelInterpreter','tex');
  c.Label.Interpreter = 'latex';
  c.Label.String = 'ErrTol, $\varepsilon$';
  % axis tight; not required
  
  axis([errVecLimits(1) errVecLimits(2) nptsLimits(1) nptsLimits(2)])
  set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
    'YTick',(10.^floor(log10(nptsLimits(1)) :1:log10(nptsLimits(2)))))
  if testFunArg.arbMean
    mType = '\(m \neq 0\)'; % arb mean
  else
    mType = '\(m = 0\)'; % zero mean
  end
%   title(sprintf('%s d %d r %d %s %s', testFunArg.fName, ...
%                 testFunArg.dim, testFunArg.order, testFunArg.varTx, mType));
  
  figSavePathName = sprintf('%s%s_guaranteed_npts_%s_%s_%s.png', ...
    figSavePath, fName,stopCrit,vartx,timeStamp );
  saveas(figHn, figSavePathName)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  figH = figure();
  set(figH, 'units', 'inches', 'Position', [1 1 9 6])
  timeLimits = [1E-3, 1E1];
  maxTimeLimit = 2*max(max(timeVec(:)), timeLimits(2));
  plot([1, 1], [timeLimits(1) maxTimeLimit], 'r', 'LineWidth',1)
  hold on
  
  for i=1:size(errVec,2)
    scatter(errVec(:,i),timeVec(:,i),pointSize,log10(tolVec(:,i)),pointShapes{i},'filled')
  end
  
  set(gca,'xscale','log')
  set(gca,'yscale','log')
  xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
  ylabel('Time (secs)')
  c = colorbar('Direction','reverse', 'Ticks',log10ErrVec, ...
    'TickLabels',errTolVecText, 'TickLabelInterpreter','tex');
  c.Label.Interpreter = 'latex';
  c.Label.String = 'ErrTol, $\varepsilon$';
  % axis tight; not required
  
  %assert(max(timeVec(:)) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))
  
  axis([errVecLimits(1) errVecLimits(2) timeLimits(1) maxTimeLimit ])
  set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
    'YTick',(10.^(log10(timeLimits(1)) :2:log10(timeLimits(2)))))
  %title(sprintf('%s d %d r %d %s %s', testFunArg.fName, ...
  %              testFunArg.dim, testFunArg.order, testFunArg.varTx, mType));  
  figSavePathName = sprintf('%s%s_guaranteed_time_%s_%s_%s.png', ...
    figSavePath, fName,stopCrit,vartx,timeStamp );
  saveas(figH, figSavePathName)
  
end

diary off

error 'finished'
