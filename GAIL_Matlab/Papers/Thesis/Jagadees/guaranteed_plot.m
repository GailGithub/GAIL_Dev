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

for fName={'MVN','Keister','Exp(cos)', }
  
  fName = fName{1};
  tstart=tic;
  muhatVec = [];
  errVec = [];
  nptsVec = [];
  timeVec = [];
  tolVec = [];
  outStructVec = {};
  indx = 1;
  pdTx = {'Baker', 'C0', 'C1','C1sin', 'C2sin', };  %, 'none'
  arbMeanType = [true,false];
  samplingMethod = {'Lattice',}; %'Sobol',
  log10ErrVec = [-5,-4,-3,-2];
  errTolVecText = {'1e-5','1e-4','1e-3','1e-2',};
  errTolVec = 10.^log10ErrVec;
  for sampling=samplingMethod
    sampling = sampling{1};
    
    if strcmp(sampling,'Sobol')
      transforms={'none'}; % no periodization used for Sobol points based algorithm
      bernOrder=[2];
    else
      transforms=pdTx;
      bernOrder=[4 2];
    end
    for errTol=errTolVec
      errTol;
      for arbMean=arbMeanType
        if arbMean==true
          newPath = strcat(figSavePath, sampling, '/', 'arbMean/');
        else
          newPath = strcat(figSavePath, sampling, '/', 'zeroMean/');
        end
        
        for tx=transforms
          vartx=tx{1};
          dimVec=[4 3 2];
          if strcmp(fName,'MVN')
            dimVec=[3 2];
          end
          for dim=dimVec
            for bern=bernOrder
              
              inputArgs = {'dim',dim, 'absTol',errTol, 'order',bern, ...
                'ptransform',vartx, 'stopAtTol',stopAtTol, ...
                'figSavePath',newPath, 'arbMean',arbMean, ...
                'samplingMethod',sampling, 'visiblePlot',visiblePlot};
              switch fName
                case 'Exp(cos)'
                  [muhat,err,time,out] = TestExpCosBayesianCubature(inputArgs{:});
                case 'Keister'
                  [muhat,err,time,out] = TestKeisterBayesianCubature(inputArgs{:});
                case 'MVN'
                  if dim~=4
                    [muhat,err,time,out] = TestMVN_BayesianCubature(inputArgs{:});
                  end
                otherwise
                  error('Unknown Integrand !');
              end
              
              errVec(indx) = err/errTol;
              if errVec(indx) > 1
                error('bug')
              end
              
              nptsVec(indx) = median([out.n]);
              timeVec(indx) = time;
              tolVec(indx) = errTol;
              outStructVec{indx} = out;
              muhatVec(indx) = muhat;
              indx = indx + 1;
              
            end
          end
        end
      end
    end
  end
  toc(tstart)
  timeStamp = datetime('now','Format','d-MMM-y HH-mm-ss');
  
  save(sprintf('%sGuaranteed_plot_data_%s_%s.mat',figSavePath,fName,timeStamp),...
    'errVec','timeVec','tolVec', 'errTolVec',...
    'outStructVec','muhatVec','log10ErrVec','errTolVecText','fName',...
    'timeStamp','figSavePath','nptsVec');
  
  
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
  errVecLimits = [1E-16, 1E2];
  nptsLimits = [2^9, 2^21];
  plot([1, 1], nptsLimits, 'r', 'LineWidth',1)
  hold on
  pointSize=30; %point size
  
  pointShapes = {'o','s','d','^','v','<','>','p','h'};
  offset = length(errVec)/length(errTolVec);
  a=1; b=offset;
  for i=1:length(errTolVec)
    scatter(errVec(a:b),nptsVec(a:b),pointSize,log10(tolVec(a:b)),pointShapes{i},'filled')
    a=a+offset; b=b+offset;
  end
  
  assert(max(nptsVec) <= nptsLimits(2), sprintf('nume samples greater than max limit %d', nptsLimits(2)))
  
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
  %title(sprintf('Guaranteed cubature n : %s', fName));
  
  figSavePathName = sprintf('%s%s guaranteed npts %s.png', ...
    figSavePath, fName, timeStamp );
  saveas(figHn, figSavePathName)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  figH = figure();
  set(figH, 'units', 'inches', 'Position', [1 1 9 6])
  errVecLimits = [1E-16, 1E2];
  timeLimits = [1E-3, 1E1];
  plot([1, 1], timeLimits, 'r', 'LineWidth',1)
  hold on
  pointSize=30; %point size
  
  pointShapes = {'o','s','d','^','v','<','>','p','h'};
  offset = length(errVec)/length(errTolVec);
  a=1; b=offset;
  for i=1:length(errTolVec)
    scatter(errVec(a:b),timeVec(a:b),pointSize,log10(tolVec(a:b)),pointShapes{i},'filled')
    a=a+offset; b=b+offset;
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
  
  assert(max(timeVec) <= timeLimits(2), sprintf('time val greater than max limit %d', timeLimits(2)))
  
  axis([errVecLimits(1) errVecLimits(2) timeLimits(1) timeLimits(2)])
  set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
    'YTick',(10.^(log10(timeLimits(1)) :2:log10(timeLimits(2)))))
  %title(sprintf('Guaranteed cubature : %s', fName));
  
  figSavePathName = sprintf('%s%s guaranteed time %s.png', ...
    figSavePath, fName, timeStamp );
  saveas(figH, figSavePathName)
  
end

diary off

error 'finished'
