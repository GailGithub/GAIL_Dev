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

%fName = 'Exp(cos)';
fName = 'Keister';
%fName = 'MVN';

tstart=tic;
errVec = [];
timeVec = [];
tolVec = [];
indx = 1;
pdTx = {'C0', 'C1','C1sin', 'C2sin', 'none', 'Baker', };  %
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
    bernOrder=[2 4];
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
        for dim=[2 3 4]
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
            timeVec(indx) = time;
            tolVec(indx) = errTol;
            indx = indx + 1;
            
          end
        end
      end
    end
  end
end
toc(tstart)

figH = figure();
set(figH, 'units', 'inches', 'Position', [1 1 9 6])
errVecLimits = [1E-16, 1E2];
timeLimits = [1E-3, 1E1];
plot([1, 1], timeLimits, 'r', 'LineWidth',1)
hold on
pointSize=30; %point size

pointShapes = {'o','+','x','d'};
offset = length(errVec)/length(errTolVec);
a=1; b=offset;
for i=1:length(errTolVec)
  scatter(errVec(a:b),timeVec(a:b),pointSize,log10(tolVec(a:b)),pointShapes{i})
  a=a+offset; b=b+offset;
end

set(gca,'xscale','log')
set(gca,'yscale','log')
%set(gca,'zscale','log')
xlabel('\(\frac{\vert\mu-\widehat{\mu} \vert}{\varepsilon}\)')
ylabel('Time (secs)')
c = colorbar('Direction','reverse', 'Ticks',log10ErrVec, ...
  'TickLabels',errTolVecText, 'TickLabelInterpreter','tex');
c.Label.Interpreter = 'latex';
c.Label.String = 'ErrTol, $\varepsilon$';
% axis tight; not required

axis([errVecLimits(1) errVecLimits(2) timeLimits(1) timeLimits(2)])
set(gca,'Xtick',(10.^(log10(errVecLimits(1)):4:log10(errVecLimits(2)))), ...
  'YTick',(10.^(log10(timeLimits(1)) :2:log10(timeLimits(2)))))
title(sprintf('Guaranteed cubature : %s', fName));
figSavePathName = sprintf('%s%s guaranteed %s.png', ...
  figSavePath, fName, datetime('now','Format','d-MMM-y HH-mm-ss') );
saveas(figH, figSavePathName)

diary off

error 'finished'
