
if isunix
  !synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
end
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

if isunix
  figSavePath = '/home/jagadees/MyWriteup/';
else
  figSavePath = 'D:/MyWriteup/';
end
figSavePath = strcat(figSavePath, 'Dec_3rdweek/');

if exist(figSavePath,'dir')==false
  mkdir(figSavePath);
end

% log the results
completereport = strcat(figSavePath,...
  '_tests-logs-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

visiblePlot=true;

%
% https://www.mathworks.com/matlabcentral/answers/98969-how-can-i-temporarily-avoid-figures-to-be-displayed-in-matlab
%
if visiblePlot==false
  set(0,'DefaultFigureVisible','off')
else
  set(0,'DefaultFigureVisible','on')
end

if false
  %sampling='Lattice';
  sampling='Sobol';
  figSavePath = strcat(figSavePath, sampling, '/');
  [muhat,aMLE,err,out] = TestExpCosBayesianCubature(1,2,'none',...
    strcat(figSavePath, 'zeroMean/'),visiblePlot,false,false,sampling)
  
  [muhat,aMLE,err,out] = TestExpCosBayesianCubature(2,2,'none',...
    strcat(figSavePath, 'zeroMean/'),visiblePlot,false,false,sampling)
  
  [muhat,aMLE,err,out] = TestMVN_BayesianCubature(2,2,'Baker',...
    strcat(figSavePath, 'zeroMean/'),visiblePlot,false,false,sampling)
  fprintf('done')
end

stopAtTol = false;

tstart=tic;
pdTx = {'C1','C1sin', 'C2sin', 'C0', 'none', 'Baker'};
arbMeanType = [true,false];
samplingMethod = {'Lattice', 'Sobol', };
for sampling=samplingMethod
  for arbMean=arbMeanType
    if arbMean==true
      newPath = strcat(figSavePath, sampling{1}, '/', 'arbMean/');
    else
      newPath = strcat(figSavePath, sampling{1}, '/', 'zeroMean/');
    end

    if strcmp(sampling,'Sobol')
      transforms={'none'};
    else
      transforms=pdTx;
    end
    for tx=transforms
      for dim=[2 3 4]
        for bern=[4 2]
          if dim~=4
            TestMVN_BayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol,sampling)
          end
          TestKeisterBayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol,sampling)
          TestExpCosBayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol,sampling)
        end
      end
    end
  end
end
toc(tstart)


if 0
  cp MyWriteup/Apr1stweek/MVN/C2sin/* Dropbox/writeup/BeamerPresent/figures/ &&
  cp MyWriteup/Apr1stweek/Keister/C2sin/* Dropbox/writeup/BeamerPresent/figures/ &&
  cp MyWriteup/Apr1stweek/Exp\(cos\)/C1sin/* Dropbox/writeup/BeamerPresent/figures/
  
  cp MyWriteup/Apr1stweek/Exp\(cos\)/C2sin/* Dropbox/writeup/BeamerPresent/figures/
end

diary off

error 'finished'
