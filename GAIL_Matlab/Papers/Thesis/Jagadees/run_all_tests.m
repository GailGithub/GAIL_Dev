
!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

figSavePath = '/home/jagadees/MyWriteup/Sep_2ndweek/';

if exist(figSavePath,'dir')==false
    mkdir(figSavePath);
end

% log the results
completereport = strcat(figSavePath,...
    '_tests-logs-', datestr(now,'yyyy-mm-dd-HH-MM-SS'),'.txt');
diary(completereport)

visiblePlot=false;

%
% https://www.mathworks.com/matlabcentral/answers/98969-how-can-i-temporarily-avoid-figures-to-be-displayed-in-matlab
%
if visiblePlot==false
  set(0,'DefaultFigureVisible','off')
else
  set(0,'DefaultFigureVisible','on')
end

if false
  [muhat,aMLE,err,out] = TestExpCosBayesianCubature(2,2,'none',...
    strcat(figSavePath, 'zeroMean/'),visiblePlot,false)
  
  [muhat,aMLE,err,out] = TestMVN_BayesianCubature(2,2,'Baker',...
    strcat(figSavePath, 'zeroMean/'),visiblePlot,false)
  fprintf('done')
end

stopAtTol = false;

tstart=tic;
pdTx = {'C1','C1sin', 'C2sin', 'C0', 'none', 'Baker'};
arbMeanType = [true,false];
for arbMean=arbMeanType
  if arbMean==true
    newPath = strcat(figSavePath, 'arbMean/');
  else
    newPath = strcat(figSavePath, 'zeroMean/');
  end
  for tx=pdTx
    for dim=[2 3 4]
      for bern=[4 2]
        TestExpCosBayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol)
        TestKeisterBayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol)
        if dim~=4
          TestMVN_BayesianCubature(dim,bern,tx{1},newPath,visiblePlot,arbMean,stopAtTol)
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
