% plots Error vs time bound of Bayesian cubature for test funtion MVN


!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

visiblePlot = true;
testAll = false;

figSavePath = '/home/jagadees/MyWriteup/Sep_2ndweek/';

if exist(figSavePath,'dir')==false
    mkdir(figSavePath);
end

errVecAll = {};
outAll = {};

tstart=tic;
pdTx = {'C1sin', 'C2sin', 'Baker'}; % {'C1','C1sin', 'C2sin', 'C0', 'none', 'Baker'};
arbMeanType = [true];  %[true,false]
for arbMean=arbMeanType
  if arbMean==true
    newPath = strcat(figSavePath, 'arbMean/');
  else
    newPath = strcat(figSavePath, 'zeroMean/');
  end
  
  for tx=pdTx
    for dim=[2 3]
      for bern=[4 2]
        [muhat,aMLE,errVecAll{end+1},outAll{end+1}] = TestMVN_BayesianCubature(...
          dim,bern,tx{1},newPath,visiblePlot,arbMean,testAll);
      end
    end
  end
end
toc(tstart)

errVecForPlot = [];
timeVecForPlot = [];
for i = 1:length(errVecAll)
  errVecForPlot = [errVecForPlot errVecAll{i}];
  for j = 1:length(outAll{i})
    timeVecForPlot = [timeVecForPlot; outAll{i}{j}.time];
   end
end

fName = 'MVN';
errorTol = 1e-4;

figSavePathName = sprintf('%s%s errorVtime 1.png', figSavePath, fName);
plot_error_vs_time(errVecForPlot, timeVecForPlot, errorTol, fName, ...
  figSavePathName, 10.^[-3, 0], 10.^[-15, 2], visiblePlot) 


