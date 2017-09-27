% plots Error vs time of Bayesian cubature for test funtion AsianArithmeticMeanOptionAutoExample

!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

visiblePlot = true;

figSavePath = '/home/jagadees/MyWriteup/Sep2ndweek_optprice/';

if exist(figSavePath,'dir')==false
    mkdir(figSavePath);
end

errVecAll = {};
outAll = {};

tstart=tic;
pdTx = {'C1sin', 'Baker'}; % {'C1','C1sin', 'C2sin', 'C0', 'none', 'Baker'};
arbMeanType = [true,false];
for arbMean=arbMeanType
  if arbMean==true
    newPath = strcat(figSavePath, 'arbMean/');
  else
    newPath = strcat(figSavePath, 'zeroMean/');
  end
  
  for tx=pdTx
    for dim=[2]
      for bern=[4 2]
        [muhat,aMLE,errVecAll{end+1},outAll{end+1}] = TestAsianArithmeticMeanOptionAutoExample(dim,bern,tx{1},newPath,visiblePlot,arbMean);
      end
    end
  end
end
toc(tstart)

errVecForPlot = [];
timeVecForPlot = [];
for i = 1:length(errVecAll)
  errVecForPlot = [errVecForPlot; errVecAll{i}];
  for j = 1:length(outAll{i})
    timeVecForPlot = [timeVecForPlot; outAll{i}{j}.time];
   end
end

fName = 'AsianArithmeticMeanOption';
errorTol = 1e-2;

figSavePathName = sprintf('%s%s errorVtime.png', ...
  figSavePath, fName);
plot_error_vs_time(errVecForPlot, timeVecForPlot, errorTol, fName, figSavePathName, visiblePlot)


