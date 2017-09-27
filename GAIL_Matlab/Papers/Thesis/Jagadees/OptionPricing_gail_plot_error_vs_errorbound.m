% plots Error vs Error bound of Bayesian cubature for test funtion AsianArithmeticMeanOptionAutoExample


!synclient HorizEdgeScroll=0 HorizTwoFingerScroll=0 # disable horizontal scrolling
gail.InitializeWorkspaceDisplay %initialize the workspace and the display parameters
format long

visiblePlot = false;
testAll = false;

figSavePath = '/home/jagadees/MyWriteup/Sep_2ndweek/';

if exist(figSavePath,'dir')==false
  mkdir(figSavePath);
end

errVecAll = {};
outAll = {};

absTolValues = [1E-2, 1E-3, 1E-4];

tstart=tic;
pdTx = {'Baker'}; % {'C1','C1sin', 'C2sin', 'C0', 'none', 'Baker'};
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
        for absTol = absTolValues
          [muhat,aMLE,errVecAll{end+1},outAll{end+1}] = TestAsianArithmeticMeanOptionAutoExample(...
            dim,bern,tx{1},newPath,visiblePlot,arbMean,testAll,absTol);
        end
      end
    end
  end
end
toc(tstart)

errVecForPlot = [];
errBdVecForPlot = [];
for i = 1:length(errVecAll)
  errVecForPlot = [errVecForPlot ;errVecAll{i}];
  for j = 1:length(outAll{i})
    errBdVecForPlot = [errBdVecForPlot; outAll{i}{j}.ErrBd];
  end
end


fName = 'AsianArithmeticMeanOption';

figSavePathName = sprintf('%s%s error_vs_errbd 1.png', figSavePath, fName);
plot_error_vs_errbd(errVecForPlot, errBdVecForPlot, fName, ...
  figSavePathName, 10.^[-8, -1], 10.^[-8, -1], visiblePlot)

mat_file_name = sprintf('%s error_vs_errbd data.mat', fName);
save(mat_file_name, 'absTolValues', 'errVecForPlot', 'errBdVecForPlot')
