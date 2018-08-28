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

stopAtTol = true;
alpha = 0.01;

testFunArgs(2)=struct('fName','MVN','dim',2,'order',2,'varTx','C2sin',...
  'sampling','Lattice','arbMean',true,'stopCriterion','GCV');
testFunArgs(1)=struct('fName','Keister','dim',4,'order',2,'varTx','C1sin',...
  'sampling','Lattice','arbMean',true,'stopCriterion','GCV');
testFunArgs(3)=struct('fName','optPrice','dim',12,'order',1,'varTx','Baker',...
  'sampling','Lattice','arbMean',true,'stopCriterion','GCV');

for i=1:3
  testFunArgs(i+3)=testFunArgs(i);
  testFunArgs(i+3).stopCriterion='full';
end
for i=1:3
  testFunArgs(i+6)=testFunArgs(i);
  testFunArgs(i+6).stopCriterion='MLE';
end

for testFunArg=testFunArgs(1:end)
  
  stopCrit=testFunArg.stopCriterion;
  fName = testFunArg.fName;
  
  if ~strcmp(fName, 'Keister')
    continue
  end
  tstart=tic;
  muhatVec = [];
  errVec = [];
  nptsVec = [];
  timeVec = [];
  tolVec = [];
  outStructVec = {};
  indx = 1;
  if strcmp(fName, 'MVN')
    log10ErrVec = -7:1:-4; 
  elseif strcmp(fName, 'Keister')
    log10ErrVec = -5:1:-2; 
  else
    log10ErrVec = -4:1:-1; 
  end
  errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
  errTolVec = 10.^log10ErrVec;
  sampling = testFunArg.sampling;
  
  for errTol=errTolVec(1:end)
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
    
    inputArgs = {'dim',dim, 'absTol',errTol, 'order',bern, 'ptransform',vartx, ....
      'stopAtTol',stopAtTol, 'stopCriterion',testFunArg.stopCriterion...
      'figSavePath',newPath, 'arbMean',arbMean, 'alpha',alpha ...
      'samplingMethod',sampling, 'visiblePlot',visiblePlot};
    testFun = '';
    switch fName
      case 'Exp(cos)'
        testFun = @()TestExpCosBayesianCubature(inputArgs{:});
      case 'Keister'
        testFun = @()TestKeisterBayesianCubature(inputArgs{:});
      case 'MVN'
        if dim==2 || dim==3
          testFun = @()TestMVN_BayesianCubature(inputArgs{:});
        else
          continue
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
  % timeStamp = datetime('now','Format','d-MMM-y HH-mm-ss');
  timeStamp = datetime('now','Format','y-MMM-d');
  
  datFileName=sprintf('%sGuaranteed_plot_data_%s_%s_%s_d%d_r%d_%s.mat',...
    figSavePath,fName,stopCrit,vartx,testFunArg.dim,testFunArg.order,timeStamp);
  save(datFileName,...
    'errVec','timeVec','tolVec', 'errTolVec',...
    'outStructVec','testFunArg','log10ErrVec','fName',...
    'timeStamp','figSavePath','nptsVec','stopCrit');
  
  if exist(datFileName, 'file') == 2
    guaranteed_plots('', {datFileName});
  end
  
end
diary off


error 'finished'
