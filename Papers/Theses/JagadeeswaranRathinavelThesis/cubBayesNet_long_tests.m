%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Workouts for cubBayesNet_g
fprintf('Starting Workouts for cubBayesNet_g')

GAIL_path = GAILstart(0);
logSavePath=strcat([GAIL_path,'OutputFiles',filesep], 'Paper_cubBayesLattice_g');

if exist(logSavePath,'dir')==false
  mkdir(logSavePath);
end

% path to save/read the .mat files
figSavePath = logSavePath;
if exist(figSavePath,'dir')==false
  mkdir(figSavePath);
end

visiblePlot=true;

% temporarily avoid figures to be displayed in matlab
if visiblePlot==false
  set(0,'DefaultFigureVisible','off')
else
  set(0,'DefaultFigureVisible','on')
end

rng(202326) % initialize random number generation to enable reproducability
stopAtTol = true;
alpha = 0.01;
nRepAuto = 100;

% initialize with a template
clear testFunArgs
testFunArgs(9)=struct('fName','','dim',2,'order',1,...
  'sampling','','arbMean',true,'stopCriterion','');
  
% template of input arguments
testFunArgs(1)=struct('fName','MVN','dim',2,'order',1,...
  'sampling','Net','arbMean',true,'stopCriterion','GCV');  % MVN d=3 transformed to Genz d=2
testFunArgs(2)=struct('fName','Keister','dim',4,'order',1,...
  'sampling','Net','arbMean',true,'stopCriterion','GCV');
testFunArgs(3)=struct('fName','optPrice','dim',12,'order',1,...  
  'sampling','Net','arbMean',true,'stopCriterion','GCV');

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
  fprintf('Integrand: %s', fName)
  
  tstart=tic;
  muhatVec = [];
  errVec = [];
  nptsVec = [];
  timeVec = [];
  tolVec = [];
  exitflagVec = [];
  outStructVec = {};
  indx = 1;
  
  if strcmp(fName, 'optPrice')
    % continue
  end
  if strcmp(fName, 'MVN')
    log10ErrVec = -5:1:-2; 
  elseif strcmp(fName, 'Keister')
    log10ErrVec = -4:1:-1; 
  else
    log10ErrVec = -4:1:-1; 
  end
  errTolVecText = arrayfun(@(x){sprintf('1e%d', x)}, log10ErrVec);
  errTolVec = 10.^log10ErrVec;
  sampling = testFunArg.sampling;
  
  for iter=1:length(log10ErrVec)
    
    arbMean=testFunArg.arbMean;
    if arbMean==true
      newPath = strcat(figSavePath, sampling, '/', 'arbMean/');
    else
      newPath = strcat(figSavePath, sampling, '/', 'zeroMean/');
    end
    
	vartx = '';
    dim=testFunArg.dim;
    bern=testFunArg.order;
    
    inputArgs = {'dim',dim, 'order',bern, .... 
      'stopAtTol',stopAtTol, 'stopCriterion',testFunArg.stopCriterion...
      'figSavePath',newPath, 'arbMean',arbMean, 'alpha',alpha ...
      'nRepAuto', nRepAuto, 'log10ErrVec',log10ErrVec,...
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
      if strcmp(fName,'optPrice') % && errTol < 1e-3
        warning('off','GAIL:cubBayesNet_g:maxreached')
        [muhat,err,time,out,tolVec_] = testFun();
        warning('on','GAIL:cubBayesNet_g:maxreached')
      else
        [muhat,err,time,out,tolVec_] = testFun();
      end
      % errVec = [errVec err/errTol];
      errVec = [errVec err./tolVec_];
      
      if quantile(err./tolVec_, 1-alpha/2) > 1 && sum((err./tolVec_) > 1) > 1
        qval = quantile(err./tolVec_, 1-alpha/2);
        fail_count = sum((err./tolVec_) > 1);
        if ~strcmp(fName,'optPrice') || (strcmp(fName,'optPrice') && sum((err./tolVec_) > 1) > 3)
            warning('Error exceeded given threshold: test failed for function %s, qval %f, fail_count %d',...
              fName, qval, fail_count);
        end

        ME = MException('cubBayesNet_g_longtests:errorExceeded', ...
          'Error exceeded given threshold: test failed for function %s',fName);
        % throw(ME)
      end
      
      exitflagVec = [exitflagVec [out.exitflag]'];
      nptsVec = [nptsVec [out.n]'];
      timeVec = [timeVec [out.time]'];
      % tolVec = [tolVec repmat(errTol,size(err))];
      tolVec = [tolVec tolVec_];
      outStructVec{indx} = out;
      muhatVec(indx) = muhat;
      indx = indx + 1;
    end
  end
  
  toc(tstart)
  % suffix this timestamp to all the files stored
  timeStamp = datetime('now','Format','y-MMM-d');
  
  datFileName=sprintf('Sobol_Guaranteed_plot_data_%s_%s_%s_d%d_r%d_%s.mat',...
    fName,stopCrit,vartx,testFunArg.dim,testFunArg.order,timeStamp);
  save([figSavePath filesep datFileName],...
    'errVec','timeVec','tolVec', 'errTolVec','exitflagVec',...
    'outStructVec','testFunArg','log10ErrVec','fName',...
    'timeStamp','figSavePath','nptsVec','stopCrit');
  
  if exist([figSavePath filesep datFileName], 'file') == 2
    guaranteed_plots(figSavePath, {datFileName});
  end
  
end

fprintf('Workouts for cubBayesNet_g: finished')
