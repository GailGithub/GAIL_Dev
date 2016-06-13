addpath('/Users/qiantianpei/Desktop/GAIL_Dev/DevelopOnly/GAIL_Matlab/Algorithms')

% distance 
distfun2 = @(n) sqrt(sum((rand(n,2)  - rand(n,2)).^2,2));

tic, tmu = meanMCCV_g(@distfun,repmat(0.5,1,4), distfun2, 2e-4, 0), toc
tic, tmu = meanMC_g(distfun2,2e-4,0), toc

tic, tmu = meanMCCV_g('YXrand',@distfun, 'muX',[0.5 0.5 0.5 0.5], ...   
    'Yrand',distfun2, 'abstol',0.0002, 'reltol',0), toc

% exponential integral 
expr2 = @(n) exp(rand(n, 1));
tic, tmu = meanMCCV_g(@expr,[0.5],expr2, 0.0002, 0), toc
tic, tmu = meanMC_g(expr2,0.0002,0), toc

% Asian geometric call
r = 0.01;
initP = 11;
muX = [exp(r*1)*initP exp(r*2)*initP exp(r*3)*initP];

tic, tmu = meanMCCV_g(@gmean_cv,muX,@gmean, 0.01, 0), toc
tic, tmu = meanMC_g(@gmean,0.01,0), toc;

% exact price
inp.payoffParam.optType = {'gmean'};
inp.payoffParam.putCallType = {'call'};
inp.payoffParam.strike = 12;

inp.assetParam.initPrice = 11;
inp.priceParam.relTol = 0;
inp.priceParam.absTol = 0.01;

GMean = optPrice(inp);

GMean.exactPrice