%% Generate Examples of American Option Pricing
function [muhat,aMLE,err,outAll] = TestAmericanOptionExample(dim,BernPolyOrder,...
          ptransform,figSavePath,visiblePlot,arbMean)
        
% gail.InitializeWorkspaceDisplay %clean up
format long

nvec = 2.^(7:15)';
nmax = max(nvec);
nRep = 100;
nlarge = nmax*2;
nn = numel(nvec);
alpha = 0.1;

if exist('AmerPutExampleAllData.mat','file')
    load AmerPutExampleAllData
    ArchEuroPut = AmerPut;
    ArchAmerPut = AmerPut;
    ArchAmerPutCV = AmerPut;
    Archnvec = nvec;
end


%% Parameters for the American option
absTol = 1e-3;
relTol = 0;
inp.timeDim.timeVector = 1/16:1/16:1/4; %weekly monitoring for one quarter
%inp.timeDim.timeVector = 1/52:1/52:1/4; %weekly monitoring for one quarter
%inp.timeDim.timeVector = 1/4:1/4:1/4; %weekly monitoring for one quarter
%inp.timeDim.timeVector = 1/8:1/8:1/4; %weekly monitoring for one quarter
inp.assetParam.initPrice = 100; %initial stock price
inp.assetParam.interest = 0.05; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 100; %strike price
inp.priceParam.absTol = absTol; %absolute tolerance of a penny
inp.priceParam.relTol = relTol; %zero relative tolerance
inp.priceParam.cubMethod = 'Sobol'; %Sobol sampling
inp.bmParam.assembleType = 'PCA';
inp.payoffParam.putCallType = {'put'};

%% Construct some different options
EuroPut = optPrice(inp); %construct a European optPrice object
AmerPut = optPrice(EuroPut); %construct an American optPrice object
AmerPut.payoffParam = struct( ...
    'optType',{{'american'}},...
    'putCallType',{{'put'}});
AmerPutCV = optPayoff(EuroPut); %construct an American and European optPayoff object for CV
AmerPutCV.payoffParam = struct( ...
    'optType',{{'american','euro'}},...
    'putCallType',{{'put','put'}});

if 0
%% Construct a very accurate answer
compGold = true;
nGoldRep = 100;
if exist('putPriceExact','var') && ...
        all(ArchAmerPut.timeDim.timeVector == AmerPut.timeDim.timeVector) && ...
        ArchAmerPut.assetParam.initPrice == AmerPut.assetParam.initPrice && ... %initial stock price
        ArchAmerPut.payoffParam.strike == ArchAmerPut.payoffParam.strike, %strike price
    compGold = false;
    disp('Already have gold standard American Put')
end

if compGold
    putPriceGold(nGoldRep,1) = 0;
    nGold = 2^21;
    tic
    for ii = 1:nGoldRep
        gail.TakeNote(ii,10) %print out every 10th ii
        x = net(scramble(sobolset(AmerPut.timeDim.nSteps), ...
            'MatousekAffineOwen'),nGold);
        payoffAmerEuro = genOptPayoffs(AmerPut,x);
        putPriceGold(ii) = mean(payoffAmerEuro(:,1));
    end
    putPriceExact = mean(putPriceGold);
end
disp(['mu  = ' num2str(putPriceExact,15) ' +/- ' num2str(2*std(putPriceGold),10)])
end


%% Bayesian Cubature using Lattice points
compMLE_lattice = true;
if exist('AmerPutExampleAllData.mat','file')
    if exist('muAmerPutMLELattice','var') && all(nvec == Archnvec)
        compMLE_lattice = false;
        disp('Already have Bayesian cubature MLE lattice American Put')
    end
end
if compMLE_lattice
    tic
    %muAmerPutUSobol = zeros(nn,1);
    fun = @(x)genOptPayoffs(AmerPut,x);
    order=2;
    dim=AmerPut.timeDim.nSteps;
    ptransform= 'C1sin';
    testAll=true;
    figSavePath='/home/jagadees/MyWriteup/May4thweek_optprice/';
    fName='optPrice';
    arbMean = false;
    [muAmerPutUSobol, outtemp] = cubMLELattice(fun, ...
        dim, AmerPut.priceParam.absTol, AmerPut.priceParam.relTol,...
        order,ptransform,testAll,figSavePath,fName,arbMean);
    out.nPaths=outtemp.n;
    
    errvecAmerPutMLELattice = abs(putPriceExact - muAmerPutUSobol);
    toc
end







%% IID sampling
AmerPutIID = optPayoff(AmerPut);
AmerPutIID.wnParam = struct('sampleKind','IID','xDistrib','Gaussian');
AmerPutIID.bmParam.assembleType = 'diff';
AmerPutIID.inputType = 'n';
compIID = true;
if exist('AmerPutExampleAllData.mat','file')
    if exist('muAmerPutIID','var') && all(nvec == Archnvec)
        compIID = false;
        disp('Already have IID American Put')
    end
end
if compIID
    tic
    muAmerPutIID = zeros(nn,nRep);
    for i = 1:nRep
        temp = cumsum(genOptPayoffs(AmerPutIID,nmax));
        muAmerPutIID(:,i) = temp(nvec)./nvec;
    end
    errvecAmerPutIID = abs(putPriceExact - muAmerPutIID);
    errmedAmerPutIID = median(errvecAmerPutIID,2);
    errtopAmerPutIID = quantile(errvecAmerPutIID,1-alpha,2);
    toc
end

%% Unscrambled Sobol sampling
compUSobol = true;
if exist('AmerPutExampleAllData.mat','file')
    if exist('muAmerPutUSobol','var') && all(nvec == Archnvec)
        compUSobol = false;
        disp('Already have unscrambled Sobol American Put')
    end
end
if compUSobol
    tic
    muAmerPutUSobol = zeros(nn,1);
    x = net(sobolset(AmerPut.timeDim.nSteps),nmax);
    temp = cumsum(genOptPayoffs(AmerPut,x));
    muAmerPutUSobol(:) = temp(nvec)./nvec;
    errvecAmerPutUSobol = abs(putPriceExact - muAmerPutUSobol);
    toc
end

%% Scrambled Sobol sampling
compSobol = true;
if exist('AmerPutExampleAllData.mat','file')
    if exist('muAmerPutSobol','var') && all(nvec == Archnvec)
        compSobol = false;
        disp('Already have scrambled Sobol American Put')
    end
end
if compSobol
    tic
    muAmerPutSobol = zeros(nn,nRep);
    for i = 1:nRep
        x = net(scramble(sobolset(AmerPut.timeDim.nSteps), ...
            'MatousekAffineOwen'),nmax);
        temp = cumsum(genOptPayoffs(AmerPut,x));
        muAmerPutSobol(:,i) = temp(nvec)./nvec;
    end
    errvecAmerPutSobol = abs(putPriceExact - muAmerPutSobol);
    errmedAmerPutSobol = median(errvecAmerPutSobol,2);
    errtopAmerPutSobol = quantile(errvecAmerPutSobol,1-alpha,2);
    toc
end



%% Save output
save AmerPutExampleAllData.mat
return

f = @(x) genOptPayoffs(AmerPut,x);
AmerOptFourierCoeffDecay(f,AmerPut.timeDim.nSteps)

AmerPutCV = optPayoff(EuroPut); %construct an American and European optPayoff object for CV
AmerPutCV.payoffParam = struct( ...
    'optType',{{'american','euro'}},...
    'putCallType',{{'put','put'}});
f.func =@(x) genOptPayoffs(AmerPutCV,x);
f.cv = AmerPutCV.exactPrice(2:end);
d = AmerPutCV.timeDim.nSteps;
n = 2^20;
y = f.func(net(scramble(sobolset(d),'MatousekAffineOwen'),n));
% figure
% plot(y(:,2),y(:,1),'.')
% xlabel('European payoffs')
% ylabel('American payoffs')
% corrmat = corr(y)
% covmat = cov(y);
% beta = covmat(1,2)/covmat(2,2)
% stdErrNotCorrectY = 2*std(y(:,1))/sqrt(n)
% stdErrNotCorrectYCV = 2*std(y(:,1) + beta * (AmerPutCV.exactPrice(2:end) - y(:,2)))/sqrt(n)
AmerPriceCV(nRep,1) = 0;
%outAmerPutCV(nRep,1) = 0;
for ii = 1:nRep
    [AmerPriceCV(ii),outAmerPutCV(ii)] ...
        = cubSobol_american_g(f, ...
        [zeros(1,AmerPutCV.timeDim.nSteps); ones(1,AmerPutCV.timeDim.nSteps)], ...
        absTol, relTol);
end
comparePrice = [AmerPrice AmerPriceCV]
range(comparePrice,1)
range(comparePrice(:))
compareN = [outAmerPut(:).nPaths outAmerPutCV(:).n]


