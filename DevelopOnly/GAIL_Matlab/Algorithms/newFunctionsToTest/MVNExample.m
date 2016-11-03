%% Generate Examples of Multivariate Normal Probabilities

InitializeWorkspaceDisplay %clean up 
format long

nvec = 2.^(7:20)';
nlarge = nvec(end)*2;
nn = numel(nvec);
C = [4 1 1; 0 1 0.5; 0 0 0.25];
Cov = C'*C
mu = 0;
a = [-6 -2 -2];
b = [5 2 1];
nRep = 100;
alpha = 0.1;

%% First compute a high accuracy answer
nGold = 2^27;
nRepGold = nRep;
MVNProbBest = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nGold, ...
   'errMeth','n','cubMeth','Sobol','intMeth','Genz');
compGold = true;
if exist('MVNProbExampleData.mat','file')
   load MVNProbExampleData
   if all(MVNProbBest.a == MVNProbBestArch.a) && ...
      all(MVNProbBest.b == MVNProbBestArch.b) && ...
      all(MVNProbBest.mu == MVNProbBestArch.mu) && ...
      all(all(MVNProbBest.Cov == MVNProbBestArch.Cov)) && ...
      MVNProbBest.absTol == MVNProbBestArch.absTol && ...
      MVNProbBest.relTol == MVNProbBestArch.relTol && ...
      all(MVNProbBest.n == MVNProbBestArch.n) && ...
      strcmp(MVNProbBest.intMeth,MVNProbBestArch.intMeth) && ...
      strcmp(MVNProbBest.cubMeth,MVNProbBestArch.cubMeth) && ...
      strcmp(MVNProbBest.errMeth,MVNProbBestArch.errMeth) && ...
      nRepGoldArch == nRepGold
      disp('Already have gold standard answer')
      compGold = false;
   end
end
if compGold
   disp('(Re-)computing gold standard answer')
   muBestvec = zeros(1,nRepGold);
   tic 
   for i = 1:nRepGold
      i
      muBestvec(1,i) = compProb(MVNProbBest); 
   end
   toc
   muBest = mean(muBestvec);
   MVNProbBestArch = MVNProbBest;
   nRepGoldArch = nRepGold;
   save MVNProbExampleData MVNProbBestArch muBestvec muBest nRepGoldArch %save important stuff
end
disp(['mu  = ' num2str(muBest,15) ' +/- ' num2str(2*std(muBestvec),10)])

%% Next try IID sampling
MVNProbIIDGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','IID','intMeth','Genz');
muMVNProbIIDGn = zeros(nn,nRep);
tic 
for i = 1:nRep
   muMVNProbIIDGn(:,i) = compProb(MVNProbIIDGn); 
end
errvecMVNProbIIDGn = abs(muBest - muMVNProbIIDGn);
errmedMVNProbIIDGn = median(errvecMVNProbIIDGn,2);
errtopMVNProbIIDGn = quantile(errvecMVNProbIIDGn,1-alpha,2);
toc


%% Then try scrambled Sobol sampling
MVNProbSobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','Sobol','intMeth','Genz');
muMVNProbSobolGn = zeros(nn,nRep);
tic 
for i = 1:nRep
   muMVNProbSobolGn(:,i) = compProb(MVNProbSobolGn); 
end
errvecMVNProbSobolGn = abs(muBest - muMVNProbSobolGn);
errmedMVNProbSobolGn = median(errvecMVNProbSobolGn,2);
errtopMVNProbSobolGn = quantile(errvecMVNProbSobolGn,1-alpha,2);
toc

%% Now try unscrambled Sobol sampling
MVNProbuSobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvec, ...
   'errMeth','n','cubMeth','uSobol','intMeth','Genz');
tic
muMVNProbuSobolGn = compProb(MVNProbuSobolGn); 
errMVNProbuSobolGn = abs(muBest - muMVNProbuSobolGn);
toc

%% Try optimal weights
tic
nvecOPT = 2.^(1:10)';
nnOPT = numel(nvecOPT);
MVNProbOptSobolGn = multivarGauss('a',a,'b',b,'Cov',Cov,'n',nvecOPT, ...
   'errMeth','n','cubMeth','SobolOpt','intMeth','Genz');
muMVNProbOptSobolGn = zeros(nnOPT,nRep);
tic 
for i = 1:nRep
   if i/10 == floor(i/10), i, end
   muMVNProbOptSobolGn(:,i) = compProb(MVNProbOptSobolGn); 
end
errvecMVNProbOptSobolGn = abs(muBest - muMVNProbOptSobolGn);
errmedMVNProbOptSobolGn = median(errvecMVNProbOptSobolGn,2);
errtopMVNProbOptSobolGn = quantile(errvecMVNProbOptSobolGn,1-alpha,2);
toc

save MVNProbExampleAllData.mat

