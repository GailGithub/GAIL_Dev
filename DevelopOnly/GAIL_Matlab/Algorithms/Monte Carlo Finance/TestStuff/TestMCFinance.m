%Testing Monte Carlo Finance stuff
format compact
clearvars

%% Setting up parameters
inp.inputType = 'n';
inp.timeDim.timeVector = 0.25:0.25:1;
inp.assetParam.volatility = 0.5;
inp.assetParam.initPrice = 100;
inp.payoffParam.strike = 100;
inp.priceParam.absTol = 0.1;
inp.priceParam.relTol = 0.01;
priceObj = optPrice(inp);
disp(priceObj)

optTypeVec = {'euro','euro','gmean','amean'};
putCallTypeVec = {'call','put','call','call'};
nOptType = numel(optTypeVec);
cubMethVec = {'IID_MC', 'IID_MC_abs', 'IID_MC_new', 'IID_MC_newtwo', 'Sobol', 'lattice'};
ncubMeth = numel(cubMethVec);
errfac = NaN(ncubMeth,nOptType);
time = errfac;
nPaths = errfac;
exactPrice = NaN(1,nOptType);

for ii = 1:ncubMeth
   priceParam.cubMethod = cubMethVec{ii};
   priceObj.priceParam = priceParam;
   disp(cubMethVec{ii});
   if priceObj.priceParam.absTol > 0 || ~strcmp(cubMethVec{ii}, 'IID_MC_abs')
      for jj = 1:nOptType
         payoffParam.optType = cellstr(optTypeVec{jj});
         payoffParam.putCallType = cellstr(putCallTypeVec{jj});
         priceObj.payoffParam = payoffParam;
         if ii == 1
            exactPrice(jj) = priceObj.exactPrice;
         end
         [price, out] = genOptPrice(priceObj);
         errfac(ii,jj) = abs(priceObj.exactPrice - price) ...
            /max(inp.priceParam.absTol,inp.priceParam.relTol*abs(priceObj.exactPrice));
         time(ii,jj) = out.time;
         nPaths(ii,jj) = out.nPaths;
      end     
   end
end

%% Display results
fprintf('             ')
for jj=1:nOptType
   fprintf('%12s',optTypeVec{jj})
end
fprintf('\n             ')
for jj=1:nOptType
   fprintf('%12s',putCallTypeVec{jj})
end
fprintf('\nExact Price    ')
fprintf('%12.6f',exactPrice)
fprintf('\n\nTime in seconds')
for ii = 1:ncubMeth
   fprintf('\n%15s',cubMethVec{ii})
   fprintf('%12.6f',time(ii,:))
   %fprintf('\n')
end
fprintf('\n\nTrue Error/Tolerance')
for ii = 1:ncubMeth
   fprintf('\n%15s',cubMethVec{ii})
   fprintf('%12.6f',errfac(ii,:))
   %fprintf('\n')
end
fprintf('\n\nNumber of Paths')
for ii = 1:ncubMeth
   fprintf('\n%15s',cubMethVec{ii})
   fprintf('%12.4g',nPaths(ii,:))
   %fprintf('\n')
end
fprintf('\n')
