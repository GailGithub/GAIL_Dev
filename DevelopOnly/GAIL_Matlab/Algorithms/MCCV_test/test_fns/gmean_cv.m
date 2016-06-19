function YX = gmean_cv(n)
    p.assetParam.initPrice = 11;
    paths = assetPath(p);
    strike = 12;
    
    gbm = genPaths(paths,n);
    meanstock = prod(gbm,2).^(1./paths.timeDim.nSteps);
    
    YX = [max(meanstock - strike, 0) ...
               .* exp(- paths.assetParam.interest .* paths.timeDim.endTime),gbm];
end