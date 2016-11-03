function []= testingCovariance ()

   clearvars
format short
n = 1e8;
X = gpuArray.rand(n,1);
d = 10;
X = reshape(X,n/d,d);
meanX = mean(X)
covX = cov(X)
        
    

end
