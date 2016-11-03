m=assetPath
mpaths=genPaths(m,1e5);
mean(mpaths) ...
   - m.assetParam.initPrice*exp(m.timeDim.timeVector ...
   *m.assetParam.interest)

in.assetParam.drift=0.2;
mdrift=assetPath(in)

mdriftpaths=genPaths(mdrift,1e5);
mean(mdriftpaths) ...
   - mdrift.assetParam.initPrice*exp(mdrift.timeDim.timeVector ...
   *(mdrift.assetParam.interest ...
   +mdrift.assetParam.volatility*mdrift.assetParam.drift))

m=optPayoff
mpayoffs=genOptPayoffs(m,1e5);
mean(mpayoffs)

mdrift=optPayoff(in)

mdriftpayoffs=genOptPayoffs(mdrift,1e5);
mean(mdriftpayoffs)

