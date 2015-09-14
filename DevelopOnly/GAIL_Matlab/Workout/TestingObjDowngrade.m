% Test object inputs
format compact
close all
clearvars

%% Start with a stochProcess and copy and modify it
fred = stochProcess, 
elaine = stochProcess(fred), 
john = stochProcess(elaine,struct('inputType','x')), 
chris = stochProcess(john,struct('timeDim',struct('timeVector',1:5)))

%% Start with a whiteNoise and downgrade to a stocProcess
fjh = whiteNoise
ech = stochProcess(fjh,struct('inputType','x'))
jch = whiteNoise(ech,struct('inputType','n', ...
   'wnParam',struct('distribName','Gaussian')))
cih = whiteNoise(ech,struct('wnParam',struct('sampleKind','Sobol')))

%% Start with a brownMotion
fhick = brownianMotion(struct('timeDim',struct('timeVector',0.1:0.1:1)))
ehick = stochProcess(fhick,struct('inputType','x'))
jhick = brownianMotion(ehick,struct('inputType','n'))

%% Start with a assetParam
derf = assetPath
eniale = brownianMotion(derf)

%% Test plotting
figure
plot(fred)
figure
plot(jch)
figure
plot(jch,'yt.',50,'markersize',100)
figure
plot(jch,'yy',50)
figure
plot(cih,rand(10,3))
try 
   figure
   plot(cih,rand(10,5))
catch
end
figure
plot(jhick,'hist',2000)
figure
plot(derf)
figure
plot(eniale)
   
