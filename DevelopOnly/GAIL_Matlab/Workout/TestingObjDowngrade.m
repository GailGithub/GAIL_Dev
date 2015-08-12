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
fhick = brownianMotion
ehick = stochProcess(fhick,struct('inputType','x'))
jhick = brownianMotion(ehick,struct('inputType','n'))

%% Test plotting
figure
plot(fred)
figure
plot(jch)
figure
plot(jch,'yt.',50)
figure
plot(cih,rand(10,3))
try 
   figure
   plot(cih,rand(10,5))
catch
end
figure
plot(jhick,'yth',200)
   
