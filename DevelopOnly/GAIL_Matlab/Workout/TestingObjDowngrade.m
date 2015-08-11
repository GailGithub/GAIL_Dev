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
jch = whiteNoise(fjh,struct('wnParam',struct('sampleKind','Sobol')))
cih = whiteNoise(ech,struct('wnParam',struct('sampleKind','Sobol')))

%% Test plotting
figure
plot(fred)
figure
plot(jch)
figure
plot(cih)
