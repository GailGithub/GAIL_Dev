%% M.S. Thesis of IIT at July,2016
% --------------------------------
% Title: Reliable Quasi-Monte Carlo with Control Variates
% Author: Da Li   
% Contact: dli37@hawk.iit.edu 
% Github: https://github.com/TatsuLee/IIT_MS_Thesis
% --------------------------------
% 
%% Files
% --------------------------------
% geoAsianAcc.m        - accuracy test of gmean Asian call option pricing with CV 
%                        generate results in table 4.2 (page 20)
% geoAsianEff.m        - efficiency test I of geometric mean Asian call option pricing with poor CV 
%                        generate results in table 4.4 (page 21) 
% ariAsianEff.m        - efficiency test II of arithematic mean Asian call option pricing with good CV 
%                        generate results in table 4.5 (page 22)
% ariAsianWalshDecay.m - Walsh coefficents decay of amean Asian option payoff with/without CV  
%                        plot figure 4.1 (page 23)
% barrierEff.m         - efficiency test II of up and in barrier call option pricing with good CV 
%                        generate results in table 4.6 (page 24)
% mvnProbEff.m         - efficiency test of multivariate normal probability  
%                        generate results in table 4.7 (page 27)
% --------------------------------
% 
%% Packages required 
% --------------------------------
% ALL option pricing tests require the 'Monte_Carlo_Finance' packcage.
% The package can be found at github:
% https://github.com/GailGithub/GAIL_Dev/tree/develop/DevelopOnly/GAIL_Matlab/Algorithms/Monte_Carlo_Finance
% Please add package folder to your matlab path before runing the tests.
