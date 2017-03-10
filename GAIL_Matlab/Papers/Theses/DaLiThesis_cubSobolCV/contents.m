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
% table4.2 @ P20 -- geoAsianAcc.m 
%	         -- accuracy test of geometric mean Asian call option pricing with CV 
% table4.4 @ P21 -- geoAsianEff.m 
%		 -- efficiency test I of geometric mean Asian call option pricing with poor CV 
% table4.5 @ P22 -- ariAsianEff.m 
%		 -- efficiency test II of arithmetic mean Asian call option pricing with good CV 
% table4.6 @ P24 -- barrierEff.m 
%		 -- efficiency test II of up and in barrier call option pricing with good CV 
% table4.7 @ P27 -- mvnProbEff.m 
%		 -- efficiency test of multivariate normal probability  
% figure2.1 @ P5  -- sobolPlot.m 
%		  -- plots of 256 random points and Sobol points in [0,1)^2 
% figure2.2 @ P9  -- conePlot.m 
%		  -- plots of Walsh coefficients in cone conditions 
% figure4.1 @ P23 -- walshDecayCVPlot.m 
%		  -- Walsh coefficients decay of arithmetic mean Asian option payoff with/without CV  
% --------------------------------
% 
%% Packages requirements 
% --------------------------------
% ALL option pricing tests require the 'Monte_Carlo_Finance' package.
% The package can be found at github:
% https://github.com/GailGithub/GAIL_Dev/tree/develop/DevelopOnly/GAIL_Matlab/Algorithms/Monte_Carlo_Finance
% Please add package folder to your Matlab path before running the tests.
