% Scripts of the three examples in MCQMC paper

% Folders
%
% MCQMC_Paper_Scripts

% Files
%
% AdaptiveMCProbErrAnal_rev_m.pdf --- Monte Carlo and Quasi-Monte Carlo Methods 2012 paper 

% Scripts
%
% choosetestfun.m --- choose and set up a test function
% cubMC.m --- main algorithm to evaluates a multidimensional integral
% cubMCparam.m --- parameter checking
% cubMCerr.m --- handle errors in cubMC and cubMCparam
% DisplayTestResults_BlacknColor.m --- display the plots in black and color format
% geomMeanAsianCall.m --- the test function of geometric asian mean option
% plotTestcubMCblack.m --- plot the mat file in black
% plotTestcubMCcolor.m --- plot the mat file in color
% randchoicegaussian.m --- the function to random choose paramters in
% gaussian test function with multiple dimensions
% randchoicegaussiand1.m --- the function to random choose paramters in
% gaussian test function with multiple dimension one.
% randchoiceGeoCall.m --- the function to random choose paramters in
% geometric asian mean option test function
% RunTestcubMConGaussian.m --- driver file to run the test on gaussian test
% function with multiple dimension
% RunTestcubMConGaussiand1.m --- driver file to run the test on gaussian test
% function with dimension one
% RunTestcubMConGeoAsianCall.m --- driver file to run the test on geometric
% asian mean option test function
% TestcubMCDiffSettings.m --- using diffent random sampling methods to test
% verifyparam.m --- Make sure the parameters defining the functions exist
% and are the right sizes

% Mat file
%
% TestcubMCon-gaussian-uniform-Out-17-Aug-2012_13.12.36N500d1tol0.001.mat
% --- the numerical integration results when using the gaussian test function