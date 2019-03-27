% MCQMC2012PAPER   
%
% Files
%
%   MCQMC2012GuaranteedFixedWidthConfidenceIntervals.pdf
%                                  - Monte Carlo and Quasi-Monte Carlo Methods 2012 paper 
%   choosetestfun                  - choose and set up a test function
%   cubMC                          - main algorithm to evaluate a
%                                    multidimensional integral
%   cubMCerr                       - handle errors in cubMC and cubMCparam
%   cubMCparam                     - parameter checking
%   DisplayTestResults_BlacknColor - display the plots in black or color
%                                    format using the mat files in dir
%                                    OutputFiles/MCQMC2012PaperOutput, load
%                                    different MAT files to generate
%                                    figures 3-6 in paper
%   FoolAutomaticAlgorithms        - plot fooling function of figure 2 in paper
%   geomMeanAsianCall              - the test function of geometric Asian mean option
%   MCQMC2012Figs                  - plot figure 1 in paper 
%   peakyfunction                  - construct peaky function
%   plotTestcubMCblack             - plot the mat file in black
%   plotTestcubMCcolor             - plot the mat file in color
%   randchoicegaussian             - randomly choose parameters in Gaussian test function with multiple dimensions
%   randchoicegaussiand1           - randomly choose parameters in Gaussian test function with multiple dimension one
%   randchoiceGeoCall              - randomly choose parameters in geometric asian mean option test function
%   RunTestcubMConGaussian         - driver file to run the test on 
%                                    Gaussian test function with dimension
%                                    2-8,to generate a MAT file used to
%                                    plot figure 5
%   RunTestcubMConGaussiand1       - driver file to run the test on
%                                    Gaussian test function with dimension
%                                    one to generate a MAT file used to
%                                    plot figure 3 and 4
%   RunTestcubMConGeoAsianCall     - driver file to run the test on
%                                    geometric Asian mean option test
%                                    function to generate a MAT file used to
%                                    plot figure 6
%   saveMCQMC2012peakyfundir       - same files
%   snooper                        - function record all the x data points
%   TestcubMCDiffSettings          - using different random sampling methods to test
%   verifyparam                    - make sure the parameters defining the functions exist
