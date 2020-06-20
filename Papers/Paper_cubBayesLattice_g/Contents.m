% PAPER_CUBBAYESLATTICE_G
%
% Files
%
%   cubBayesLattice_guaranteed_plots - produces the error plots used in the paper.
%      Also saves the required data for plots in .mat files.
%
%   cubBayesLattice_long_tests - master script to call the workout tests. 
%      Thows exception if the integration error is higher than expected.
%
%   cubBayesMLE_Matern_g - Implements the basic Bayesian cubature algorithm 
%      without any speedup achieved in the paper, used to compare the speedup.
%
%   GenzFunc - Implements Genz function.
%
%   guaranteed_plots - Generates plots from the save .mat files.
%   IntegrandPlots - Produces other plots used in the paper.
%
%   keisterFunc - Implements the Keister function.
%
%   matern_guaranteed_plots - produces the error plots for the integration example
%           using cubBayesMLE_Matern_g
%
%   plot_fourier_kernel - generates the plots of the shift invariant kenerl
%
%   TestAsianArithmeticMeanOptionAutoExample - main script to estimate Asian 
%          arithmetic option price
%
%   TestKeisterBayesianCubature - main script to integrate Keister function
%
%   TestMVN_BayesianCubature - main script to integrate Multivariate normal probability
%   
%	   
%   Fast automatic Bayesian cubature using Lattice sampling.pdf - Paper 
