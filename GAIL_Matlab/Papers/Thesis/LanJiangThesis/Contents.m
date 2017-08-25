% Lanjiangthesis
%
% Files
%
%
%   PlotTwoHumps - This script plots Figure 1.1 in Lan Jiang's thesis
%
%   twohumps - two humps function used in Section 1.1
% 
%   RatioBEnsigtol_fixed_nsig - produces Figure 3.1 in Lan Jiang's thesis
%
%   RatioBEnsigtol_fixed_kmax - produces Figure 3.2 in Lan Jiang's thesis
%
%   BE_CLT_n_ratio_diff_fudge - This script produce Figure 3.3 in Lan
%   Jiang's thesis
%
%   RunTestcubMConProduct - generates the product test function in Section
%   3.5.1 & 4.5.1, and produces Figure 3.4 and 4.2
%
%   producttestfun - produces the product function used to section 3.5.1 in
%   Lan Jiang's thesis.
%
%   productfun - produces the product function used to section 4.5.1 in Lan
%   Jiang's thesis.
%
%   randchoiceproduct - randomly generates the parameters in product test
%   functions in Section 3.5.1 and 4.5.1
%
%   RunTestcubMConGaussiand1 - the Gaussian test function in Section 3.5.2,
%   producing Figure 3.5
%
%   randchoicegaussiand1 - randomly generates the parameters in test functions
%   used by choosetestfun for univariate integration cases.
%
%   RunTestcubMConGaussian- the Gaussian test function in Section 3.5.2,
%   producing Figure 3.6
%
%   randchoicegaussian - randomly generates the parameters in test functions
%   used by choosetestfun for multivariate integration cases.
%
%   TestcubMCgRELTOLDiffSettings - function to test different routines when
%   given a relative error tolerance.
%
%   TestcubMCgDiffSettings - function to test different routines when
%   given an absolute error tolerance.
%
%   RunTestcubMCgonKeister - generates the Keister test function in Section
%   4.5.2 and produces Figure 4.3.
%
%   Keistertrue - computes the true value of the Keister integral in dimension d
%   accuracy might degrade as d increases due to round-off error
%
%   RunTestcubMConMVNP - generates the product test function in Section 4.5.3,
%   and produces Figure 4.4.
%
%   MultivarNorProb - This function computes the cumulative distribution
%   function of the multivariate normal with mean mu, covariance matrix cov
%   and within the region defined by box.
%
%   MVNPexact - This function computes the true solution of multivariate
%   normal probability when the covariance matrix is in a special form:
%   diagonal is 1 and off diagonal elements are all the same.
%
%   choosetestfun - The function chooses and sets up a test function from
%   the parameters input by the user and contained in the structures.
%
%   RunMVNPdeltaplot - produces Figure 4.5 in Lan Jiang's thesis
%
%   PlotmeanMCBernoulli_gResults - generates Figure 5.1 & 5.2
%
%   PlotRatioHoeffCLT -  generates Figure 5.3
%   
%   meanMCBernoulli_g - Monte Carlo to estimate means of Bernoulli random
%   variables
%
%   Test_meanMCBernoulli_g - test function of meanMCBernoulli_g
%
%   DisplayTestResults* - used to plot the examples in Section 3 and 4
%
%   RunTestcubMConGaussian  - This is the Gaussian test function in Section 3.5.2, producing Figure 3.6
%   randchoiceKeister - randchoiceKeister - randomly generates the parameters in Keister test
%   kurtosismax - compute the maximum kurtosis
%   
%   LanJiangThesis.pdf - Thesis
