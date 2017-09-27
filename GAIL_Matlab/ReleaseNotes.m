% Release Notes 
%
% GAIL Version 2.2, 2017.
% 
% (1) Major changes in algorithms
% 
% We have a new algorithm called meanMC_CLT. In addition, all algorithms in
% version 2.1 are improved with reduced computational complexity and some
% new features. For example, funappx_g and funmin_g are locally adaptive;
% our three multidimensional integration algorithms (cubMC_g, cubLattice_g,
% and cubSobol_g) may take the probability measure uniform sphere in
% addition to uniform hyperbox and Gaussian. The directory Algorithms is
% reorganized with new subfolders to group algorithms and object classes by
% their shared functionalities. We deprecated funappx_g, funmin_g, and
% integral_g in version 2.1 and renamed them in this version as
% funappxPenalty_g, funminPenalty_g, and integralPenalty_g, respectively.
% Lastly, meanMCBer_g in version 2.1 is removed.
% 
% 
% (2) Major changes in publications
% 
% In the folder "Papers", we include currently working and recently
% published research articles related to our core algorithms. In addition,
% we have included nine PhD or MS theses relevant to our research areas.
% For theses and published papers, we included the PDF files of the
% articles as well as MATLAB scripts for reproducing key figures, tables,
% and numerical results in the papers.
% 
% 
% (3) Major changes in documentation
% 
% We continue to release both PDF and searchable HTML documentation for
% GAIL in MATLAB.  In addition, we have added more than a half dozen demos
% for our algorithms.
% 
% (4) Major changes in tests
% 
% We continue to conduct automated nightly tests and weekly long tests on
% our own server.
