%% Guaranteed Automatic Integration Library (GAIL) Version 2.3.1 Documentation
%
% GAIL (Guaranteed Automatic Integration Library) is created, developed,
% and maintained by Fred Hickernell (Illinois Institute of Technology),
% Sou-Cheng Choi (IIT and Kamakura Corporation), Yuhan Ding (IIT), Lan
% Jiang (IIT alum), Jagadeeswaran Rathinavel (IIT), Lluis Antoni Jimenez
% Rugama (IIT alum), Xin Tong (University of Illinois at Chicago), Kan
% Zhang (IIT), Yizhi Zhang (IIT alum), and Xuan Zhou (IIT alum).
%
% GAIL is a suite of algorithms for integration problems in one, many, and
% infinite dimensions, and whose answers are guaranteed to be correct.
%
%% Introduction
%
% <html>
% <a href="help_license.html">GAIL License</a>
% <p><a href="help_readme.html">README</a>
% <p><a href="help_ReleaseNotes.html">Release Notes</a>
% </html>
%
%% Functions
%
% <html>
%    <a href="help_funappx_g.html">funappx_g</a>
% <p><a href="help_integral_g.html">integral_g</a>
% <p><a href="help_funmin_g.html">funmin_g</a>
% <p><a href="help_meanMC_g.html">meanMC_g</a>
% <p><a href="help_cubMC_g.html">cubMC_g</a>
% <p><a href="help_cubSobol_g.html">cubSobol_g</a>
% <p><a href="help_cubLattice_g.html">cubLattice_g</a>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% <p><a href="help_cubBayesNet_g.html">cubBayesNet_g</a>
% <p><a href="help_meanMC_CLT.html">meanMC_CLT</a>
% </html>
%
%% Demos
%
% <html>
% <a href="demo_funappx_g.html">Demos for funappx_g</a>
% <p><a href="demo_funmin_g.html">Demos for funmin_g</a>
% <p><a href="demo_integral_g.html">Demos for integral_g</a>
% <p><a href="demo_meanMC_g.html">Demos for meanMC_g</a>
% <p><a href="demo_cubMC_g.html">Demos for cubMC_g</a>
% <p><a href="demo_cubSobol_g.html">Demos for cubSobol_g</a>
% <p><a href="demo_cubLattice_g.html">Demos for cubLattice_g</a>
% <p><a href="demo_cubBayesLattice_g.html">Demos for cubBayesLattice_g</a>
% <p><a href="demo_cubBayesNet_g.html">Demos for cubBayesNet_g</a>
% <p><a href="demo_meanMC_CLT.html">Demos for meanMC_CLT</a>
% </html>
%
%% Website
% For more information about GAIL, visit
% <http://gailgithub.github.io/GAIL_Dev/ Gailteam>
%



%% GAIL License
% Copyright © 2020, Illinois Institute of Technology. All rights reserved.
%  
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%   * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%     
%   * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in the 
%     documentation and/or other materials provided with the distribution.
%     
%   * Neither the name of Illinois Institute of Technology nor the names of
%     its contributors may be used to endorse or promote products derived 
%     from this software without specific prior written permission.
%  
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS 
% "AS IS" AND WITHOUT ANY WARRANTY OF ANY KIND, WHETHER EXPRESS, IMPLIED, 
% STATUTORY OR OTHERWISE, INCLUDING WITHOUT LIMITATION WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR USE AND NON-INFRINGEMENT, ALL 
% OF WHICH ARE HEREBY EXPRESSLY DISCLAIMED. MOREOVER, THE USER OF THE 
% SOFTWARE UNDERSTANDS AND AGREES THAT THE SOFTWARE MAY CONTAIN BUGS, 
% DEFECTS, ERRORS AND OTHER PROBLEMS THAT COULD CAUSE SYSTEM FAILURES, AND 
% ANY USE OF THE SOFTWARE SHALL BE AT USER'S OWN RISK. THE COPYRIGHT 
% HOLDERS AND CONTRIBUTORS MAKES NO REPRESENTATION THAT THEY WILL ISSUE 
% UPDATES OR ENHANCEMENTS TO THE SOFTWARE.  
%  
% IN NO EVENT WILL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY 
% DIRECT, INDIRECT, SPECIAL, INCIDENTAL, CONSEQUENTIAL, EXEMPLARY OR 
% PUNITIVE DAMAGES, INCLUDING, BUT NOT LIMITED TO, DAMAGES FOR INTERRUPTION 
% OF USE OR FOR LOSS OR INACCURACY OR CORRUPTION OF DATA, LOST PROFITS, OR 
% COSTS OF PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES, HOWEVER CAUSED 
% (INCLUDING BUT NOT LIMITED TO USE, MISUSE, INABILITY TO USE, OR 
% INTERRUPTED USE) AND UNDER ANY THEORY OF LIABILITY, INCLUDING BUT NOT 
% LIMITED TO CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE AND WHETHER OR NOT THE 
% COPYRIGHT HOLDER AND CONTRIBUTORS WAS OR SHOULD HAVE BEEN AWARE OR 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE OR FOR ANY CLAIM ALLEGING 
% INJURY RESULTING FROM ERRORS, OMISSIONS, OR OTHER INACCURACIES IN THE 
% SOFTWARE OR DESTRUCTIVE PROPERTIES OF THE SOFTWARE.  TO THE EXTENT THAT 
% THE LAWS OF ANY JURISDICTIONS DO NOT ALLOW THE FOREGOING EXCLUSIONS AND 
% LIMITATION, THE USER OF THE SOFTWARE AGREES THAT DAMAGES MAY BE 
% DIFFICULT, IF NOT IMPOSSIBLE TO CALCULATE, AND AS A RESULT, SAID USER HAS 
% AGREED THAT THE MAXIMUM LIABILITY OF THE COPYRIGHT HOLDER AND 
% CONTRIBUTORS SHALL NOT EXCEED US$100.00.
%  
% THE USER OF THE SOFTWARE ACKNOWLEDGES THAT THE SOFTWARE IS BEING PROVIDED 
% WITHOUT CHARGE, AND AS A RESULT, THE USER, ACKNOWLEDGING THAT HE OR SHE 
% HAS READ THE SAME, AGREES THAT THE FOREGOING LIMITATIONS AND RESTRICTIONS 
% REPRESENT A REASONABLE ALLOCATION OF RISK.%% Guaranteed Automatic Integration Library (GAIL)
% 
% GAIL Version 2.3.1, 2020.
% 
% See LICENSE.m for copyright and disclaimer. Refer to ReleaseNotes.m for
% what is new in this version.
% 
% GAIL is a suite of algorithms for integration problems in one and many
% dimensions, and whose answers are guaranteed to be correct.
% 
% 
%% Developed by
% 
% 
% Fred Hickernell, Sou-Cheng Choi, and their collaborators including Yuhan
% Ding, Lan Jiang, Lluis Antoni Jimenez Rugama, Da Li, Jagadeeswaran
% Rathinavel, Kan Zhang, Yizhi Zhang, and Xuan Zhou, Department of Applied
% Mathematics, Illinois Institute of Technology (IIT) and Xin Tong,
% Department of Mathematics, Statistics, and Computer Science, University
% of Illinois at Chicago. 
% 
% We thank the contributions of Aleksei Sorokin, Noah Grudowski, Francisco
% Hernandez, Cu Hauw Hung, Yueyi Li, Xincheng Sheng, Xiaoyang Zhao, Tianci
% Zhu, and the IIT classes of SCI 498 Adaptive Monte Carlo Algorithms with
% Applications to Financial Risk Management, Summer 2016; MATH 491 Reading
% & Research, Summer 2015; SCI 498/MATH 491 Computational Social Sciences,
% Summer 2016; MATH 491-195 Solving Problems in the Social Sciences Using
% Tools from Computational Mathematics and Statistics, Summer 2015; Math
% 573 Reliable Mathematical Software, Fall 2013.
% 
% 
% 
% Please cite the following software, papers, and materials:
% 
% Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available from
% http://gailgithub.github.io/GAIL_Dev/
% (this software)
% 
% Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible Research via
% Supportable Scientific Software," Journal of Open Research Software, Volume 2,
% Number 1, e22, pp. 1-7, 2014.
% (describes principles of Reliable Reproducible Research and Supportable
% Scientific Software)
% 
% Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Xin Tong, "Local
% Adaption for Approximation and Minimization of Univariate Functions,"
% Journal of Complexity 40, pp. 17-33, 2017.
% (describes *funappx_g.m* and *funmin_g.m*)
% 
% Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable Mathematical
% Software" [Course Slides], Illinois Institute of Technology, Chicago, IL, 2013.
% Available from http://gailgithub.github.io/GAIL_Dev/
% (develops practices of Reliable Reproducible Research and Supportable
% Scientific Software)
% 
% Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
% (describes deprecated *integralTrap_g.m* and deprecated *funappxtau_g.m*)
% 
% Yuhan Ding, "Guaranteed Adaptive Univariate Function Approximation," PhD
% thesis, Illinois Institute of Technology, 2015.
% (describes deprecated *funappxPenalty_g.m*)
% 
% Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
% "Guaranteed conservative fixed width confidence intervals via Monte
% Carlo sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012
% (J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
% Springer-Verlag, Berlin, pp. 105-128, 2014.
% (describes *meanMC_g.m* and *cubMC_g.m*)
% 
% Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable adaptive
% cubature using digital sequences", Monte Carlo and Quasi-Monte Carlo
% Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D. Nuyens, eds.),
% Springer Proceedings in Mathematics and Statistics, vol. 163,
% Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp. 367-383.
% (describes *cubSobol_g.m*)
% 
% Fred J. Hickernell, Lluis Antoni Jimenez Rugama, and Da Li, "Adaptive
% quasi-Monte Carlo methods, 2017+, submitted for publication, 
% arXiv:1702.01491 [math.NA].
% 
% Yizhi Zhang, Guaranteed Adaptive Automatic Algorithms for 
% Univariate Integration: Methods, Costs and Implementations, PhD
% Thesis, Illinois Institute of Technology, 2018.
% (describes *integral_g.m*)
% 
% Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
% Means of Random Variables, PhD Thesis, Illinois Institute of
% Technology, 2016.
% (describes *meanMC_g.m* and *cubMC_g.m*)
% 
% Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive
% multidimensional integration based on rank-1 lattices," Monte Carlo
% and Quasi-Monte Carlo  Methods: MCQMC, Leuven, Belgium, April 2014
% (R. Cools and D. Nuyens, eds.), Springer Proceedings in Mathematics
% and Statistics, vol. 163, Springer-Verlag, Berlin, 2016, arXiv:1411.1966,
% pp. 407-422.
% (describes *cubLattice_g.m*)
% 
% Fred J. Hickernell, Sou-Cheng T. Choi, Lan Jiang, and Lluis Antoni Jimenez
% Rugama, "Monte Carlo simulation, automatic stopping criteria for", Wiley
% StatsRef: Statistics Reference Online, pp. 1-7, 2018.
% (review of *cubMC_g*, *cubLattice_g*, and *cubSobol_g*)
%
% Jagadeeswaran Rathinavel and Fred J. Hickernell, "Fast automatic Bayesian
% cubature using lattice sampling," Statistics and Computingt 29,
% pp. 1215-1229, 2019. Available from https://doi.org/10.1007/s11222-019-09895-9
% (describes *cubBayesLattice_g*)
% 
% Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using matching
% kernels and designs," PhD thesis, Illinois Institute of Technology, 2019.
% (describes *cubBayesLattice_g* and *cubBayesNet_g*)
% 
% Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
% (discusses practice and challenges in Sustainable Software for Science)
% 
% Daniel S. Katz, Sou-Cheng T. Choi, Nancy Wilkins-Diehr, Neil Chue Hong,
% Colin C. Venters, James Howison, Frank J. Seinstra, Matthew Jones, Karen
% Cranston, Thomas L. Clune,  Miguel de Val-Borro, and Richard Littauer,
% "Report on the Second Workshop on Sustainable Software for Science:
% Practice and Experiences (WSSSPE2)," Journal of Open Research Software,
% Volume 4, Number 1, e7, 2016.
% (discusses practice and challenges in Sustainable Software for Science)
% 
% Daniel S. Katz, Sou-Cheng T. Choi, Kyle E. Niemeyer, James
% Hetherington, Frank Loffler, Dan Gunter, Ray Idaszak, Steven R. Brandt,
% Mark A. Miller, Sandra Gesing, Nick D. Jones, Nic Weber, Suresh Marru,
% Gabrielle Allen, Birgit Penzenstadler, Colin C. Venters, Ethan Davis,
% Lorraine Hwang, Ilian Todorov, Abani Patra, and Miguel de Val-Borro,
% "Report on the Third Workshop on Sustainable Software for Science:
% Practice and Experiences (WSSSPE3)," Journal of Open Research Software,
% Volume 4, Number 1, e37, 2016.
% (discusses practice and challenges in Sustainable Software for Science)
% 
% Da Li, "Reliable Quasi-Monte Carlo with Control Variates," Master's thesis,
% Illinois Institute of Technology, 2016.
% (describes *cubSobol_g* for control variates)
% 
% Arfon M. Smith, Daniel S. Katz, Kyle E. Niemeyer, and FORCE11 Software
% Citation Working Group, "Software citation principles," PeerJ Computer
% Science, Volume 2, e86, 2016. 
% (motivates research software citation)
% 
% Xin Tong, "A Guaranteed, Adaptive, Automatic Algorithm for Univariate
% Function Minimization," MS thesis, Illinois Institute of Technology, 2014.
% (describes deprecated *funmin01_g*)
% 
% 
%% Downloads
% 
% GAIL can be downloaded from http://gailgithub.github.io/GAIL_Dev/.
% 
% Alternatively, you can get a local copy of the GAIL repository with
% this command:
% 
%   git clone https://github.com/GailGithub/GAIL_Dev.git
% 
% 
%% Requirements
% 
% 
% You will need to install MATLAB 7 or a later version.
%
% GAIL is developed in MATLAB versions R2016a to R2020a. In particular, three
% of our core algorithms, cubSobol_g, cubBayesNet_g, and cubBayesLattice_g
% require the following MATLAB add-on toolboxes: Signal Processing Toolbox,
% Optimization Toolbox, Statistics and Machine Learning Toolbox. In MATLAB,
% we could use the following command to find out toolbox dependencies of an
% algorithm:
% 
%   names = dependencies.toolboxDependencyAnalysis({'cubBayesNet_g'}).
% 
% For development and testing purposes, we use the third-party toolboxes,
% Chebfun version 5.7.0 and Doctest for MATLAB, version 2010.
% 
%% Documentation
% 
% 
% Detailed documentation is available in the folder, GAIL_Matlab/Documentation.
%
% You can also go to MATLAB's Help. Under the section of Supplemental Software,
% you will find GAIL Toolbox's searchable HTML documentation.
% 
% A PDF version of GAIL's documentation with selected examples is available at
% https://github.com/GailGithub/GAIL_Dev/blob/master/Documentation/gail_ug_2_3_1.pdf.
% 
%% General Usage Notes
% 
% 
% GAIL Version 2.3.1 includes the following ten algorithms:
% 
% 1.  funappx_g: One-dimensional function approximation on bounded interval
% 
% 2.  funmin_g: global minimum value of univariate function on a closed interval
% 
% 3.  integral_g: One-dimensional integration on bounded interval
% 
% 4.  meanMC_g: Monte Carlo method for estimating mean of a random variable
% 
% 5.  cubMC_g: Monte Carlo method for numerical multiple integration
% 
% 6.  cubSobol_g: Quasi-Monte Carlo method using Sobol' cubature for
% d-dimensional integration
% 
% 7.  cubLattice_g: Quasi-Monte Carlo method using rank-1 Lattices cubature
% for d-dimensional integration
% 
% 8. cubBayesLattice_g: Bayesian cubature method for d-dimensional integration
%    using lattice points
% 
% 9. CubBayesNet_g: Bayesian cubature method for d-dimensional integration
%    using Sobol points
% 
% 10.  meanMC_CLT: Monte Carlo method with Central Limit Theorem (CLT)
% confidence intervals for estimating mean of a random variable
% 
% Each one of our key GAIL algorithms can parse inputs with the following three patterns of
% APIs, where f is a real-valued MATLAB function or function handle; in_param and out_param
% are MATLAB structure arrays, or in some cases, object instances; and x is an estimated
% output:
% 
% 1. Ordered input values: [x, out_param] = algo(f, inputVal1, inputVal2, inputVal3,...)
% 
% 2. Input structure array: [x, out_param] = algo(f, in_param)
% 
% 3. Ordered input values, followed by optional name-value pairs:
% 	[x, out_param] = algo(f, 'input2', inputVal2, 'input3', inputVal3,...)
%  
% 
%% Installation Instruction
% 
% 
% 1.  Unzip the contents of the zip file to a directory and maintain the
%     existing directory and subdirectory structure. (Please note: If you
%     install into the "toolbox" subdirectory of the MATLAB program
%     hierarchy, you will need to click the button "Update toolbox path
%     cache" from the File/Preferences... dialog in MATLAB.)
% 
% 2.  In MATLAB, add the GAIL directory to your path. This can be done
%     by running "GAIL_Install.m".  Alternatively, this can be done by
%     selecting "File/Set Path..." from the main or Command window
%     menus, or with the command "pathtool". We recommend that you
%     select the "Save" button on this dialog so that GAIL is on the
%     path automatically in future MATLAB sessions.
% 
% 3.  To check if you have installed GAIL successfully, type "help
%     funappx_g" to see if its documentation shows up.
% 
% Alternatively, you could do this:
% 
% 1.  Download DownloadInstallGail_2_3_1.m and put it where you want
%     GAIL to be installed.
% 
% 2.  Execute it in MATLAB.
% 
% To uninstall GAIL, execute "GAIL_Uninstall".
% 
% To reinstall GAIL, execute "GAIL_Install".
%  
%% Tests
% 
% 
% We provide quick doctests for each of the functions above. To run
% doctests in funappx_g, for example, issue the command doctest
% funappx_g.
% 
% We also provide unit tests for MATLAB version 8 or later. To run unit
% tests for funmin_g, for instance, execute run(ut_funmin_g).
% 
% We execute automated nightly fast tests and weekly long tests
% on our server. Moreover, these tests are now conducted for all MATLAB
% versions from R2016a to R2019b. The test reports are available on Mega
% cloud storage at https://mega.nz/. More specifically, fast and long test
% reports are archived at
%      https://mega.nz/#F!4olmWa6L!vYuscSnGqkvkZrGJXW5Umw
% and
%      https://mega.nz/#F!I0cAEKJD!AyQ_8tmxkknfIsuEW0_jnA
% respectively
% 
% 
%% Known Issues
% 
% During our documentation development with MATLAB release 2019a, the
% software's internal HTML viewer is found to display LaTeX expression in
% larger font size than it is intended to be. This is an aesthetic issue with
% no impact on the content accuracy.  Users may use a web browser to view our
% HTML documentation instead. The main page to GAIL's HTML documentation is
% GAIL.html, located in the subfolder, Documentation/html/.
% 
% 
%% Contact Information
% 
% 
% Please send any queries, questions, or comments to
% gail-users@googlegroups.com or visit our project website:
% http://gailgithub.github.io/GAIL_Dev/
% 
% 
%% Acknowledgements
% 
% 
% Our work was supported in part by grants from the National Science
% Foundation under grant NSF-DMS-1115392, and the Office of Advanced
% Scientific Computing Research, Office of Science, U.S. Department of
% Energy, under contract DE-AC02-06CH11357.
%
% 
%% Organizations
% 
% 
% <<illinois-institute-of-technology-vector-logo.jpg>>
% 
% <<kamakura-corporation-vector-logo.png>>
%% Release Notes 
%
% GAIL Version 2.3.1, 2020.
% 
%% Major changes in algorithms
%
% In this release, we have a new algorithm called *CubBayesNet_g*. Similar to
% cubBayesLattice_g, it is an automatic Bayesian cubature that considers 
% the integrand a realization of a Gaussian process. CubBayesNet_g uses
% Sobol points whereas cubBayesLattice_g uses lattice points.
% 
% 
%% Major changes in publications
% 
% In the folder "Papers", we have added a few recently published research
% articles and theses related to our core algorithms.
% 
% First, we have Rathinavel's PhD thesis that develops the theory
% behind cubBayesLattice_g and CubBayesNet_g.
% 
% In addition, we have included Ding, Hickernell, and Jimenez Rugama's 
% recent paper, An Adaptive Algorithm Employing Continuous Linear 
% Functionals.
 %% Functions
%
%% 1-D approximation
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
%
%% 1-D integration
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%
%% 1-D minimization
%
% <html>
% <a href="help_funmin_g.html">funmin_g</a>
% </html>
%
%
%% Higher dimensional integration
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
% 
% <html>
% <p><a href="help_cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_CLT.html">meanMC_CLT</a>
% </html>
% 


%% GAIL Demos
%
%% 1-D approximation
%
% <html>
% <a href="demo_funappx_g.html">funappx_g</a>
% </html>
%
%% 1-D minimization
%
% <html>
% <a href="demo_funmin_g.html">funmin_g</a>
% </html>
%
%% 1-D integration
%
% <html>
% <a href="demo_integral_g.html">integral_g</a>
% </html>
%
%% Higher dimensional integration
%
% <html>
% <a href="demo_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="demo_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="demo_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="demo_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="demo_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <a href="demo_cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
% <html>
% <a href="demo_meanMC_CLT.html">meanMC_CLT</a>
% </html>

%% funappx_g
% 1-D guaranteed locally adaptive function approximation (or
%   function recovery) on [a,b]
%% Syntax
% fappx = *funappx_g*(f)
%
% fappx = *funappx_g*(f,a,b,abstol)
%
% fappx = *funappx_g*(f,'a',a,'b',b,'abstol',abstol)
%
% fappx = *funappx_g*(f,in_param)
%
% [fappx, out_param] = *funappx_g*(f,...)
%% Description
%
% fappx = *funappx_g*(f) approximates function f on the default interval
%  [0,1] by an approximated function handle fappx within the guaranteed
%  absolute error tolerance of 1e-6. When Matlab version is higher or
%  equal to 8.3, fappx is an interpolant generated by griddedInterpolant.
%  When Matlab version is lower than 8.3, fappx is a function handle
%  generated by ppval and interp1. Input f is a function handle. The
%  statement y = f(x) should accept a vector argument x and return a
%  vector y of function values that is of the same size as x.
%
% fappx = *funappx_g*(f,a,b,abstol) for a given function f and the ordered
%  input parameters that define the finite interval [a,b], and a
%  guaranteed absolute error tolerance abstol.
%
% fappx = *funappx_g*(f,'a',a,'b',b,'abstol',abstol) approximates function
%  f on the finite interval [a,b], given a guaranteed absolute error
%  tolerance abstol. All four field-value pairs are optional and can be
%  supplied in different order.
%
% fappx = *funappx_g*(f,in_param) approximates function f on the finite
%  interval [in_param.a,in_param.b], given a guaranteed absolute error
%  tolerance in_param.abstol. If a field is not specified, the default
%  value is used.
%
% [fappx, out_param] = *funappx_g*(f,...) returns an approximated function
%  fappx and an output structure out_param.
%
% *Properties*
%
% * fappx can be used for linear extrapolation outside [a,b].
%
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end point of interval, default value is 0.
%
% * in_param.b --- right end point of interval, default value is 1.
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6.
%
% *Optional Input Arguments*
%
% * in_param.ninit --- initial number of subintervals. Default to 20.
%
% * in_param.nmax --- when number of points hits the value, iteration
%  will stop, default value is 1e7.
%
% * in_param.maxiter --- max number of iterations, default value is 1000.
%
% *Output Arguments*
%
% * fappx --- approximated function handle (Note: When Matlab version is
%  higher or equal to 8.3, fappx is an interpolant generated by
%  griddedInterpolant. When Matlab version is lower than 8.3, fappx is a
%  function handle generated by ppval and interp1.)
%
% * out_param.f --- input function.
%
% * out_param.a --- left end point of interval.
%
% * out_param.b --- right end point of interval.
%
% * out_param.abstol --- guaranteed absolute error tolerance.
%
% * out_param.maxiter --- max number of iterations.
%
% * out_param.ninit --- initial number of subintervals.
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- this is a vector with two elements, for
%  tracking important warnings in the algorithm. The algorithm is
%  considered successful (with out_param.exitflag == [0 0]) if no other
%  flags arise warning that the results are not guaranteed. The initial
%  value is [0 0] and the final value of this parameter is encoded as
%  follows:</li>
%   <ul type="circle">
%    <li>[1 0]:  If reaching overbudget. It states whether
%                the max budget is attained without reaching the
%                guaranteed error tolerance.</li>
%    <li>[0 1]:  If reaching overiteration. It states whether
%                the max iterations is attained without reaching the
%                guaranteed error tolerance.</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.iter --- number of iterations.
%
% * out_param.npoints --- number of points we need to reach the
%  guaranteed absolute error tolerance.
%
% * out_param.errest --- an estimation of the absolute error for the
%  approximation.
%
%% Guarantee
%
% *Please check the details of the guarantee in [1].*
%
%% Examples
% *Example 1*
%
% Approximate function $x^2$ on $[-2,2]$ with error tolerance $10^{-7}$, default
% cost budget and initial number of subintervals 18.
f = @(x) x.^2; [~, out_param] = funappx_g(f,-2,2,1e-7,18)


%%
% *Example 2*
%
% Approximate function $x^2$ on $[-2,2]$ with default error tolerance, default
% cost budget and initial number of subintervals 17.
f = @(x) x.^2;
[~, out_param] = funappx_g(f,'a',-2,'b',2,'ninit',17)


%%
% *Example 3*
%
% Approximate function $x^2$ on $[-5,5]$ with error tolerance $10^{-6}$, default
% cost budget and initial number of subintervals 18.
clear in_param; in_param.a = -5; in_param.b = 5; f = @(x) x.^2;
in_param.abstol = 10^(-6); in_param.ninit=18;
[~, out_param] = funappx_g(f,in_param)


%% See Also
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/interp1.html">interp1</a>
% </html>
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/griddedinterpolant-class.html">griddedinterpolant</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_funmin_g.html">funmin_g</a>
% </html>
%
% <http://gailgithub.github.io/GAIL_Dev/ GAIL_Dev>
%
%% References
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
% Adaption for Approximation and Minimization of Univariate Functions,"
% Journal of Complexity 40, pp. 17-33, 2017.
%
% [2] Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic
% Algorithms: Cones, Not Balls," Journal of Complexity 30, pp. 21-45,
% 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% funmin_g
% 1-D guaranteed locally adaptive function optimization on [a,b]
%% Syntax
% fmin = *funmin_g*(f)
%
% fmin = *funmin_g*(f,a,b,abstol)
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol)
%
% fmin = *funmin_g*(f,in_param)
%
% [fmin, out_param] = *funmin_g*(f,...)
%% Description
%
% fmin = *funmin_g*(f) finds minimum value of function f on the
%  default interval [0,1] within the guaranteed absolute error tolerance
%  of 1e-6. Input f is a function handle.
%
% fmin = *funmin_g*(f,a,b,abstol) finds minimum value of
%  function f with ordered input parameters that define the finite
%  interval [a,b], and a guaranteed absolute error tolerance abstol.
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol) finds minimum
%  value of function f on the interval [a,b] with a guaranteed absolute
%  error tolerance. All three field-value pairs are optional and can be
%  supplied in different order.
%
% fmin = *funmin_g*(f,in_param) finds minimum value of function f
%  on the interval [in_param.a,in_param.b] with a guaranteed absolute
%  error tolerance in_param.abstol. If a field is not specified, the
%  default value is used.
%
% [fmin, out_param] = *funmin_g*(f,...) returns minimum value fmin
%  of function f and an output structure out_param.
%
% *Input Arguments*
%
% * f --- input function.
%
% * in_param.a --- left end point of interval, default value is 0.
%
% * in_param.b --- right end point of interval, default value is 1.
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6.
%
% *Optional Input Arguments*
%
% * in_param.ninit --- initial number of subintervals. Default to 20.
%
% * in_param.nmax --- cost budget, default value is 1e7.
%
% * in_param.maxiter --- max number of iterations, default value is 1000.
%
% *Output Arguments*
%
% * out_param.f --- input function
%
% * out_param.a --- left end point of interval
%
% * out_param.b --- right end point of interval
%
% * out_param.abstol --- guaranteed absolute error tolerance
%
% * out_param.nmax --- cost budget
%
% * out_param.ninit --- initial number of points we use
%
% * out_param.npoints --- number of points needed to reach the guaranteed
%  absolute error tolerance
%
% <html>
% <ul type="square">
%  <li>out_param.exit --- this is a vector with two elements, for
%  tracking important warnings in the algorithm. The algorithm is
%  considered successful (with out_param.exit == [0 0]) if no flags arise
%  warning that the results are not guaranteed. The initial value is [0 0]
%  and the final value of this parameter is encoded as follows:</li>
%   <ul type="circle">
%    <li>[1 0]:  If reaching overbudget. It states whether
%                the max budget is attained without reaching the
%                guaranteed error tolerance.</li>
%    <li>[0 1]:  If reaching overiteration. It states whether
%                the max iterations is attained without reaching the
%                guaranteed error tolerance.</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.errest --- estimation of the absolute error bound
%
% * out_param.iter --- number of iterations
%
% * out_param.intervals --- the intervals containing point(s) where the
%  minimum occurs. Each column indicates one interval where the first raw
%  is the left point and the second row is the right point
%
%% Guarantee
%
% *Please check the details of the guarantee in [1].*
%
%% Examples
% *Example 1*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ with default input parameters.
  f = @(x) exp(0.01*(x-0.5).^2); [fmin,out_param] = funmin_g(f)


%%
% *Example 2*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-2,2]$ with error tolerance
% $10^{-7}$, cost budget $1000000$, initial number of points $10$.
f = @(x) exp(0.01*(x-0.5).^2);
[fmin,out_param] = funmin_g(f,-2,2,1e-7,10,1000000)



%%
% *Example 3*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-13,8]$ with error tolerance
% $10^{-7}$, cost budget $1000000$, initial number of points $100$.
clear in_param; in_param.a = -13; in_param.b = 8;
in_param.abstol = 1e-7;
in_param.ninit = 100;
in_param.nmax = 10^6;
[fmin,out_param] = funmin_g(f,in_param)


%%
% *Example 4*
%
% Minimize function $\exp(0.01 (x-0.5)^2)$ on $[-2,2]$ with error tolerance $10^{-5}$,
% cost budget $1000000$, initial number of points $64$.
f=@(x) exp(0.01*(x-0.5).^2);
[fmin,out_param] = funmin_g(f,'a',-2,'b',2,'ninit',64,'nmax',1e6,'abstol',1e-5)


%% See Also
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/fminbnd.html">fminbnd</a>
% </html>
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%% References
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
% Adaption for Approximation and Minimization of Univariate Functions,"
% Journal of Complexity 40, pp. 17-33, 2017.
%
% [2] Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
% Univariate Function Minimization," MS thesis, Illinois Institute of
% Technology, 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [5] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% integral_g
% 1-D guaranteed function integration using Simpson's rule
%% Syntax
% q = *integral_g*(f)
%
% q = *integral_g*(f,a,b,abstol)
%
% q = *integral_g*(f,'a',a,'b',b,'abstol',abstol)
%
% q = *integral_g*(f,in_param)
%
% [q, out_param] = *integral_g*(f,...)
%% Description
%
% q = *integral_g*(f) computes q, the definite integral of function f on
%   the interval [a,b] by Simpson's rule with in a guaranteed absolute
%   error of 1e-6. Default starting number of sample points taken is 100
%   and default cost budget is 1e7. Input f is a function handle. The
%   function y = f(x) should accept a vector argument x and return a vector
%   result y, the integrand evaluated at each element of x.
%
% q = *integral_g*(f,a,b,abstol)  computes q, the definite integral of
%   function f on the finite interval [a,b] by Simpson's rule with the
%   ordered input parameters, and guaranteed absolute error tolerance
%   abstol.
%
% q = *integral_g*(f,'a',a,'b',b,'abstol',abstol) computes q, the definite
%   integral of function f on the finite interval [a,b] by Simpson's rule
%   within a guaranteed absolute error tolerance abstol. All four
%   field-value pairs are optional and can be supplied.
%
% q = *integral_g*(f,in_param) computes q, the definite integral of
%   function f by Simpson's rule within a guaranteed absolute error
%   in_param.abstol. If a field is not specified, the default value is
%   used.
%
% [q, out_param] = *integral_g*(f,...) returns the approximated integration
%   q and output structure out_param.
%
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end of the integral, default value is 0
%
% * in_param.b --- right end of the integral, default value is 1
%
% * in_param.abstol --- guaranteed absolute error tolerance, default value
%  is 1e-6
%
% *Optional Input Arguments*
%
% * in_param.nlo --- lowest initial number of function values used, default
%  value is 10
%
% * in_param.nhi --- highest initial number of function values used,
%  default value is 1000
%
% * in_param.nmax --- cost budget (maximum number of function values),
%  default value is 1e7
%
% * in_param.maxiter --- max number of iterations, default value is 1000
%
% *Output Arguments*
%
% * q --- approximated integral
%
% * out_param.f --- input function
%
% * out_param.a --- low end of the integral
%
% * out_param.b --- high end of the integral
%
% * out_param.abstol --- guaranteed absolute error tolerance
%
% * out_param.nlo --- lowest initial number of function values
%
% * out_param.nhi --- highest initial number of function values
%
% * out_param.nmax --- cost budget (maximum number of function values)
%
% * out_param.maxiter --- max number of iterations
%
% * out_param.ninit --- initial number of points we use, computed by nlo
%  and nhi
%
% * out_param.hcut --- cut off value of the largest width between points
%  used to estimate the third derivative of the function. See [1] for
%  details.
%
% * out_param.exceedbudget --- it is true if the algorithm tries to use
%  more points than cost budget, false otherwise.
%
% * out_param.conechange --- it is true if the cone constant has been
%  changed, false otherwise. See [1] for details. If true, you may wish to
%  change the input in_param.ninit to a larger number.
%
% * out_param.npoints --- number of points we need to
%  reach the guaranteed absolute error tolerance abstol.
%
% * out_param.errest --- approximation error defined as the differences
%  between the true value and the approximated value of the integral.
%
%
%% Guarantee
%
% *Please check the details of the guarantee in [1].*
%
%
%% Examples
% *Example 1*
%
% Integrate function $x^2$ with default input parameter to make the error less
% than $10^{-6}$.
[q, out_param] = integral_g(@(x) x.^2)



%%
% *Example 2*
%
% Integrate function $\exp(-x^2)$ on $[1,2]$ with lowest initial number of function
% values $100$ and highest initial number of function values $10000$, absolute
% error tolerance $10^{-5}$ and cost budget $10000000$.
f = @(x) exp(-x.^2); [q, out_param] = integral_g(f,'a',1,'b',2,'nlo',100,'nhi',10000,...
    'abstol',1e-5,'nmax',1e7)

 
%% See Also
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/integral.html">integral</a>
% </html>
%
% <html>
% <a href="http://www.mathworks.com/help/matlab/ref/quad.html">quad</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_funmin_g.html">funmin_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Yizhi Zhang, "Guaranteed Adaptive Automatic Algorithms for 
% Univariate Integration: Methods, Costs and Implementations," PhD
% Thesis, Illinois Institute of Technology, 2018.
%
% [2] Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% meanMC_g
% Monte Carlo method to estimate the mean of a random variable.
%% Syntax
% tmu = *meanMC_g*(Yrand)
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha,
%   'fudge',fudge,'nSig',nSig,'n1',n1,'tbudget',tbudget,'nbudget',nbudget)
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param)
%% Description
%
% tmu = *meanMC_g*(Yrand) estimates the mean, mu, of a random variable Y to
%  within a specified generalized error tolerance, tolfun :=
%  max(abstol,reltol*|mu|), i.e., |mu - tmu| <= tolfun with probability at
%  least (1 - alpha), where abstol is the absolute error tolerance, and
%  reltol is the relative error tolerance. Usually the reltol determines
%  the accuracy of the estimation, however, if |mu| is rather small, then
%  abstol determines the accuracy of the estimation. Input Yrand is a
%  function handle that accepts a positive integer input n and returns an
%  n x 1 vector of IID instances of the random variable Y.
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha,fudge,nSig,n1,tbudget,nbudget)
%  estimates the mean of a random variable Y to within a specified 
%  generalized error tolerance tolfun with guaranteed confidence 
%  level 1-alpha using all ordered parsing inputs abstol, reltol, alpha, 
%  fudge, nSig, n1, tbudget, nbudget.
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha)
%  estimates the mean of a random variable Y to within a specified
%  generalized error tolerance tolfun with guaranteed confidence level
%  1-alpha. All the field-value pairs are optional and can be supplied in
%  different order, if a field is not supplied, the default value is used.
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param) estimates the mean of a
%  random variable Y to within a specified generalized error tolerance
%  tolfun with the given parameters in_param and produce the estimated
%  mean tmu and output parameters out_param. If a field is not specified,
%  the default value is used.
%
% *Input Arguments*
%
% * Yrand --- he function for generating n IID instances of a random
%  variable Y whose mean we want to estimate. Y is often defined as a
%  function of some random variable X with a simple distribution. The
%  input of Yrand should be the number of random variables n, the output
%  of Yrand should be n function values. For example, if Y = X.^2 where X
%  is a standard uniform random variable, then one may define Yrand =
%  @(n) rand(n,1).^2.
%
% * in_param.abstol --- the absolute error tolerance, which should be
%  positive, default value is 1e-2.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  between 0 and 1, default value is 1e-1.
%
% * in_param.alpha --- the uncertainty, which should be a small positive
%  percentage, default value is 1%.
%
% * in_param.fudge --- standard deviation inflation factor, which should
%  be larger than 1, default value is 1.2.
%
% * in_param.nSig --- initial sample size for estimating the sample
%  variance, which should be a moderately large integer bigger than or
%  equal to 30, the default value is 1e4.
%
% * in_param.n1 --- initial sample size for estimating the sample mean,
%  which should be a moderate large positive integer at least 30, the
%  default value is 1e4.
%
% * in_param.tbudget --- the time budget in seconds to do the two-stage
%  estimation, which should be positive, the default value is 100 seconds.
%
% * in_param.nbudget --- the sample budget to do the two-stage
%  estimation, which should be a large positive integer, the default
%  value is 1e9.
%
% *Output Arguments*
%
% * tmu --- the estimated mean of Y.
%
% * out_param.tau --- the iteration step.
%
% * out_param.n --- the sample size used in each iteration.
%
% * out_param.nremain --- the remaining sample budget to estimate mu. It was
%  calculated by the sample left and time left.
%
% * out_param.ntot --- total sample used.
%
% * out_param.hmu --- estimated mean in each iteration.
%
% * out_param.tol --- the reliable upper bound on error for each iteration.
%
% * out_param.var --- the sample variance.
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0   Successs</li>
%    <li>1   Not enough samples to estimate the mean</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.kurtmax --- the upper bound on modified kurtosis.
%
% * out_param.time --- the time elapsed in seconds.
%
%
%%  Guarantee
% This algorithm attempts to calculate the mean, mu, of a random variable
% to a prescribed error tolerance, tolfun:= max(abstol,reltol*|mu|), with
% guaranteed confidence level 1-alpha. If the algorithm terminates without
% showing any warning messages and provides an answer tmu, then the follow
% inequality would be satisfied:
%
% Pr(|mu-tmu| <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% defined in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta. And
% the following inequality holds:
%
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
%% Examples
%
%%
% *Example 1*
%
% Calculate the mean of $x^2$ when $x$ is uniformly distributed in
% $[0, 1]$, with the absolute error tolerance = $10^{-3}$ and uncertainty $5\%$.

  in_param.reltol = 0; in_param.abstol = 1e-3;
  in_param.alpha = 0.05; Yrand = @(n) rand(n,1).^2;
  tmu = meanMC_g(Yrand,in_param); exactsol = 1/3;
  check = double(abs(exactsol-tmu) < 1e-3)

%%
% *Example 2*
%
% Calculate the mean of $\exp(x)$ when $x$ is uniformly distributed in
% $[0, 1]$, with the absolute error tolerance $10^{-3}$.

  tmu = meanMC_g(@(n)exp(rand(n,1)),1e-3,0); exactsol = exp(1)-1;
  check = double(abs(exactsol-tmu) < 1e-3)

%%
% *Example 3*
%
% Calculate the mean of $\cos(x)$ when $x$ is uniformly distributed in
% $[0, 1]$, with the relative error tolerance $10^{-2}$ and uncertainty $0.05$.

  tmu = meanMC_g(@(n)cos(rand(n,1)),'reltol',1e-3,'abstol',1e-4,'alpha',0.01);
  exactsol = sin(1);
  check = double(abs(exactsol-tmu) < max(1e-3,1e-2*abs(exactsol)))

%% See Also
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, "Guaranteed
% conservative fixed width confidence intervals via Monte Carlo
% sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
% Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), pp. 105-128,
% Springer-Verlag, Berlin, 2014. DOI: 10.1007/978-3-642-41095-6_5
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
% Statistics. London, UK: CHAPMAN & HALL
%
% [7] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
% Means of Random Variables, PhD Thesis, Illinois Institute of
% Technology, 2016.
%
% If you find GAIL helpful in your work, please support us by citing the
% above paper and software.
%% meanMC_CLT
% Monte Carlo method to estimate the mean of a random variable
%% Syntax
% sol = *MEANMC_CLT*(Y,absTol,relTol,alpha,nSig,inflate)
%
%% Description
%
% sol = *MEANMC_CLT*(Y,absTol,relTol,alpha,nSig,inflate) estimates the
%  mean, mu, of a random variable to within a specified error tolerance,
%  i.e., | mu - tmu | <= max(absTol,relTol|mu|) with probability at least
%  1-alpha, where abstol is the absolute error tolerance.  The default
%  values are abstol=1e-2 and alpha=1%. Input Y is a function handle that
%  accepts a positive integer input n and returns an n x 1 vector of IID
%  instances of the random variable.
%
% This is a heuristic algorithm based on a Central Limit Theorem
%  approximation.
%   
% *Input Arguments*
%
% * Y --- the function or structure for generating n IID instances of a
%     random variable Y whose mean we want to estimate. Y is often defined
%     as a function of some random variable X with a simple distribution.
%     The input of Yrand should be the number of random variables n, the
%     output of Yrand should be n function values. For example, if Y = X.^2
%     where X is a standard uniform random variable, then one may define
%     Yrand = @(n) rand(n,1).^2.
%
% * absTol --- the absolute error tolerance, which should be
%  non-negative --- default = 1e-2
%
% * relTol --- the relative error tolerance, which should be
%  non-negative and no greater than 1 --- default = 0
%
% * alpha --- the uncertainty, which should be a small positive
%  percentage --- default = 1%
%
% * nSig --- the number of samples used to compute the sample variance
%  --- default = 1000
%
% * inflate --- the standard deviation inflation factor --- default = 1.2
%
% *Output Arguments*
%
% * Y --- the random generator
%
% * absTol --- the absolute error tolerance
%
% * relTol --- the relative error tolerance
%
% * alpha --- the uncertainty
%
% * mu --- the estimated mean of Y.
%
% * stddev --- sample standard deviation of the random variable
%
% * nSample --- total sample used.
%
% * time --- the time elapsed in seconds.
%
% * errBd --- the error bound.
%
%% Examples
%
%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = x_1 x_2$ in the interval
% $[0,1]^2$ with absolute tolerance $10^{-3}$ and relative tolerance $0$:

  [mu,out] = meanMC_CLT(@(n) prod(rand(n,2),2), 0.001);
  exact = 1/4;
  check = double(abs(exact - mu) < 2e-3)
  
%%
% *Example 2*
%
% Estimate the integral $f(x)=\exp(-x^2)$ in the interval $[0,1]$
% using $x$ as a control variate and relative error $10^{-3}$:

  f = @(x)[exp(-x.^2), x];
  YXn = @(n)f(rand(n,1));
  s = struct('Y',YXn,'nY',1,'trueMuCV',1/2);
  exact = erf(1)*sqrt(pi)/2;
  success = 0; runs = 1000; tol = 1e-3;
  for i=1:runs, success = success + double(abs(exact-meanMC_CLT(s,0,tol)) < tol*exact); end
  check = success >= 0.99 * runs  
  
%%
% *Example 3*
%
% Estimate the Keister's integration in dimension 1 with $a=1$, $\frac{1}{\sqrt 2}$ and
% using $\cos(x)$ as a control variate:

  normsqd = @(x) sum(x.*x,2);
  f = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)).* exp((1/2-a^2)*normt);
  f1 = @(x,a,d) f(normsqd(x),a,d);
  f2 = @(x)[f1(x,1,1),f1(x,1/sqrt(2),1),cos(x)];
  YXn = @(n)f2(randn(n,1));
  s = struct('Y',YXn,'nY',2,'trueMuCV',1/sqrt(exp(1)));
  [hmu,out] = meanMC_CLT(s,0,1e-3);
  exact = 1.380388447043143;
  check = double(abs(exact-hmu) < max(0,1e-3*abs(exact)))

%%
% *Example 4*
%
% Estimate the integral with integrand $f(x) = x_1^3 x_2^3 x_3^3$ in the
% interval $[0,1]^3$ with pure absolute error $10^{-3}$ using $x_1 x_2 x_3$ as
% a control variate:

  f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, x(:,1).*x(:,2).*x(:,3)];
  s = struct('Y',@(n)f(rand(n,3)),'nY',1,'trueMuCV',1/8);
  [hmu,out] = meanMC_CLT(s,1e-3,0);
  exact = 1/64;
  check = double(abs(exact-hmu) < max(1e-3,1e-3*abs(exact)))

%%
% *Example 5*
%
% Estimate the integrals with integrands $f_1(x) = x_1^3 x_2^3 x_3^3$ and 
% $f_2(x)= x_1^2 x_2^2 x_3^2- \frac{1}{27}+\frac{1}{64}$ in the interval $[0,1]^3$
% using $x_1 x_2 x_3$ and $x_1+x_2+x_3$ as control variates:

  f = @(x) [x(:,1).^3.*x(:,2).^3.*x(:,3).^3, ...
            x(:,1).^2.*x(:,2).^2.*x(:,3).^2-1/27+1/64, ...
            x(:,1).*x(:,2).*x(:,3), ...
            x(:,1)+x(:,2)+x(:,3)];
  s = struct('Y',@(n)f(rand(n,3)),'nY',2,'trueMuCV',[1/8 1.5]);
  [hmu,out] = meanMC_CLT(s,1e-4,1e-3);
  exact = 1/64;
  check = double(abs(exact-hmu) < max(1e-4,1e-3*abs(exact)))
  
%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.%% cubMC_g
% Monte Carlo method to evaluate a multidimensional integral
%% Syntax
% [Q,out_param] = *cubMC_g*(f,hyperbox)
%
% Q = *cubMC_g*(f,hyperbox,measure,abstol,reltol,alpha)
%
% Q = *cubMC_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol,'alpha',alpha)
%
% [Q out_param] = *cubMC_g*(f,hyperbox,in_param)
%% Description
%
% [Q,out_param] = *cubMC_g*(f,hyperbox) estimates the integral of f over
%  hyperbox to within a specified generalized error tolerance, tolfun =
%  max(abstol, reltol*| I |), i.e., | I - Q | <= tolfun with probability
%  at least (1 - alpha), where abstol is the absolute error tolerance, and
%  reltol is the relative error tolerance. Usually the reltol determines
%  the accuracy of the estimation, however, if | I | is rather small,
%  then abstol determines the accuracy of the estimation. Input f is a
%  function handle that accepts an n x d matrix input, where d is the
%  dimension of the hyperbox, and n is the number of points being
%  evaluated simultaneously. 
%
% When measure is 'uniform', 'uniform box', 'normal' or 'Gaussian', the
%  input hyperbox is a 2 x d matrix, where the first row corresponds to the
%  lower limits and the second row corresponds to the upper limits. When
%  measure is 'uniform ball' or 'uniform sphere', the input hyperbox is a
%  vector with d+1 elements, where the first d values correspond to the
%  center of the ball and the last value corresponds to the radius of the
%  ball. For these last two measures, a user can optionally specify what
%  transformation should be used in order to get a uniform distribution on
%  a ball of sphere. When measure is 'uniform ball_box', the box-to-ball
%  transformation, which gets a set of points uniformly distributed on a
%  ball from a set of points uniformly distributed on a box, will be used.
%  When measure is 'uniform ball_normal', the normal-to-ball
%  transformation, which gets a set of points uniformly distributed on a
%  ball from a set of points normally distributed on the space, will be
%  used. Similarly, the measures 'uniform sphere_box' and 'uniform
%  sphere_normal' can be defined. The default transformations are the
%  box-to-ball and the box-to-sphere transformations, depending on the
%  region of integration.
%
% Q = *cubMC_g*(f,hyperbox,measure,abstol,reltol,alpha)
%  estimates the integral of function f over hyperbox to within a
%  specified generalized error tolerance tolfun with guaranteed confidence
%  level 1-alpha using all ordered parsing inputs f, hyperbox, measure,
%  abstol, reltol, alpha, fudge, nSig, n1, tbudget, nbudget, flag. The
%  input f and hyperbox are required and others are optional.
%
% Q = *cubMC_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol,'alpha',alpha)
%  estimates the integral of f over hyperbox to within a specified
%  generalized error tolerance tolfun with guaranteed confidence level
%  1-alpha. All the field-value pairs are optional and can be supplied in
%  different order. If an input is not specified, the default value is used.
%
% [Q out_param] = *cubMC_g*(f,hyperbox,in_param) estimates the integral of
%  f over hyperbox to within a specified generalized error tolerance
%  tolfun with the given parameters in_param and produce output parameters
%  out_param and the integral Q.
%
% *Input Arguments*
%
% * f --- the integrand.
%
% * hyperbox --- the integration hyperbox. The default value is
%  [zeros(1,d); ones(1,d)], the default d is 1.
%
% * in_param.measure --- the measure for generating the random variable,
%  the default is 'uniform'. The other measures could be handled are
%  'uniform box', 'normal'/'Gaussian', 'uniform ball'/'uniform
%  ball_box'/'uniform ball_normal' and 'uniform sphere'/'uniform
%  sphere_box'/'uniform sphere_normal'. The input should be
%  a string type, hence with quotes.
%
% * in_param.abstol --- the absolute error tolerance, the default value
%  is 1e-2.
%
% * in_param.reltol --- the relative error tolerance, the default value
%  is 1e-1.
%
% * in_param.alpha --- the uncertainty, the default value is 1%.
%
% *Optional Input Arguments*
%
% * in_param.fudge --- the standard deviation inflation factor, the
%  default value is 1.2.
%
% * in_param.nSig --- initial sample size for estimating the sample
%  variance, which should be a moderate large integer at least 30, the
%  default value is 1e4.
%
% * in_param.n1 --- initial sample size for estimating the sample mean,
%  which should be a moderate large positive integer at least 30, the
%  default value is 1e4.
%
% * in_param.tbudget --- the time budget to do the estimation, the
%  default value is 100 seconds.
%
% * in_param.nbudget --- the sample budget to do the estimation, the
%  default value is 1e9.
%
% <html>
% <ul type="square">
%  <li>in_param.flag --- the value corresponds to parameter checking status:</li>
%   <ul type="circle">
%    <li>0   not checked</li>
%    <li>1   checked by meanMC_g</li>
%    <li>2   checked by cubMC_g</li>
%   </ul>
% </ul>
% </html>
%
% *Output Arguments*
%
% * Q --- the estimated value of the integral.
%
% * out_param.n --- the sample size used in each iteration.
%
% * out_param.ntot --- total sample used, including the sample used to
%  convert time budget to sample budget and the sample in each iteration
%  step.
%
% * out_param.nremain --- the remaining sample budget to estimate I. It was
%  calculated by the sample left and time left.
%
% * out_param.tau --- the iteration step.
%
% * out_param.hmu --- estimated integral in each iteration.
%
% * out_param.tol --- the reliable upper bound on error for each iteration.
%
% * out_param.kurtmax --- the upper bound on modified kurtosis.
%
% * out_param.time --- the time elapsed in seconds.
%
% * out_param.var --- the sample variance.
%
% <html>
% <ul type="square">
%  <li>out_param.exit --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0   success</li>
%    <li>1   Not enough samples to estimate the mean</li>
%    <li>10  hyperbox does not contain numbers</li>
%    <li>11  hyperbox is not 2 x d when measure is 'uniform'
%                           or 'normal'</li>
%    <li>12  hyperbox is only a point in one direction when
%                           measure is 'uniform' or 'normal'</li>
%    <li>13  hyperbox is infinite when measure is 'uniform'</li>
%    <li>14  hyperbox is not doubly infinite when measure
%                           is 'normal'</li>
%    <li>15  hyperbox has an infinite coordinate for the
%            center of the ball or sphere or a infinite radius
%                           for the ball or sphere</li>
%    <li>16  The radius of the ball or sphere is a non-positive
%                           real number</li>
%    <li>18  Hyperbox not 1 x (d+1) when measure is 'uniform
%                           ball' or 'uniform sphere'</li>
%    <li>19  The dimension of the sphere is smaller than 2</li>
%   </ul>
% </ul>
% </html>
%
%%  Guarantee
% This algorithm attempts to calculate the integral of function f over a
% hyperbox to a prescribed error tolerance |tolfun| = max(|abstol|, |reltol| ||I||)
% with guaranteed confidence level |1-alpha|. If the algorithm terminates
% without showing any warning messages and provides an answer |Q|, then the
% following inequality would be satisfied:
%
% |Pr(| Q - I | <= tolfun) >= 1-alpha|.
%
% The cost of the algorithm, |N_tot|, is also bounded above by |N_up|, which is
% a function in terms of |abstol|, |reltol|, |nSig|, |n1|, |fudge|, |kurtmax|, |beta|. And
% the following inequality holds:
%
% |Pr (N_tot <= N_up) >= 1-beta|.
%
% Please refer to our paper for detailed arguments and proofs.
%
%%  Examples

%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = \sin(x)$ over the interval
% $[1,2]$ with default parameters.

 f = @(x) sin(x); interval = [1;2];
 Q = cubMC_g(f,interval,'uniform',1e-3,1e-2);
 exactsol = 0.9564;
 check = double(abs(exactsol-Q) < max(1e-3,1e-2*abs(exactsol)))

%%
% *Example 2*
%
% Estimate the integral with integrand $f(x) = \exp(-x_1^2-x_2^2)$ over the
% hyperbox $[0, 0; 1, 1]$, where $x= [x_1, x_2]$ is a vector.

 f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [0 0;1 1];
 Q = cubMC_g(f,hyperbox,'uniform',1e-3,0);
 exactsol = 0.5577;
 check = double(abs(exactsol-Q) < 1e-3)

%%
% *Example 3*
%
% Estimate the integral with integrand $f(x) = 2^d \prod(x_1 x_2 \cdots
% x_d)+0.555$ over the hyperbox |[zeros(1,d); ones(1,d)]|, where $x =
% [x_1, x_2, \ldots, x_d]$ is a vector.

  d = 3; f = @(x) 2^d*prod(x,2)+0.555; hyperbox = [zeros(1,d); ones(1,d)];
  in_param.abstol = 1e-3;in_param.reltol=1e-3;
  Q = cubMC_g(f,hyperbox,in_param);
  exactsol = 1.555;
  check = double(abs(exactsol-Q) < max(1e-3,1e-3*abs(exactsol)))

%%
% *Example 4*
%
% Estimate the integral with integrand $f(x) = \exp(-x_1^2-x_2^2)$ in
% $R^2$, where $x = [x_1, x_2]$ is a vector.

 f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-inf -inf;inf inf];
 Q = cubMC_g(f,hyperbox,'normal',0,1e-2);
 exactsol = 1/3;
 check = double(abs(exactsol-Q) < max(0,1e-2*abs(exactsol)))

%%
% *Example 5*
%
% Estimate the integral with integrand $f(x) = x_1^2+x_2^2$ in the disk with
% center $(0,0)$ and radius 1, where $x = [x_1, x_2]$ is a vector.

 f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1];
 Q = cubMC_g(f,hyperbox,'uniform ball','abstol',1e-3,'reltol',1e-3);
 exactsol = pi/2;
 check = double(abs(exactsol-Q) < max(1e-3,1e-3*abs(exactsol)))

%% See Also
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen, "Guaranteed
% conservative fixed width confidence intervals via Monte Carlo
% sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
% Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), pp. 105-128,
% Springer-Verlag, Berlin, 2014. DOI: 10.1007/978-3-642-41095-6_5
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
% Statistics. London, UK: CHAPMAN & HALL
%
% [7] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
% Means of Random Variables, PhD Thesis, Illinois Institute of
% Technology, 2016.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% cubLattice_g
% Quasi-Monte Carlo method using rank-1 Lattices cubature
% over a d-dimensional region to integrate within a specified generalized
% error tolerance with guarantees under Fourier coefficients cone decay
% assumptions.
%% Syntax
% [q,out_param] = *cubLattice_g*(f,hyperbox)
%
% q = *cubLattice_g*(f,hyperbox,measure,abstol,reltol)
%
% q = *cubLattice_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%
% q = *cubLattice_g*(f,hyperbox,in_param)
%% Description
%
% [q,out_param] = *cubLattice_g*(f,hyperbox) estimates the integral of f
%  over the d-dimensional region described by hyperbox, and with an error
%  guaranteed not to be greater than a specific generalized error
%  tolerance, tolfun:=max(abstol,reltol*| integral(f) |). Input f is a
%  function handle. f should accept an n x d matrix input, where d is the
%  dimension and n is the number of points being evaluated simultaneously.
%
% When measure is 'uniform', the input hyperbox is a 2 x d matrix, where
%  the first row corresponds to the lower limits and the second row
%  corresponds to the upper limits of the integral. When measure is
%  'uniform ball' or 'uniform sphere', the input hyperbox is a vector with
%  d+1 elements, where the first d values correspond to the center of the
%  ball and the last value corresponds to the radius of the ball. For these
%  last two measures, a user can optionally specify what transformation
%  should be used in order to get a uniform distribution on a ball. When
%  measure is 'uniform ball_box', the box-to-ball transformation, which
%  gets a set of points uniformly distributed on a ball from a set of
%  points uniformly distributed on a box, will be used. When measure is
%  'uniform ball_normal', the normal-to-ball transformation, which gets a
%  set of points uniformly distributed on a ball from a set of points
%  normally distributed on the space, will be used. Similarly, the measures
%  'uniform sphere_box' and 'uniform sphere_normal' can be used to specify
%  the desired transformations. The default transformations are the
%  box-to-ball and the box-to-sphere transformations, depending on the
%  region of integration. Given the construction of our Lattices, d must be
%  a positive integer with 1 <= d <= 600.
% 
% q = *cubLattice_g*(f,hyperbox,measure,abstol,reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either. Inputs f and hyperbox 
%  are required. The other optional inputs are in the correct order:
%  measure,abstol,reltol,shift,mmin,mmax,fudge, and transform.
%
% 
% q = *cubLattice_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All the field-value
%  pairs are optional and can be supplied in any order. If an input is not
%  specified, the default value is used.
% 
% q = *cubLattice_g*(f,hyperbox,in_param) estimates the integral of f over the
%  hyperbox. The answer is given within the generalized error tolerance tolfun.
% 
% *Input Arguments*
%
% * f --- the integrand whose input should be a matrix n x d where n is
%  the number of data points and d the dimension, which cannot be
%  greater than 600. By default f is f=@ x.^2.
%
% * hyperbox --- the integration region defined by its bounds. When
%  measure is 'uniform' or 'normal', hyperbox must be a 2 x d matrix,
%  where the first row corresponds to the lower limits and the second
%  row corresponds to the upper limits of the integral. When measure is
%  'uniform ball' or 'uniform sphere', the input hyperbox is a vector
%  with d+1 elements, where the first d values correspond to the center
%  of the ball and the last value corresponds to the radius of the ball.
%  The default value is [0;1].
%
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in the hyperbox or
%  normally distributed with covariance matrix I_d. The possible values
%  are 'uniform', 'normal', 'uniform ball', 'uniform ball_box', 'uniform
%  ball_normal', 'uniform sphere', 'uniform sphere_box' and 'uniform
%  sphere_normal'. For 'uniform', the hyperbox must be a finite volume,
%  for 'normal', the hyperbox can only be defined as (-Inf,Inf)^d and,
%  for 'uniform ball' or 'uniform sphere', hyperbox must have finite
%  values for the coordinates of the center and a finite positive value
%  for the radius. By default it is 'uniform'.
%
% * in_param.abstol --- the absolute error tolerance, abstol>=0. By 
%  default it is 1e-4. For pure absolute tolerance, set in_param.reltol
%  = 0.
%
% * in_param.reltol --- the relative error tolerance, which should be in
%  [0,1]. Default value is 1e-2. For pure absolute tolerance, set
%  in_param.abstol = 0.
% 
% *Optional Input Arguments*
%
% * in_param.shift --- in_param.shift --- the Rank-1 lattices can be shifted to avoid the
%  origin or other particular points. The shift is a vector in [0,1]^d.
%  By default we consider a shift uniformly sampled from [0,1]^d.
% 
% * in_param.mmin --- the minimum number of points to start is 2^mmin.
%  The cone condition on the Fourier coefficients decay requires a
%  minimum number of points to start. The advice is to consider at least
%  mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%  default it is 10.
% 
% * in_param.mmax --- tthe maximum budget is 2^mmax. By construction of
%  our Lattices generator, mmax is a positive integer such that
%  mmin<=mmax. mmax should not be bigger than the gail.lattice_gen
%  allows. The default value is 20.
% 
% * in_param.fudge --- the positive function multiplying the finite 
%  sum of Fast Fourier coefficients specified in the cone of functions.
%  This input is a function handle. The fudge should accept an array of
%  nonnegative integers being evaluated simultaneously. For more
%  technical information about this parameter, refer to the references.
%  By default it is @(m) 5*2.^-m.
% 
% <html>
% <ul type="square">
%  <li>in_param.transform --- the algorithm is defined for continuous
%     periodic functions. If the input function f is not, there are 5
%     types of transform to periodize it without modifying the result. 
%     By default it is the Baker's transform. The options are:</li>
%   <ul type="circle">
%    <li>id : no transformation.</li>
%    <li>Baker : Baker's transform or tent map in each coordinate. Preserving
%            only continuity but simple to compute. Chosen by default.</li>
%    <li>C0 : polynomial transformation only preserving continuity.</li>
%    <li>C1 : polynomial transformation preserving the first derivative.</li>
%    <li>C1sin : Sidi's transform with sine, preserving the first derivative.
%              This is in general a better option than 'C1'.</li>
%   </ul>
%  </ul>
% </html>
%
% *Output Arguments*
%
% * q --- the estimated value of the integral.
%
% * out_param.d --- dimension over which the algorithm integrated.
% 
% * out_param.n --- number of Rank-1 lattice points used for computing
%  the integral of f.
% 
% * out_param.bound_err --- predicted bound on the error based on the
%  cone condition. If the function lies in the cone, the real error will
%  be smaller than generalized tolerance.
% 
% * out_param.time --- time elapsed in seconds when calling cubLattice_g.
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- this is a binary vector stating whether
%  warning flags arise. These flags tell about which conditions make the
%  final result certainly not guaranteed. One flag is considered arisen
%  when its value is 1. The following list explains the flags in the
%  respective vector order:</li>
%   <ul type="circle">
%                    <li>1 : If reached overbudget, meaning the max budget
%                     is attained without reaching the guaranteed error
%                     tolerance.</li> 
%                    <li>2 : If the function lies outside the cone, results
%                     are not guaranteed to be accurate. Note that this
%                     parameter is computed on the transformed function,
%                     not the input function. For more information on the
%                     transforms, check the input parameter
%                     in_param.transform; for information about the cone
%                     definition, check the article mentioned below.</li>
%   </ul>
%  </ul>
% </html>
%
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in $[0,1]^d$
% with a prescribed generalized error tolerance. The Fourier coefficients
% of the integrand are assumed to be absolutely convergent. If the
% algorithm terminates without warning messages, the output is given with
% guarantees under the assumption that the integrand lies inside a cone of
% functions. The guarantee is based on the decay rate of the Fourier
% coefficients. For integration over domains other than $[0,1]^d$, this cone
% condition applies to $f \circ \psi$ (the composition of the
% functions) where $\psi$ is the transformation function for $[0,1]^d$ to
% the desired region. For more details on how the cone is defined, please
% refer to the references below.
%
%% Examples
%
%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = x_1 x_2$ in the interval $[0,1]^2$:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)]; 
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','C1sin'); 
  exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)
  
%%
% *Example 2*
%
% Estimate the integral with integrand $f(x) = x_1^2  x_2^2 x_3^2$ in the
% interval $R^3$ where $x_1$, $x_2$ and $x_3$ are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
  q = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,...
      'transform','C1sin','shift',2^(-25)*ones(1,3)); 
  exactsol = 1;
  check = double(abs(exactsol-q) < max(1e-3,1e-3*abs(exactsol)))
  
%%
% *Example 3*
%
% Estimate the integral with integrand $f(x) = \exp(-x_1^2-x_2^2)$ in the
% interval $[-1,2]^2$:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
  q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1'); 
  exactsol = 1/9*(sqrt(pi)/2*(erf(2)+erf(1)))^2;
  check = double(abs(exactsol-q) < max(1e-3,1e-2*abs(exactsol)))

%%
% *Example 4*
%
% Estimate the price of an European call with $S_0=100$, $K=100$, $r=\sigma^2/2$,
% $\sigma=0.05$, and $T=1$.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); 
  hyperbox = [-inf(1,1);inf(1,1)];
  q = cubLattice_g(f,hyperbox,'normal',1e-4,1e-2,'transform','C1sin'); 
  price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
  check = double(abs(price-q) < max(1e-4,1e-2*abs(price)))

%%
% *Example 5*
%
% Estimate the integral with integrand $f(x) = 8 x_1 x_2 x_3 x_4 x_5$ in the
% interval $[0,1]^5$ with pure absolute error $10^{-5}$.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 6*
%
% Estimate the integral with integrand $f(x) = 3/(5-4 (\cos(2 \pi x)))$ in the
% interval $[0,1]$ with pure absolute error $10^{-5}$.

  f = @(x) 3./(5-4*(cos(2*pi*x))); hyperbox = [0;1];
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','id'); 
  exactsol = 1;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 7*
%
% Estimate the integral with integrand $f(x) = x_1^2+x_2^2$ over the disk
% with center $(0,0)$ and radius 1 with pure absolute error $10^{-4}$,
% where  $x = [x_1 x_2]$ is a vector.

  f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1];
  q = cubLattice_g(f,hyperbox,'uniform ball','abstol',1e-4,'reltol',0); 
  exactsol = pi/2;
  check = double(abs(exactsol-q) < 1e-4)
  
%% See Also
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive
% multidimensional integration based on rank-1 lattices," Monte Carlo and
% Quasi-Monte Carlo  Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools
% and D. Nuyens, eds.), Springer Proceedings in Mathematics and Statistics,
% vol. 163, Springer-Verlag, Berlin, 2016, arXiv:1411.1966, pp. 407-422.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% [6] Kai-Tai Fang and Yuan Wang, Number-theoretic Methods in 
% Statistics, Chapman & Hall, London, 1994.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.

%% cubSobol_g
% Quasi-Monte Carlo method using Sobol' cubature over the
% d-dimensional region to integrate within a specified generalized error
% tolerance with guarantees under Walsh-Fourier coefficients cone decay
% assumptions
%% Syntax
% [q,out_param] = *cubSobol_g*(f,hyperbox)
%
% q = *cubSobol_g*(f,hyperbox,measure,abstol,reltol)
%
% q = *cubSobol_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%
% q = *cubSobol_g*(f,hyperbox,in_param)
%% Description
%
% [q,out_param] = *cubSobol_g*(f,hyperbox) estimates the integral of f
%  over the d-dimensional region described by hyperbox, and with an error
%  guaranteed not to be greater than a specific generalized error tolerance,
%  tolfun:=max(abstol,reltol*| integral(f) |). Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension and n is the
%  number of points being evaluated simultaneously. 
%
% When measure is 'uniform', the input hyperbox is a 2 x d matrix, where
% the first row corresponds to the lower limits and the second row
% corresponds to the upper limits of the integral. When measure is 'uniform
% ball' or 'uniform sphere', the input hyperbox is a vector with d+1
% elements, where the first d values correspond to the center of the ball
% and the last value corresponds to the radius of the ball. For these last
% two measures, a user can optionally specify what transformation should be
% used in order to get a uniform distribution on a ball. When measure is
% 'uniform ball_box', the box-to-ball transformation, which gets a set of
% points uniformly distributed on a ball from a set of points uniformly
% distributed on a box, will be used. When measure is 'uniform
% ball_normal', the normal-to-ball transformation, which gets a set of
% points uniformly distributed on a ball from a set of points normally
% distributed on the space, will be used. Similarly, the measures 'uniform
% sphere_box' and 'uniform sphere_normal' can be used to specify the
% desired transformations. The default transformations are the box-to-ball
% and the box-to-sphere transformations, depending on the region of
% integration. Given the construction of Sobol' sequences, d must be a
% positive integer with 1 <= d<= 1111.
%
% q = *cubSobol_g*(f,hyperbox,measure,abstol,reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either. Inputs f and hyperbox
%  are required. The other optional inputs are in the correct order:
%  measure,abstol,reltol,mmin,mmax,and fudge.
%
% q = *cubSobol_g*(f,hyperbox,'measure',measure,'abstol',abstol,'reltol',reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All the field-value
%  pairs are optional and can be supplied in any order. If an input is not
%  specified, the default value is used.
%
% q = *cubSobol_g*(f,hyperbox,in_param) estimates the integral of f over the
%  hyperbox. The answer is given within the generalized error tolerance tolfun.
%
% *Input Arguments*
%
% <html><ul type=square>
%  <li>f --- the integrand whose input should be a matrix n x d where n is
%  the number of data points and d the dimension, which cannot be
%  greater than 1111. By default f is f=@ x.^2.</li>
%  --- if using control variates, f needs to be a structure with two fields:
%  First field: 'func', need to be a function handle with n x (J+1)
%  dimension outputs, where J is the number of control variates.
%  First column is the output of target function, next J columns are
%  the outputs of control variates.
%  Second field: 'cv', need to be a 1 x J vector that stores the
%  exact means of control variates in the same order from
%  the function handle. For examples of how to use control variates,
%  please check Example 7 below.</li>
%  </ul>
% </html>
%
% * hyperbox --- the integration region defined by its bounds. When measure
%  is 'uniform' or 'normal', hyperbox must be a 2 x d matrix, where the
%  first row corresponds to the lower limits and the second row corresponds
%  to the upper limits of the integral. When measure is 'uniform ball'
%  or 'uniform sphere', the input hyperbox is a vector with d+1 elements,
%  where the first d values correspond to the center of the ball and the
%  last value corresponds to the radius of the ball. The default value
%  is [0;1].
%
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in the hyperbox
%  or normally distributed with covariance matrix I_d. The possible
%  values are 'uniform', 'normal', 'uniform ball', 'uniform ball_box',
%  'uniform ball_normal', 'uniform sphere', 'uniform sphere_box' and
%  'uniform sphere_normal'. For 'uniform', the hyperbox must be a
%  finite volume, for 'normal', the hyperbox can only be defined as
%  (-Inf,Inf)^d and, for 'uniform ball' or 'uniform sphere', hyperbox
%  must have finite values for the coordinates of the center and a
%  finite positive value for the radius. By default it is 'uniform'.
%
% * in_param.abstol --- the absolute error tolerance, abstol>=0. By
%  default it is 1e-4. For pure absolute tolerance, set in_param.reltol
%  = 0.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  in [0,1]. Default value is 1e-2. For pure absolute tolerance, set
%  in_param.abstol = 0.
%
% *Optional Input Arguments*
%
% * in_param.mmin --- the minimum number of points to start is 2^mmin.
%  The cone condition on the Fourier coefficients decay requires a
%  minimum number of points to start. The advice is to consider at least
%  mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%  default it is 10.
%
% * in_param.mmax --- the maximum budget is 2^mmax. By construction of
%  the Sobol' generator, mmax is a positive integer such that
%  mmin<=mmax<=53. The default value is 24.
%
% * in_param.fudge --- the positive function multiplying the finite
%  sum of Fast Walsh Fourier coefficients specified in the cone of functions.
%  This input is a function handle. The fudge should accept an array of
%  nonnegative integers being evaluated simultaneously. For more
%  technical information about this parameter, refer to the references.
%  By default it is @(m) 5*2.^-m.
%
% *Output Arguments*
%
% * q --- the estimated value of the integral.
%
% * out_param.d --- dimension over which the algorithm integrated.
%
% * out_param.n --- number of Sobol' points used for computing the
%  integral of f.
%
% * out_param.bound_err --- predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error will be
%  smaller than generalized tolerance.
%
% * out_param.time --- time elapsed in seconds when calling cubSobol_g.
%
% * out_param.beta --- the value of beta when using control variates
%     as in f-(h-Ih)beta, if using 'betaUpdate' option, beta is a vector
%     storing value of each iteration.
%
% * y --- fast transform coefficients of the input function.
%
% * kappanumap --- wavenumber mapping used in the error bound.
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- this is a binary vector stating whether
%  warning flags arise. These flags tell about which conditions make the
%  final result certainly not guaranteed. One flag is considered arisen
%  when its value is 1. The following list explains the flags in the
%  respective vector order:</li>
%   <ul type="circle">
%                <li>1 : If reaching overbudget. It states whether
%                the max budget is attained without reaching the
%                guaranteed error tolerance.</li>
%                <li>2 : If the function lies outside the cone. In
%                this case, results are not guaranteed. For more
%                information about the cone definition, check the
%                article mentioned below.</li>
%   </ul>
%  </ul>
% </html>
%
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in $[0,1]^d$
% with a prescribed generalized error tolerance. The Walsh-Fourier
% coefficients of the integrand are assumed to be absolutely convergent. If
% the algorithm terminates without warning messages, the output is given
% with guarantees under the assumption that the integrand lies inside a
% cone of functions. The guarantee is based on the decay rate of the
% Walsh-Fourier coefficients. For integration over domains other than
% $[0,1]^d$, this cone condition applies to $f \circ \psi$ (the
% composition of the functions) where $\psi$ is the transformation
% function for $[0,1]^d$ to the desired region. For more details on how the
% cone is defined, please refer to the references below.
%
%% Examples
%
%%
% *Example 1*
%
% Estimate the integral with integrand $f(x) = x_1 x_2$ in the hyperbox $[0,1]^2$:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2); ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 2*
%
% Estimate the integral with integrand $f(x) = x_1^2  x_2^2 x_3^2$
% in the hyperbox $R^3$ where $x_1$, $x_2$ and $x_3$ are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; 
  hyperbox = [-inf(1,3);inf(1,3)]; exactsol = 1;
  q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3, 'fudge',@(m) 10*2.^-m); 
  check = double(abs(exactsol-q) < max(1e-3,1e-3*abs(exactsol)))

%%
% *Example 3*
%
% Estimate the integral with integrand $f(x) = exp(-x_1^2-x_2^2)$ in the
% hyperbox $[-1,2]^2$:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2); 2*ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-2); 
  exactsol = 1/9*(sqrt(pi)/2*(erf(2)+erf(1)))^2;
  check = double(abs(exactsol-q) < max(1e-3,1e-2*abs(exactsol)))

%%
% *Example 4*
%
% Estimate the price of an European call with $S_0=100$, $K=100$, $r=\sigma^2/2$,
% $\sigma=0.05$, and $T=1$.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); 
  hyperbox = [-inf(1,1);inf(1,1)];
  q = cubSobol_g(f,hyperbox,'normal',1e-4,1e-2); 
  price = normcdf(0.05)*100 - 0.5*100*exp(-0.05^2/2);
  check = double(abs(price-q) < max(1e-4,1e-2*abs(price)))

%%
% *Example 5*
%
% Estimate the integral with integrand $f(x) = 8 x_1 x_2 x_3 x_4 x_5$ in the interval
% $[0,1)^5$ with pure absolute error $10^{-5}$.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0); exactsol = 1/4;
  check = double(abs(exactsol-q) < 1e-5)

%%
% *Example 6*
%
% Estimate the integral with integrand $f(x) = x_1^2+x_2^2$ over the disk
% with center $(0,0)$ and radius $1$ with pure absolute error
% $10^{-5}$, where $x = [x_1, x_2]$ is a vector.

  f = @(x) x(:,1).^2+x(:,2).^2; hyperbox = [0,0,1]; exactsol = pi/2;
  q = cubSobol_g(f,hyperbox,'uniform ball','abstol',1e-4,'reltol',0); 
  check = double(abs(exactsol-q) < 1e-4)

%%
% *Example 7*
%
% Estimate the integral with integrand $f(x) = 10 x_1 - 5 x_2^2 + x_3^3$
% in the interval $[0,2)^3$ with pure absolute error $10^{-5}$ using
% two control variates $h_1(x) = x_1$ and $h_2(x) = x_2^2$.

  g.func = @(x) [10*x(:,1)-5*x(:,2).^2+1*x(:,3).^3, x(:,1), x(:,2).^2];
  g.cv = [1,4/3]; hyperbox= [zeros(1,3);2*ones(1,3)];
  q = cubSobol_g(g,hyperbox,'uniform',1e-6,0); exactsol = 16/3;
  check = double(abs(exactsol-q) < 1e-6)


%% See Also
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <p><a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
% <html>
% <p><a href="cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama "Reliable
% adaptive cubature using digital sequences", Monte Carlo and Quasi-Monte
% Carlo Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D.
% Nuyens, eds.), Springer Proceedings in Mathematics and Statistics, vol.
% 163, Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp.
% 367-383.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
% Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
% Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
% from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software," Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable
% Mathematical Software" [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://gailgithub.github.io/GAIL_Dev/
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice and Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% [6] Fang, K.-T., & Wang, Y. (1994). Number-theoretic Methods in
% Statistics. London, UK: CHAPMAN & HALL
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% cubBayesLattice_g
% Bayesian cubature method to estimate the integral of a random variable 
% using rank-1 Lattices over a d-dimensional region within a 
% specified generalized error tolerance with guarantees under Bayesian
% assumptions.
%% Syntax
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%           'order',order,'ptransform',ptransform,'arbMean',arbMean)
%
% [OBJ] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%         'order',order,'ptransform',ptransform,'arbMean',arbMean)
%
% [Q,OutP] = *compInteg*(OBJ)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,absTol,relTol)
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,inParams)
%
%% Description
% 
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,'absTol',absTol,'relTol',relTol,
%   'order',order,'ptransform',ptransform,'arbMean',arbMean) Initializes
%   the object with the given parameters and also returns an
%   estimate of integral Q.
%
% [Q,OutP] = *compInteg*(OBJ) estimates the integral of f over hyperbox
%   $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling to within a specified generalized
%   error tolerance, tolfun = max(abstol, reltol*| I |), i.e., | I - Q | <=
%   tolfun with confidence of at least 99%, where I is the true integral
%   value, Q is the estimated integral value, abstol is the absolute error
%   tolerance, and reltol is the relative error tolerance. Usually the
%   reltol determines the accuracy of the estimation; however, if | I | is
%   rather small, then abstol determines the accuracy of the estimation.
%   Given the construction of our Lattices, d must be a positive integer
%   with 1 <= dim <= 600. For higher dimensions, it is recommended to use 
%   simpler periodization transformation like 'Baker'.
%
% It is recommended to use *compInteg* for estimating the integral repeatedly
% after the object initialization.
% 
% OutP is the structure holding additional output params, more details
% provided below. Input f is a function handle that accepts an n x d
% matrix input, where d is the dimension of the hyperbox, and n is the
% number of points being evaluated simultaneously.
%
% The following additional input parameter passing styles also supported:
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim) estimates the integral of f over
%   hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,absTol,relTol); estimates the integral 
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. All parameters
%   should be input in the order specified above. The answer is given within 
%   the generalized error tolerance tolfun. All other input parameters
%    are initialized with default values as given below.  
%
% [OBJ,Q] = *cubBayesLattice_g*(f,dim,inParms); estimates the integral 
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using rank-1 Lattice sampling. 
%   The structure inParams shall hold the optional input parameters.
%
% *Input Arguments*
%
% * f --- the integrand.
%
% * dim --- number of dimensions of the integrand.
%
% *Optional Input Arguments*
%
% * absTol --- absolute error tolerance | I - Q | <= absTol. Default is 0.01
%
% * relTol --- relative error tolerance | I - Q | <= I*relTol. Default is 0
%
% * order --- order of the Bernoulli polynomial of the kernel r=1,2. 
%             If r==0, algorithm automatically chooses the kernel order
%             which can be a non-integer value.
%             Default is 2
%
% * ptransform --- periodization variable transform to use: 'Baker', 'C0',
%                  'C1', 'C1sin', or 'C2sin'. Default is 'C1sin'
%
% * arbMean --- If false, the algorithm assumes the integrand was sampled
%                 from a Gaussian process of zero mean. Default is 'true'
%
% * alpha --- confidence level for a credible interval of Q. Default is 0.01
%
% * mmin --- min number of samples to start with: 2^mmin. Default is 10
%
% * mmax --- max number of samples allowed: 2^mmax. Default is 22
%
% <html>
% <ul type="square">
% <li>stopCriterion --- stopping criterion to use. Supports three options: </li>
%   <ul type="circle">
%                <li>1) MLE: Empirical Bayes</li>
%                <li>2) GCV: Generalized Cross Validation</li>
%                <li>3) full: Full Bayes</li>
%   </ul>
%    Default is MLE: Empirical Bayes
%  </ul>
% </html>
%
% * useGradient -- If true uses gradient descent in parameter search.
%                   Default is false
%
% * oneTheta -- If true uses common shape parameter for all dimensions,
%                 else allow shape parameter vary across dimensions.
%                 Default is true
%
% *Output Arguments*
%
% * n --- number of samples used to compute the integral of f.
%
% * time --- time to compute the integral in seconds.
%
% <html>
% <ul type="square">
%  <li>exitFlag --- indicates the exit condition of the
%  algorithm:</li>
%   <ul type="circle">
%                <li>1 - integral computed within the error tolerance and
%                      without exceeding max sample limit 2^mmax </li>
%                <li>2 - used max number of samples and yet not met the
%                      error tolerance</li>
%   </ul>
%  </ul>
% </html>
%
% * ErrBd  --- estimated integral error | I - Q |
%
% * optParams --- optional parameters useful to debug and get better
%                  understanding of the algorithm
%
% * optParams.aMLEAll ---- returns the shape parameters computed
%
%
%
%%  Guarantee
%
% This algorithm attempts to calculate the integral of function f over the
% hyperbox $[0,1]^{\mathrm{dim}}$ to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level, e.g., 99% when alpha=0.5%. If the
% algorithm terminates without showing any warning messages and provides
% an answer Q, then the following inequality would be satisfied:
%
% Pr(| Q - I | <= tolfun) = 99%
%
% Please refer to our paper [1] for detailed arguments and proofs.
%
%% Examples
%

%%
% *Example 1: Integrating a simple Quadratic function*
%
% Estimate the integral with integrand $f(x) = x^2$ over the interval $[0,1]$
% with default parameters: order=2, ptransform=C1sin, abstol=0.01, relTol=0

warning('off','GAIL:cubBayesLattice_g:fdnotgiven')
[~,muhat] = cubBayesLattice_g;
exactInteg = 1.0/3;
warning('on','GAIL:cubBayesLattice_g:fdnotgiven')
check = double(abs(exactInteg-muhat) < 0.01)

%%
% *Example 2: ExpCos*
%
% Estimate the integral of Exponential of Cosine function
% $f({x}) = \exp\left(\sum_{i=1}^2cos(2\pi x_i)\right)$ over the
% interval $[0,1]^2$ with parameters: order=2, C1sin variable transform, abstol=0.001,
% relTol=0.01

fun = @(x) exp(sum(cos(2*pi*x), 2));
dim=2; absTol=1e-3; relTol=1e-2;
exactInteg = besseli(0,1)^dim;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'absTol',absTol,'oneTheta',false}];
obj=cubBayesLattice_g(fun,dim,inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 3: Keister function*
%
% Estimate the integral with keister function as integrand over the
% interval $[0,1]^2$ with parameters: order=2, C1 variable transform,
% abstol=0.001, relTol=0.01

dim=3; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',2, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)

%%
% *Example 4:  Multivariate normal probability*
% 
% Estimate the multivariate normal probability for the given hyper interval
% $\left(\begin{array}{c} -6\\ -2\\ -2\end{array}\right) $ and
% $\left(\begin{array}{c} 5\\ 2\\ 1\end{array}\right)$ in $\bf{R}^3$
% having zero mean and covariance
% $\left(\begin{array}{ccc} 4& 1& 1\\ 0& 1& 0.5\\ 0& 0& 0.25 \end{array}\right)$ with
% parameters: order=1, C1sin variable transform,
% abstol=0.001, relTol=0.01

dim=2; absTol=1e-3; relTol=1e-2; fName = 'MVN';
C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
MVNParams.CovProp.C = chol(MVNParams.Cov)';
muBest = 0.676337324357787;
integrand =@(t) GenzFunc(t,MVNParams);
inputArgs={'absTol',absTol,'relTol',relTol};
inputArgs=[inputArgs {'order',1,'ptransform','C1sin','arbMean',true}];
inputArgs=[inputArgs {'useGradient',true}];
[~,muhat]=cubBayesLattice_g(integrand,dim, inputArgs{:});
check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))

%%
% *Example 5: Keister function*
%
% Estimating the Keister integral with Kernel order r chosen automatically
 
dim=2; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'order',0, 'ptransform','C1','arbMean',true}];
obj=cubBayesLattice_g(fKeister,dim,inputArgs{:});
[muhat,outParams] = compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
check = double(outParams.optParams.r > 0)

%%
% *Example 6*
%
% A simple example which uses dimension specific shape parameter

const = [1E-4 1 1E4];
fun = @(x)sum(bsxfun(@times, const, sin(2*pi*x.^2)), 2);
dim=3; absTol=1e-3; relTol=1e-2;
exactInteg = fresnels(2)*sum(const)/2;
inputArgs = {'relTol',relTol, 'order',2, 'ptransform','C1sin'};
inputArgs = [inputArgs {'absTol',absTol,'oneTheta',false,'useGradient',false}];
obj=cubBayesLattice_g(fun, dim, inputArgs{:});
[muhat,outParams]=compInteg(obj);
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))
etaDim = size(outParams.optParams.aMLEAll, 2)


%% See Also
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_cubBayesNet_g.html">cubBayesNet_g</a>
% </html>
%
%% References
%
% [1] Jagadeeswaran Rathinavel, Fred J. Hickernell, Fast automatic Bayesian 
%   cubature using lattice sampling.  Stat Comput 29, 1215-1229 (2019).
%   https://doi.org/10.1007/s11222-019-09895-9
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using
%   matching kernels and designs," PhD thesis, Illinois Institute of Technology, 2019.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% cubBayesNet_g
% Bayesian cubature method to estimate the integral of a random variable 
% using digital nets over a d-dimensional region within a 
% specified generalized error tolerance with guarantees under Bayesian
% assumptions. Currently, only Sobol points are supported.
%% Syntax
% [OBJ,Q] = *cubBayesNet_g*(f,dim,'absTol',absTol,'relTol',relTol,
%     'order',order,'arbMean',arbMean)
%
% [OBJ] = *cubBayesNet_g*(f,dim,'absTol',absTol,'relTol',relTol,...
%     'order',order,'arbMean',arbMean)
%
% [Q,OutP] = *compInteg*(OBJ)
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim)
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim,absTol,relTol)
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim,inParams)
%
%% Description
% 
% [OBJ,Q] = *cubBayesNet_g*(f,dim,'absTol',absTol,'relTol',relTol,...
%     'order',order, 'arbMean',arbMean); initializes the object with the 
%   given parameters and also returns an estimate of integral Q.
%
% [Q,OutP] = *compInteg*(OBJ) estimates the integral of f over hyperbox
%   $[0,1]^{\mathrm{dim}}$ using digital nets (Sobol points) to within a specified
%   generalized error tolerance, tolfun = max(abstol, reltol*| I |), i.e.,
%   | I - Q | <= tolfun with confidence of at least 99%, where I is the
%   true integral value, Q is the estimated integral value, abstol is the
%   absolute error tolerance, and reltol is the relative error tolerance.
%   Usually the reltol determines the accuracy of the estimation; however,
%   if | I | is rather small, then abstol determines the accuracy of the
%   estimation.
% 
% It is recommended to use *compInteg* for estimating the integral
% repeatedly after the object initialization.
%
% OutP is the structure holding additional output params, more details
% provided below. Input f is a function handle that accepts an n x d
% matrix input, where d is the dimension of the hyperbox, and n is the
% number of points being evaluated simultaneously.
%
% The following additional input parameter passing styles also supported:
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim); estimates the integral of f over
%   hyperbox $[0,1]^{\mathrm{dim}}$ using digital nets (Sobol points). All other input parameters
%   are initialized with default values as given below. Returns the initialized
%   object OBJ and the estimate of integral Q.  
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim,absTol,relTol); estimates the integral
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using digital nets (Sobol points). All parameters
%   should be input in the order specified above. The answer is given within
%   the generalized error tolerance tolfun. All other input parameters
%   are initialized with default values as given below.
%
% [OBJ,Q] = *cubBayesNet_g*(f,dim,inParms); estimates the integral 
%   of f over hyperbox $[0,1]^{\mathrm{dim}}$ using digital nets (Sobol points). 
%   The structure inParams shall hold the optional input parameters.
%
%
% *Input Arguments*
%
% * f --- the integrand
%
% * dim --- number of dimensions of the integrand
%
% *Optional Input Arguments*
%
% * absTol --- absolute error tolerance | I - Q | <= absTol. Default is 0.01
%
% * relTol --- relative error tolerance | I - Q | <= I*relTol. Default is 0
%
% * arbMean --- If false, the algorithm assumes the integrand was sampled
%                   from a Gaussian process of zero mean. Default is 'true'
%
% * alpha --- confidence level for a credible interval of Q. Default is 0.01
%
% * mmin --- min number of samples to start with: 2^mmin. Default is 8
%
% * mmax --- max number of samples allowed: 2^mmax. Default is 20
%
% *Output Arguments*
%
% * n --- number of samples used to compute the integral of f.
%
% * time --- time to compute the integral in seconds.
%
% <html>
% <ul type="square">
%  <li>exitFlag --- indicates the exit condition of the algorithm:
%  algorithm:</li>
%   <ul type="circle">
%                <li>1 - integral computed within the error tolerance and
%                      without exceeding max sample limit 2^mmax </li>
%                <li>2 - used max number of samples and yet not met the
%                      error tolerance</li>
%   </ul>
%  </ul>
% </html>
%
% * ErrBd  --- estimated integral error | I - Q |
%
% * optParams --- optional parameters useful to debug and get better
%                  understanding of the algorithm
%
% * optParams.aMLEAll ---- returns the shape parameters computed
%
%
%%  Guarantee
%
% This algorithm attempts to calculate the integral of function f over the
% hyperbox $[0,1]^{\mathrm{dim}}$ to a prescribed error tolerance 
% tolfun:= max(abstol,reltol*| I |) with guaranteed confidence level,
% e.g.,99% when alpha=0.5%. If the algorithm terminates without showing any
% warning messages and provides an answer Q, then the following inequality
% would be satisfied:
%
% Pr(| Q - I | <= tolfun) = 99%
%
% Please refer to our paper [1] for detailed arguments and proofs.
%
%% Examples
%

%%
% *Example 1: Quadratic*
%
% Estimate the integral with integrand $f(x) = x^2$ over the interval
% $[0,1]$ with default parameters: order=1, abstol=0.01, relTol=0

warning('off','GAIL:cubBayesNet_g:fdnotgiven')
[~,muhat] = cubBayesNet_g;
exactInteg = 1.0/3;
warning('on','GAIL:cubBayesNet_g:fdnotgiven')
check = double(abs(exactInteg-muhat) < 0.01)

%%
% *Example 2: ExpCos*
%
% Estimate the integral with integrand
% $f({x}) = \exp\left(\sum_{i=1}^2cos(2\pi x_i)\right)$ over the
% interval $[0,1]^2$ with parameters: order=2, abstol=0.001, relTol=0.01

fun = @(x) exp(sum(cos(2*pi*x), 2));
dim=2; absTol=1e-3; relTol=1e-2;
exactInteg = besseli(0,1)^dim;
inputArgs = {'absTol',absTol,'relTol',relTol};
[~,muhat]=cubBayesNet_g(fun, dim, inputArgs{:});
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))

%%
% *Example 3: Keister function*
%
% Estimate the Keister's integrand, a multidimensional integral  
% inspired by a physics application over the interval $[0,1]^2$ 
% with parameters: order=2, abstol=0.001, relTol=0.01

dim=2; absTol=1e-3; relTol=1e-2;
normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
replaceZeros = @(t) (t+(t==0)*eps); % to avoid getting infinity, NaN
yinv = @(t)(erfcinv( replaceZeros(abs(t)) ));
ft = @(t,dim) cos( sqrt( normsqd(yinv(t)) )) *(sqrt(pi))^dim;
fKeister = @(x) ft(x,dim); exactInteg = Keistertrue(dim);
inputArgs ={'absTol',absTol, 'relTol',relTol};
inputArgs =[inputArgs {'arbMean',true}];
[~,muhat]=cubBayesNet_g(fKeister,dim,inputArgs{:});
check = double(abs(exactInteg-muhat) < max(absTol,relTol*abs(exactInteg)))

%%
% *Example 4: Multivariate normal probability*
%
% For $\bf{X}\sim N(\bf{\mu},\Sigma)$, estimate the following
% probability:
%
% $$ P\left(\bf{a} \leq \bf{X} \leq \bf{b} \right) = \int_{\bf{a}}^{\bf{b}}
% \frac{{\rm e}^{(\bf{x}-\bf{\mu})^T {\Sigma}^{-1}(\bf{x}-\bf{\mu})}}
% {(2\pi)^{d/2}\left|{\Sigma}\right|^{1/2}}\,{\rm d}\bf{x}. $$
%
% Given $\left(\begin{array}{c} -6\\ -2\\ -2\end{array}\right)$ and 
% $\left(\begin{array}{c} 5\\ 2\\ 1\end{array}\right)$  with zero mean $\mu=0$
% and covariance  $\left(\begin{array}{ccc} 4& 1& 1\\ 0& 1& 0.5\\ 0& 0& 0.25 \end{array}\right)$.

C = [4 1 1; 0 1 0.5; 0 0 0.25]; MVNParams.Cov = C'*C; MVNParams.C = C;
MVNParams.a = [-6 -2 -2]; MVNParams.b = [5 2 1]; MVNParams.mu = 0;
MVNParams.CovProp.C = chol(MVNParams.Cov)';
muBest = 0.676337324357787;
integrand =@(t) GenzFunc(t,MVNParams);
inputArgs={'absTol',absTol,'relTol',relTol};
inputArgs=[inputArgs {'arbMean',true}];
obj=cubBayesNet_g(integrand,dim, inputArgs{:});
[muhat,outParams] = compInteg(obj);
check = double(abs(muBest-muhat) < max(absTol,relTol*abs(muBest)))
etaDim = size(outParams.optParams.aMLEAll, 2)


%% See Also
%
% <html>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% </html>
%
% <html>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
% <html>
% <a href="help_cubMC_g.html">cubMC_g</a>
% </html>
%
% <html>
% <a href="help_meanMC_g.html">meanMC_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
% <html>
% <a href="help_cubBayesLattice_g.html">cubBayesLattice_g</a>
% </html>
%
%
%% References
%
% [1] Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using
%   matching kernels and designs," PhD thesis, Illinois Institute of
%   Technology, 2019.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%   Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%   Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%   Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%   from http://gailgithub.github.io/GAIL_Dev/
%
% [3] Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using
%   matching kernels and designs," PhD thesis, Illinois Institute of Technology, 2019.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%% Demos for funappx_g
%
% <html>
% <a href="demo_funappx_g2.html">GUI of funappx_g</a>
% </html>
%
% <html>
% <a href="demo_funappx_g1.html">Approximating a highly fluctuating curve</a>
% </html>
% 

%% Approximate a highly fluctuating curve using *funappx_g*
% Author: Sou-Cheng Choi, July 2017

%% Function definition
%
% Define a highly fluctuating function as follows:
%
% $$ f(x) = x^2 \sin \biggl(\frac{2 \pi}{ x^2} \biggr). $$
%
close all; clear all; format compact; format short;
f = @(x) x.^2 .* sin((2*pi)./x.^2);

%% Function approximation
% We use *funappx_g* to approximate $f$ over the interval $[a,b]$, where
% $a = 0.1$ and $b = 2.5$:
a = 0.1;
b = 2.5;
[q,out] = funappx_g(f, a, b);

%% Plots of the function and approximant
% We plot $f(x)$ and the approximant returned by *funappx_g*, $q(x)$,
% below:
figure;
x = a:1e-6:b;
plot(x,f(x),'r.', x,q(x),'g-');
xlabel('$x$','interpreter','latex')
h_legend=legend('$f(x)$', '$q(x)$');
set(h_legend,'interpreter','latex');
axis tight

%% Plot of the approximation errors
% The following plot shows that all pointwise absolute errors are less than
% the default tolerance of $10^{-6}$.
figure;
semilogy(x,abs(f(x)-q(x)));
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight
max_abs_error = max(abs(f(x)-q(x)))

%% A slightly different example
% If we changes $a$ to a smaller number such as $10^{-2}$, then even if
% we relax the tolerance to $10^{-4}$, *funappx_g* may still return an
% approximant that fails to meet the tolerance. The reason is that $f$ on
% $[a,b]$ is no longer in the cone of functions conducive for successful
% approximation.
a = 1e-2;
abstol = 1e-4;
[q2,out2] = funappx_g(f, a, b, abstol);
figure;
x = a:1e-6:b;
semilogy(x,abs(f(x)-q2(x)));
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight
max_abs_error = max(abs(f(x)-q2(x)))

%% A workaround
% We can widen the cone by increasing the number of initial points given to
% *funappx_g*.
inparam.a = a;
inparam.b = b;
inparam.abstol = abstol;
inparam.ninit = 5e6;
inparam.nmax = inparam.ninit*10;
[q3,out3] = funappx_g(f, inparam);
x = a:1.0/(out3.npoints*2):b;
figure;
semilogy(x,abs(f(x)-q3(x)));
xlabel('$x$','interpreter','latex')
ylabel('absolute error')
axis tight
max_abs_error = max(abs(f(x)-q3(x)))

%% A better way
% Using a large value of |ninit| defeats the purpose of *funappx_g*'s
% locally adaptive design. Notice that the failure region was $[0.01,0.1]$,
% So we can use *funappx_g* with a high value of |ninit| only in this
% region.
inparam.a = a;
inparam.b = 0.1;
inparam.ninit = 2e5;
inparam.nmax =  1e7;
inparam.output_x = 1;
[q4,out4] = funappx_g(f, inparam);

% Use default value of ninit on [0.1,2.5]
inparam.a = inparam.b;
inparam.b = b;
inparam.ninit = 20;
[q5,out5] = funappx_g(f, inparam);

% Define a new approximant on [a,b]
xx = [out4.x, out5.x(2:end)];
yy = [out4.y, out5.y(2:end)];
if gail.matlab_version >= 8.3
    fappx = griddedInterpolant(xx,yy,'linear');
else
    fappx = @(t) ppval(interp1(xx,yy,'linear','pp'), t);
end;

% Evaluate the error again
x = a:1e-7:b;
max_abs_error = max(abs(f(x)-fappx(x)))

%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%% A GUI (graphical user interface) for *funappx_g*
% Author: Yuhan Ding, July 2017
%
% To approximate a peaky function with *funappx_g* and to show how
% *funappx_g* generates grid points for locally adaptive linear spline
% approximation

%% Function definition
%
% Define a peaky function as follows:
%
close all; clear all; format compact; format short;
f = @(x) exp(-1000*(x-0.2).^2);
x = 0:0.0001:1;
figure;
plot(x,f(x))
axis tight

%% Function Approximation
% We use *funappx_g* to approximate $f$ over the interval $[0,1]$ with
% error tolerance $10^{-2}$ and 15 initial subintervals:
[~,out_param] = funappx_g(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15)

% We find that to reach the error tolerance, we need 105 points to
% approximate the function.

%% Process to Generate Grid Points
%
% Step 1: start with $16$ evenly spaced points:
%
% <<localgui1.png>>
%
% Step 2: add points to the peaky part:
%
% <<localgui2.png>>
%
% Step 6: after serveral iterations, the approximation error almost meets the given tolerance:
%
% <<localgui6.png>>
%
% Step 7: the error tolerance is reached:
%
% <<localgui7.png>>
%
%
% This process can also be reproduced by the following command:
% funappx_g_gui(@(x) exp(-1000*(x-0.2).^2),0,1,1e-2,15,15);
%
%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%% Demos for funmin_g
%
% <html>
% <a href="demo_funmin_g2.html">Comparing funmin_g with fminbnd and chebfun</a>
% </html>
%
% <html>
% <a href="demo_funmin_g1.html">Finding global minimum of a highly oscillating function</a>
% </html>
%
%% Find the minimum of a highly oscillating curve using *funmin_g*

%% Function definition
%
% Define a highly oscillating function as follows:
%
% $$ f(x) = \sin (10 \pi x^4 ) - x. $$
%
close all; clearvars; format compact; format short;
f = @(x) sin(10*pi*x.^4)-x;

%% Function minimization
% We use *funmin_g* to approximate $f$ over the interval $[a,b]$, where
% $a = 0$ and $b = 2$ with default parameter values:
a = 0; b = 2;
[fmin,outmin] = funmin_g(f, a, b);

%% Plots of the function and its minimum
% We plot $f(x)$ and the approximate minimum returned by *funmin_g* below.
% It is obvious that the approximation is not satisfactory. We compute the
% error by comparing to the true minimum returned by the Mathematica
% command, |N[Minimize[{Sin[10 Pi x^4] - x, 0 <= x <= 2}, {x}],15]|.  The
% reason is probably that this function is not contained in the cone of
% functions sufficient for successful function minimization.
funmin_g_demo(fmin,outmin)
truefmin=-2.99843616266006;
truexmin=1.99843665971919;
max_abs_error = max(abs(truefmin-fmin))


%% A fix
% We can widen the cone by increasing the number of initial points given to
% *funmin_g*.
inparam.a = a;
inparam.b = b;
inparam.ninit = 1000;
inparam.nmax = inparam.ninit*10;
[fmin2,outmin2] = funmin_g(f, inparam);

funmin_g_demo(fmin2,outmin2)

max_abs_error = max(abs(truefmin-fmin2))

%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%   _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%% Compare *funmin_g* with *fminbnd* and *chebfun*
% Author: Xin Tong, May 2020

%% Function definition and minimization
%
% Define a function with two minima as follows:
%
% $$ f(x) = -5 \exp(-100(x-0.15)^2) - \exp(-80(x-0.65)^2). $$
%
% We use *funmin_g*, *fminbnd*, and *chebfun* to find the minimum of $f$ 
% over the interval $[0,1]$.


%% Set up
close all; clearvars; format compact; format short;
gail.InitializeDisplay
set(0,'defaultLineMarkerSize',15)

%% Plot function
xplot = 0:0.001:1;
fplot = fmin_ex1(xplot);
h(1) = plot(xplot,fplot,'-');
set(h(1),'color',MATLABBlue)
h_legend = legend([h(1)],{'$f(x)$'},'Location','Southeast');
set(h_legend,'interpreter','latex');
hold on

%% Plot minimum values and sample points
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
[ffmg,outfmg] = funmin_g(@fmin_ex1,0,1); 
h(2) = plot(mean(outfmg.intervals),ffmg,'.');
set(h(2),'color',MATLABGreen,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(3) = plot(xAll,fAll,'.');
set(h(3),'color',MATLABGreen)

% fminbnd
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
options = optimset('TolX',outfmg.abstol,'TolFun',outfmg.abstol);
[xfmb,ffmb] = fminbnd(@fmin_ex1,0,1,options);
h(4) = plot(xfmb,ffmb,'.');
set(h(4),'color',MATLABOrange,'MarkerSize',80)
load fmin_ex1X xAll fAll
h(5) = plot(xAll,fAll,'.');
set(h(5),'color',MATLABOrange)

% chebfun
xAll = [];
fAll = [];
save fmin_ex1X xAll fAll
chebf = chebfun(@fmin_ex1,[0,1],'chebfuneps', outfmg.abstol, 'splitting','on');
chebfval = min(chebf); 
chebxvals = roots(diff(chebf));
[v,i] = min(abs(fmin_ex1(chebxvals)-chebfval));
chebxval = chebxvals(i);
chebn = length(chebf);
h(6) = plot(chebxval,chebfval,'o');
set(h(6),'color',MATLABPurple,'MarkerSize',20)
load fmin_ex1X xAll fAll
h(7) = plot(xAll,fAll,'.');
set(h(7),'color',MATLABPurple)
h_legend = legend([h(1) h(2) h(4) h(6) h(3) h(5) h(7)],{'$f(x)$','funmin\_g''s min',...
    'fminbnd''s min','chebfun''s min','funmin\_g''s sample','fminbnd''s sample','chebfun''s sample'},...
    'Location','Southeast');
set(h_legend,'interpreter','latex');

delete('fmin_ex1X.mat');

function y = fmin_ex1(x)
if exist('fmin_ex1X.mat','file')
   load fmin_ex1X xAll fAll
else
   xAll = [];
   fAll = [];
end
xAll = [xAll; x(:)];
y = -5*exp(-100*(x-0.15).^2) - exp(-80*(x-0.65).^2);
fAll = [fAll; y(:)];
save fmin_ex1X xAll fAll
end

%% References
%
% [1] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Xin Tong, "Local
%     Adaption for Approximation and Minimization of Univariate Functions,"
%     _Journal of Complexity_ 40, pp. 17-33, 2017.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/

%% Demos for integral_g
%
% <html>
% <a href="demo_integral_g1.html">Integrate a spiky function</a>
% </html>
%
 %% Integrate a spiky function using *integral_g*
% Authors:  Fred Hickernell and Sou-Cheng Choi, August 2017


%% Function definition
%
% This example is taken from [1], where a function is defined on $[0,1]$
% with twelve spikes.
%
close all; clear all; format compact; format short e;
[~,~,MATLABVERSION] = GAILstart(false);

xquad = 0.13579; %number used by quad to split interval into three parts
xleft = [0 xquad/2 xquad 3*xquad/2 2*xquad];
xctr = [2*xquad 1/4+xquad 1/2 3/4-xquad 1-2*xquad];
xrght = [1-2*xquad 1-3*xquad/2 1-xquad 1-xquad/2 1];
xall = [xleft xctr(2:5) xrght(2:5)]';
nnode = length(xall);

fbump = @(x) 4^3*((x.*(1-x)).^3).*((x>=0)&(x<=1)); %one bump
xplot = (0:0.002:1)'; %points to plot
spikyfun = @(x) foolfunmaker(x, @(x,c) fbump((x-c(1))/c(2)),...
    ones(nnode-1,1), [xall(1:nnode-1) diff(xall)]);

%% Plot of the spiky function
% In the following, we plot $f(x)$ and show the data sampling points
% picked by MATLAB's built-in integration function *quad*, which explains
% why *quad* essentially gives the answer zero for our spiky function:
figure;
h = plot(xplot,spikyfun(xplot), 'k-', xall, zeros(nnode,1), 'k.');
axis([0 1 -0.3 1.1])
set(gca,'Ytick',-0.2:0.2:1)
legend(h,{'$f$','data'},'location','southeast')


%% Integral approximation
% We use MATLAB built-in functions and *integral_g* [2] from GAIL [3] to
% integrate $f$ over the unit interval:
a = 0;
b = 1;
abstol = 1e-11;
if MATLABVERSION >= 8,
    MATintegralspiky = integral(spikyfun,a,b,'AbsTol',abstol)
end
MATquadspiky = quad(spikyfun,a,b,abstol)
MATgailspiky = integral_g(spikyfun,a,b,abstol)


%% Compute approximation errors
% The true integral value of the spiky function is $\frac{16}{35}$. The following
% code computes absolute errors from the above approximation methods. Only
% *integral_g* achieves the required accuracy with respect to the absolute
% tolerance of $10^{-11}$ in this example.
integralspiky = 16/35;
if MATLABVERSION >= 8,
  abs_errors = abs(integralspiky - [MATintegralspiky, MATquadspiky, MATgailspiky])
else
  abs_errors = abs(integralspiky - [MATquadspiky, MATgailspiky])
end
if_meet_abstol = (abs_errors < abstol)

%% References
%
% [1] Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
%     Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic
%     Algorithms: Cones, Not Balls," Journal of Complexity 30, pp. 21-45,
%     2014.
%
% [2] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%% Demos for meanMC_g
%
% <html>
% <a href="count_success.html">Demo for meanMC_g</a>
% </html>%% Counting the success rate of meanMC_g
% Authors: Lan Jiang and Sou-Cheng Choi, July 2017
%
% Define an integration problem as follows:
%
% $$I = \int_0^1 x^2 dx.$$
% 
% The analytical solution is $\frac{1}{3}$. If we use *meanMC_g* to estimate the
% integral with $1000$ replications, we expect the success rate to be bigger
% than or equal to |(1 - alpha)|.

success = 0;
n = 1000;
in_param.reltol = 0; in_param.abstol = 1e-3;
in_param.alpha = 0.05; Yrand = @(n) rand(n,1).^2;
exactsol = 1/3;
for i = 1:n,
    tmu = meanMC_g(Yrand,in_param);
    check = abs(exactsol-tmu) < 1e-3;
    if check == 1,
        success = success + 1;
    end
end
disp(['Over ' num2str(n) ' replications, there are ' num2str(success) ' successes.'])
disp(['The success rate is ' num2str(success/n) ', which is larger than '...
    num2str(1-in_param.alpha) '.'])
%% Demos for cubMC_g
%
% <html>
% <a href="demo_normal_probabilities.html">Computing normal probabilities</a>
% </html>%% Demos for cubSobol_g
%
% <html>
% <a href="demo_normal_probabilities.html">Computing normal probabilities</a>
% </html>%% Demos for cubLattice_g
%
% <html>
% <a href="demo_normal_probabilities_small.html">Computing normal probabilities</a>
% </html>
%
%% Estimation of normal probabilities by *cubSobol_g* and *cubMC_g*
% Authors: Lluis Antoni Jimenez Rugama and Lan Jiang, August 2017

%% Introduction
% For $\bf{X} \sim N(\bf{\mu}, \Sigma)$, we will estimate the following
% probability:
%
% $$ P\left(\bf{a} \leq \bf{X} \leq \bf{b} \right) = \int_{\bf{a}}^{\bf{b}}
% \frac{{\rm e}^{(\bf{x}-\bf{\mu})^T {\Sigma}^{-1}(\bf{x}-\bf{\mu})}}
% {(2\pi)^{d/2}\left|{\Sigma}\right|^{1/2}}\,{\rm d}\bf{x}. $$
%
% We present three tests, each of which approximates the aforementioned
% probability using *cubSobol_g* and *cubMC_g*, which are quasi-Monte Carlo
% and IID Monte Carlo algorithms in GAIL, respectively. In order to
% facilitate the computations when $d$ is high (~30), we are going to apply
% a special transformation of the integrand proposed by Alan Genz.



%% Basic integration parameters set up
% For all the examples, the dimension of the problem is $d=30$.
% The user input tolerances are also set up below: |abstol| is the absolute
% error tolerance, and |reltol| the relative error tolerance. When |reltol|
% is set to 0, the algorithms use pure absolute error bound, and
% vice versa. Finally, for simplicity we define the mean of the distribution
% to be $\bf{\mu}=\bf{0}$:

function demo_normal_probabilities()
d = 30; % Dimension of the problem
abstol = 1e-3; % User input, absolute error bound
reltol = 0;  % User input, relative error bound
mu = zeros(d,1); % Mean of the distribution

%% First test: $\Sigma=I_d$
% For this first example, we consider $\Sigma=I_d$, and
% $\bf{b}=-\bf{a}=(3.5,\dots,3.5)$. In this case, the
% solution of the integral is known so we can verify that the error
% conditions are met:
Sigma = eye(d); % We set the covariance matrix to the identity
factor = 3.5;
hyperbox = [-factor*ones(1,d) ; factor*ones(1,d)]; % We define the integration limits
exactsol = (gail.stdnormcdf(factor)-gail.stdnormcdf(-factor))^d; % Exact integral solution

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 1.1: cubMC_g')
disp(['Estimated probability with cubMC_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.ntot) ' points.'])
disp(['Real error is ' ...
    num2str(abs(exactsol-approx_prob))...
    ' which is less than the user input tolerance '...
    num2str(gail.tolfun(abstol,reltol,1,exactsol,'max')) '.'])

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 1.2: cubSobol_g')
disp(['Estimated probability with cubSobol_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.n) ' points.'])
disp(['Real error is ' ...
    num2str(abs(exactsol-approx_prob))...
    ' which is less than the user input tolerance '...
    num2str(gail.tolfun(abstol,reltol,1,exactsol,'max')) '.'])

%% Second test: $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% For this second example, we consider $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% ($1$ on the diagonal, $0.6$ off the diagonal),
% $\bf{a}=(-\infty,\dots,-\infty)$, and $\bf{b}=\sqrt{d}(U_1,\dots,U_d)$
% ($\bf{b}$ is chosen randomly). The solution for this integral is known
% too so we can verify the real error:
sig = 0.6;
Sigma = sig*ones(d,d); Sigma(1:d+1:d*d) = 1; % set the covariance matrix
hyperbox = [-Inf*ones(1,d) ; sqrt(d)*rand(1,d)]; % define the integration limits
exactsol = integral(@(t)MVNPexact(t,hyperbox(2,:),sig),...
    -inf, inf,'Abstol',1e-8,'RelTol',1e-8)/sqrt(2*pi);

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 2.1: cubMC_g')
disp(['Estimated probability with cubMC_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.ntot) ' points.'])
disp(['Real error is ' ...
    num2str(abs(exactsol-approx_prob))...
    ' which is less than the user input tolerance '...
    num2str(gail.tolfun(abstol,reltol,1,exactsol,'max')) '.'])

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 2.2: cubSobol_g')
disp(['Estimated probability with cubSobol_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.n) ' points.'])
disp(['Real error is ' ...
    num2str(abs(exactsol-approx_prob))...
    ' which is less than the user input tolerance '...
    num2str(gail.tolfun(abstol,reltol,1,exactsol,'max')) '.'])

%% Third test: $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% For this last example, we consider the same covariance matrix in the
% second test but the upper and lower limits are different,
% $\bf{a}=-d/3(U_1,\dots,U_d)$, and $\bf{b}=d/3(U_{d+1},\dots,U_{2d})$
% (both $\bf{a}$ and $\bf{b}$ are chosen randomly):
hyperbox = [-(d/3)*rand(1,d) ; (d/3)*rand(1,d)]; % We define the integration limits

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 3.1: cubMC_g')
disp(['Estimated probability with cubMC_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.ntot) ' points.'])

% Solution approx_prob and integration output parameters in out_param
[approx_prob,out_param] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
disp('Test 3.2: cubSobol_g')
disp(['Estimated probability with cubSobol_g is: ' num2str(approx_prob)])
disp(['The algorithm took ' num2str(out_param.time) ' seconds and '...
    num2str(out_param.n) ' points.'])

%% Appendix: Auxiliary function definitions
% The following functions are defined for the above test examples.
% |multi_normcdf_cubSobol| and |multi_normcdf_cubMC| redefine *cubSobol_g*
% and *cubMC_g* respectively for computing normal probabilites based on
% Alan Genz's transformation. |f| is the function resulting from applying
% Alan Genz's transform that is called in either *cubSobol_g* or *cubMC_g*.

function [p,out, y, kappanumap] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol)
% Using cubSobol_g, multi_normcdf_cubMC computes the cumulative
% distribution function of the multivariate normal distribution with mean
% mu, covariance matrix Sigma and within the region defined by hyperbox.
hyperbox = bsxfun(@minus, hyperbox, mu');
C = chol(Sigma)'; d = size(C,1);
a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);
[p, out, y, kappanumap] = cubSobol_g(...
    @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
    'uniform',abstol,reltol);
end

function [Q,param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol)
% Using cubMC_g, multi_normcdf_cubMC computes the cumulative distribution
% function of the multivariate normal distribution with mean mu, covariance
% matrix Sigma and within the region defined by hyperbox.
hyperbox = bsxfun(@minus, hyperbox, mu');
C = chol(Sigma)'; d = size(C,1);
a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);
[Q,param] = cubMC_g(...
    @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
    'uniform',abstol,reltol);
end

function f_eval = f(s,e,hyperbox,w,C)
% This is the integrand resulting from applying Alan Genz's transformation,
% which is recursively defined.
f_eval = (e-s)*ones(size(w,1),1);
aux = ones(size(w,1),1);
y = [];
for i = 2:size(hyperbox,2);
    y = [y gail.stdnorminv(s+w(:,i-1).*(e-s))];
    aux = sum(bsxfun(@times,C(i,1:i-1),y),2);
    a = (hyperbox(1,i)-aux)/C(i,i);
    b = (hyperbox(2,i)-aux)/C(i,i);
    s = gail.stdnormcdf(a);
    e = gail.stdnormcdf(b);
    f_eval = f_eval .* (e-s);
end
end

function MVNPfunvalfinal = MVNPexact(t,b,sig)
% MVNPexact calculates the true solution of multivariate normal probability
% when the covariance matrix is in a special form: diagonal is 1 and off
% diagonal elements are all the same.
%
% b   - the upper limits of the integral with size 1 x d
% sig - the off diagonal element
% dim - the dimension of the integral
% t   - the variable
MVNPfunval = (gail.stdnormcdf((b(1)+sqrt(sig)*t)/sqrt(1-sig)));
dim =  length(b);
for i =2:dim
    MVNPfunval= MVNPfunval.*(gail.stdnormcdf((b(i)+sqrt(sig)*t)/sqrt(1-sig)));
    %i=i+100;
end
MVNPfunvalfinal = MVNPfunval.*exp(-t.^2/2);
end
end




%% References
%
% [1] Fred J. Hickernell, Lluis Antoni Jimenez Rugama "Reliable adaptive
%     cubature using digital sequences", Monte Carlo and Quasi-Monte Carlo
%     Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D. Nuyens,
%     eds.), Springer Proceedings in Mathematics and Statistics, vol. 163,
%     Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp.
%     367-383.
%
% [2] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
%     "Guaranteed conservative fixed width confidence intervals via Monte
%     Carlo sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012
%     (J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
%     Springer-Verlag, Berlin, pp. 105-128, 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
%     Means of Random Variables, PhD Thesis, Illinois Institute of
%     Technology, 2016.
%% Demos for cubBayesLattice_g
%
% <html>
% <a href="demo_normal_probabilities_small.html">Computing normal probabilities</a>
% </html>
%% Demos for cubBayesNet_g
%
% <html>
% <a href="demo_normal_probabilities_small.html">Computing normal probabilities</a>
% </html>
%
%% Estimation of normal probabilities by multiple integration algorithms in GAIL
% Authors: Lluis Antoni Jimenez Rugama, Lan Jiang, and Jagadeeswaran Rathinavel, April 2019

%% Introduction
% For $\bf{X}\sim N(\bf{\mu},\Sigma)$, we will estimate the following
% probability:
%
% $$ P\left(\bf{a} \leq \bf{X} \leq \bf{b} \right) = \int_{\bf{a}}^{\bf{b}}
% \frac{{\rm e}^{(\bf{x}-\bf{\mu})^T {\Sigma}^{-1}(\bf{x}-\bf{\mu})}}
% {(2\pi)^{d/2}\left|{\Sigma}\right|^{1/2}}\,{\rm d}\bf{x}. $$
%
% We present three tests, each of which approximates the aforementioned
% probability using *cubSobol_g*, *cubMC_g* and *cubBayesLattice_g*,
% *cubLattice_g*, and *cubBayesNet_g*
% which  are quasi-Monte Carlo, IID Monte Carlo and Bayesian cubature
% algorithms respectively in GAIL. In order to facilitate the computations
% when $d$ is high (>4), we are going to apply a special
% transformation of the integrand proposed by Alan Genz.


%% Basic integration parameters set up
% For all the examples, the dimension of the problem is $d=20$.
% The user input tolerances are also set up below: |abstol| is the absolute
% error tolerance, and |reltol| the relative error tolerance. When |reltol|
% is set to 0, the algorithms use pure absolute error bound, and
% vice versa. Finally, for simplicity we define the mean of the distribution
% to be $\bf{\mu}=\bf{0}$:

function demo_normal_probabilities_small(nRep)
d = 20;  % Dimension of the problem
abstol = 1e-3; % User input, absolute error bound
reltol = 0;  % User input, relative error bound
mu = zeros(d,1); % Mean of the distribution
if nargin < 1
   nRep = 10;
end
nTest = 2;
Ivec(nTest) = 0;
approx_prob_MC(nRep,nTest) = 0;
% out_param_MC(nRep,nTest) = 0;
timeMC(nRep,nTest) = 0;
nSampleMC(nRep,nTest) = 0;

approx_prob_sobol(nRep,nTest) = 0;
% out_param_sobol(nRep,nTest) = 0;
timeSob(nRep,nTest) = 0;
nSampleSob(nRep,nTest) = 0;

approx_prob_lat(nRep,nTest) = 0;
% out_param_lat(nRep,nTest) = 0;
timeLat(nRep,nTest) = 0;
nSampleLat(nRep,nTest) = 0;


approx_prob_BayLat(nRep,nTest) = 0;
% out_param_BayLat(nRep,nTest) = 0;
timeBayLat(nRep,nTest) = 0;
nSampleBayLat(nRep,nTest) = 0;

approx_prob_BaySob(nRep,nTest) = 0;
timeBaySob(nRep,nTest) = 0;
nSampleBaySob(nRep,nTest) = 0;

%% First test: $\Sigma=I_d$
% For this first example, we consider $\Sigma=I_d$, and
% $\bf{b}=-\bf{a}=(3.5,\dots,3.5)$. In this case, the
% solution of the integral is known so we can verify that the error
% conditions are met:
Sigma = eye(d); % We set the covariance matrix to the identity
factor = 3.5;
hyperbox = [-factor*ones(1,d) ; factor*ones(1,d)]; % We define the integration limits
exactsol = (gail.stdnormcdf(factor)-gail.stdnormcdf(-factor))^d; % Exact integral solution
Ivec(1) = exactsol;

% Solution approx_prob and integration output parameters in out_param
% Test 1.1: cubMC_g
for k=1:nRep
  [approx_prob_MC(k,1),out_param_MC(k,1)] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
end
timeMC(:,1) = [out_param_MC(:,1).time];
nSampleMC(:,1) = [out_param_MC(:,1).ntot];
report_integration_result('Test 1.1', 'cubMC_g',abstol,reltol,exactsol,... 
  mean(approx_prob_MC(:,1)),mean(timeMC(:,1)),mean(nSampleMC(:,1)))

% Test 1.2: cubLattice_g
for k=1:nRep
  [approx_prob_lat(k,1),out_param_lat(k,1)] = multi_normcdf_cubLat(hyperbox,mu,Sigma,abstol,reltol);
end
timeLat(:,1) = [out_param_lat(:,1).time];
nSampleLat(:,1) = [out_param_lat(:,1).n];
report_integration_result('Test 1.2', 'cubLattice_g',abstol,reltol,exactsol,... 
  mean(approx_prob_lat(:,1)),mean(timeLat(:,1)),mean(nSampleLat(:,1)))

% Test 1.3: cubSobol_g
for k=1:nRep
  [approx_prob_sobol(k,1),out_param_sobol(k,1)] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
end
timeSob(:,1) = [out_param_sobol(:,1).time];
nSampleSob(:,1) = [out_param_sobol(:,1).n];
report_integration_result('Test 1.3', 'cubSobol_g',abstol,reltol,exactsol,... 
  mean(approx_prob_sobol(:,1)),mean(timeSob(:,1)),mean(nSampleSob(:,1)))

% Test 1.4: cubBayesLattice_g
for k=1:nRep
  [approx_prob_BayLat(k,1),out_param_BayLat(k,1)] = multi_normcdf_cubBayesLat(hyperbox,mu,Sigma,abstol,reltol);
end
timeBayLat(:,1) = [out_param_BayLat(:,1).time];
nSampleBayLat(:,1) = [out_param_BayLat(:,1).n];
report_integration_result('Test 1.4', 'cubBayesLattice_g', abstol,reltol,...
   NaN,mean(approx_prob_BayLat(:,1)), (mean(timeBayLat(:,1))), (mean(nSampleBayLat(:,1))))

% Test 1.5: cubBayesNet_g
for k=1:nRep
  [approx_prob_BaySob(k,1),out_param_BaySob(k,1)] = multi_normcdf_cubBayesNet(hyperbox,mu,Sigma,abstol,reltol);
end
timeBaySob(:,1) = [out_param_BaySob(:,1).time];
nSampleBaySob(:,1) = [out_param_BaySob(:,1).n];
report_integration_result('Test 1.5','cubBayesNet_g',abstol,reltol,NaN,...
   mean(approx_prob_BaySob(:,1)),mean(timeBaySob(:,1)),mean(nSampleBaySob(:,1)))


%% Second test: $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% For this second example, we consider $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% ($1$ on the diagonal, $0.6$ off the diagonal),
% $\bf{a}=(-\infty,\dots,-\infty)$, and $\bf{b}=\sqrt{d}(U_1,\dots,U_d)$
% ($\bf{b}$ is chosen randomly). The solution for this integral is known
% too so we can verify the real error:
sig = 0.6;
Sigma = sig*ones(d,d); Sigma(1:d+1:d*d) = 1; % set the covariance matrix
hyperbox = [-Inf*ones(1,d) ; sqrt(d)*rand(1,d)]; % define the integration limits
exactsol = integral(@(t)MVNPexact(t,hyperbox(2,:),sig),...
  -inf, inf,'Abstol',1e-8,'RelTol',1e-8)/sqrt(2*pi);
Ivec(2) = exactsol;

% Solution approx_prob and integration output parameters in out_param
% Test 2.1: cubMC_g
for k=1:nRep
  [approx_prob_MC(k,2),out_param_MC(k,2)] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
end
timeMC(:,2) = [out_param_MC(:,2).time];
nSampleMC(:,2) = [out_param_MC(:,2).ntot];
report_integration_result('Test 2.1','cubMC_g',abstol,reltol,...
  exactsol,mean(approx_prob_MC(:,2)),mean(timeMC(:,2)),mean(nSampleMC(:,2)))

% Test 2.2: cubLattice_g
for k=1:nRep
  [approx_prob_lat(k,2),out_param_lat(k,2)] = multi_normcdf_cubLat(hyperbox,mu,Sigma,abstol,reltol);
end
timeLat(:,2) = [out_param_lat(:,2).time];
nSampleLat(:,2) = [out_param_lat(:,2).n];
report_integration_result('Test 2.2','cubLattice_g',abstol,reltol,...
  exactsol,mean(approx_prob_lat(:,2)),mean(timeLat(:,2)),mean(nSampleLat(:,2)))

% Test 2.3: cubSobol_g
for k=1:nRep
  [approx_prob_sobol(k,2),out_param_sobol(k,2)] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
end
timeSob(:,2) = [out_param_sobol(:,2).time];
nSampleSob(:,2) = [out_param_sobol(:,2).n];
report_integration_result('Test 2.3','cubSobol_g',abstol,reltol,...
  exactsol,mean(approx_prob_sobol(:,2)),mean(timeSob(:,2)),mean(nSampleSob(:,2)))

% Test 2.4: cubBayesLattice_g
for k=1:nRep
  [approx_prob_BayLat(k,2),out_param_BayLat(k,2)] = multi_normcdf_cubBayesLat(...
    hyperbox,mu,Sigma,abstol,reltol);
end
timeBayLat(:,2) = [out_param_BayLat(:,2).time];
nSampleBayLat(:,2) = [out_param_BayLat(:,2).n];
report_integration_result('Test 2.4','cubBayesLattice_g',abstol,reltol,...
  NaN,mean(approx_prob_BayLat(:,2)),mean(timeBayLat(:,2)),mean(nSampleBayLat(:,2)))

% Test 2.5: cubBayesNet_g
for k=1:nRep
  [approx_prob_BaySob(k,2),out_param_BaySob(k,2)] = multi_normcdf_cubBayesNet(...
    hyperbox,mu,Sigma,abstol,reltol);
end
timeBaySob(:,2) = [out_param_BaySob(:,2).time];
nSampleBaySob(:,2) = [out_param_BaySob(:,2).n];
report_integration_result('Test 2.5','cubBayesNet_g',abstol,reltol,...
  NaN,mean(approx_prob_BaySob(:,2)),mean(timeBaySob(:,2)),mean(nSampleBaySob(:,2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ivec = repmat(Ivec,nRep,1);

absErrMC = abs(Ivec-approx_prob_MC);
succMC = mean(absErrMC <= abstol)
avgAbsErrMC = mean(absErrMC)

absErrSob = abs(Ivec-approx_prob_sobol);
succSob = mean(absErrSob <= abstol)
avgAbsErrSob = mean(absErrSob)

absErrLat = abs(Ivec-approx_prob_lat);
succLat = mean(absErrLat <= abstol)
avgAbsErrLat = mean(absErrLat)

absErrBayLat = abs(Ivec-approx_prob_BayLat);
succBayLat = mean(absErrBayLat <= abstol)
avgAbsErrBayLat = mean(absErrBayLat)

absErrBaySob = abs(Ivec-approx_prob_BaySob);
succBaySob = mean(absErrBaySob <= abstol)
avgAbsErrBaySob = mean(absErrBaySob)

timeMC        = mean(timeMC);
timeLat       = mean(timeLat);
timeSob       = mean(timeSob);  
timeBayLat    = mean(timeBayLat);
timeBaySob    = mean(timeBaySob);
nSampleMC     = mean(nSampleMC);
nSampleLat    = mean(nSampleLat);
nSampleSob    = mean(nSampleSob);
nSampleBayLat = mean(nSampleBayLat);
nSampleBaySob = mean(nSampleBaySob);

outFileName = gail.save_mat('Paper_cubBayesLattice_g',['MVNCubExBayesDataNRep' int2str(nRep)],...
    true, abstol, ...
    avgAbsErrMC, avgAbsErrLat, avgAbsErrSob, avgAbsErrBayLat, avgAbsErrBaySob, ...
    succMC, succLat, succSob, succBayLat, succBaySob, ...
    timeMC, timeLat, timeSob, timeBayLat, timeBaySob, ...
    nSampleMC, nSampleLat, nSampleSob, nSampleBayLat, nSampleBaySob);

MVNCubExBayesOut(outFileName)
fprintf('')


%% Third test: $\Sigma=0.4I_d + 0.6\bf{1}\bf{1}^T$
% For this last example, we consider the same covariance matrix in the
% second test but the upper and lower limits are different,
% $\bf{a}=-d/3(U_1,\dots,U_d)$, and $\bf{b}=d/3(U_{d+1},\dots,U_{2d})$
% (both $\bf{a}$ and $\bf{b}$ are chosen randomly):
hyperbox = [-(d/3)*rand(1,d) ; (d/3)*rand(1,d)]; % We define the integration limits

% Solution approx_prob and integration output parameters in out_param
% Test 3.1: cubMC_g
[approx_prob,out_param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol);
report_integration_result('Test 3.1','cubMC_g',abstol,reltol,...
  NaN,approx_prob,out_param.time,out_param.ntot)

% Test 3.2: cubSobol_g
[approx_prob,out_param] = multi_normcdf_cubSobol(hyperbox,mu,Sigma,abstol,reltol);
report_integration_result('Test 3.2','cubSobol_g',abstol,reltol,...
  NaN,approx_prob,out_param.time,out_param.n)

% Test 3.3: cubBayesLattice_g
[approx_prob,out_param] = multi_normcdf_cubBayesLat(hyperbox,mu,...
  Sigma,abstol,reltol);
report_integration_result('Test 3.3','cubBayesLattice_g',abstol,reltol,...
  NaN,approx_prob,out_param.time,out_param.n)

% Test 3.4: cubBayesNet_g
[approx_prob,out_param] = multi_normcdf_cubBayesNet(hyperbox,mu,...
  Sigma,abstol,reltol);
report_integration_result('Test 3.4','cubBayesNet_g',abstol,reltol,...
  NaN,approx_prob,out_param.time,out_param.n)
fprintf('')



%% Appendix: Auxiliary function definitions
% The following functions are defined for the above test examples.
% |multi_normcdf_cubSobol| and |multi_normcdf_cubMC| redefine *cubSobol_g*
% and *cubMC_g* respectively for computing normal probabilites based on
% Alan Genz's transformation. |f| is the function resulting from applying
% Alan Genz's transform that is called in either *cubSobol_g* or *cubMC_g*.

  function [p,out, y, kappanumap] = multi_normcdf_cubSobol(hyperbox,mu,...
      Sigma,abstol,reltol)
    % Using cubSobol_g, multi_normcdf_cubMC computes the cumulative
    % distribution function of the multivariate normal distribution with mean
    % mu, covariance matrix Sigma and within the region defined by hyperbox.
    hyperbox = bsxfun(@minus, hyperbox, mu');
    C = chol(Sigma)'; d = size(C,1);
    a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
    s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);
    [p, out, y, kappanumap] = cubSobol_g(...
      @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
      'uniform',abstol,reltol);
  end

  function [p,out, y, kappanumap] = multi_normcdf_cubLat(hyperbox,mu,...
      Sigma,abstol,reltol)
    % Using cubLattice_g, multi_normcdf_cubLat computes the cumulative
    % distribution function of the multivariate normal distribution with mean
    % mu, covariance matrix Sigma and within the region defined by hyperbox.
    hyperbox = bsxfun(@minus, hyperbox, mu');
    C = chol(Sigma)'; d = size(C,1);
    a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
    s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);
    [p, out, y, kappanumap] = cubLattice_g(...
      @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
      'uniform',abstol,reltol);
  end

  function [p,out] = multi_normcdf_cubBayesLat(hyperbox,mu,Sigma,abstol,reltol)
    % Using cubBayesLattice_g, multi_normcdf_cubBayesLat computes the cumulative
    % distribution function of the multivariate normal distribution with mean
    % mu, covariance matrix Sigma and within the region defined by hyperbox.

    hyperbox = bsxfun(@minus, hyperbox, mu');
    C = chol(Sigma)';
    a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
    s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);

    [~,dim] = size(hyperbox);
    inputArgs = {'absTol',abstol, 'reltol',reltol, ...
      'order',1, 'ptransform','Baker', ....
      'stopAtTol',true, 'stopCriterion','full'...
      'arbMean',true, 'alpha',0.01 ...
      'optTechnique','None'};

    integrand = @(x) f(s,e,hyperbox,x,C);
    inputArgs{end+1} = 'fName'; inputArgs{end+1} = 'MVN';

    objCubBayes=cubBayesLattice_g(integrand,dim, inputArgs{:});
    [p,out]=compInteg(objCubBayes);

  end

  function [p,out] = multi_normcdf_cubBayesNet(hyperbox,mu,Sigma,abstol,reltol)
    % Using cubBayesLattice_g, multi_normcdf_cubBayes computes the cumulative
    % distribution function of the multivariate normal distribution with mean
    % mu, covariance matrix Sigma and within the region defined by hyperbox.

    hyperbox = bsxfun(@minus, hyperbox, mu');
    C = chol(Sigma)';
    a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
    s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);

    [~,dim] = size(hyperbox);
    inputArgs = {'absTol',abstol, 'reltol',reltol, ...
      'order',1, ....
      'stopAtTol',true, 'stopCriterion','full'...
      'arbMean',true, 'alpha',0.01 ...
      'optTechnique','None'};

    integrand = @(x) f(s,e,hyperbox,x,C);
    inputArgs{end+1} = 'fName'; inputArgs{end+1} = 'MVN';

    objCubBayes=cubBayesNet_g(integrand,dim,inputArgs{:});
    [p,out]=compInteg(objCubBayes);

  end

  function [Q,param] = multi_normcdf_cubMC(hyperbox,mu,Sigma,abstol,reltol)
    % Using cubMC_g, multi_normcdf_cubMC computes the cumulative distribution
    % function of the multivariate normal distribution with mean mu, covariance
    % matrix Sigma and within the region defined by hyperbox.
    hyperbox = bsxfun(@minus, hyperbox, mu');
    C = chol(Sigma)'; d = size(C,1);
    a = hyperbox(1,1)/C(1,1); b = hyperbox(2,1)/C(1,1);
    s = gail.stdnormcdf(a); e = gail.stdnormcdf(b);
    [Q,param] = cubMC_g(...
      @(x) f(s,e,hyperbox,x,C), [zeros(1,d-1);ones(1,d-1)],...
      'uniform',abstol,reltol);
  end

  function f_eval = f(s,e,hyperbox,w,C)
    % This is the integrand resulting from applying Alan Genz's transformation,
    % which is recursively defined.
    f_eval = (e-s)*ones(size(w,1),1);
    aux = ones(size(w,1),1);
    y = [];
    for i = 2:size(hyperbox,2);
      y = [y gail.stdnorminv(s+w(:,i-1).*(e-s))];
      aux = sum(bsxfun(@times,C(i,1:i-1),y),2);
      a = (hyperbox(1,i)-aux)/C(i,i);
      b = (hyperbox(2,i)-aux)/C(i,i);
      s = gail.stdnormcdf(a);
      e = gail.stdnormcdf(b);
      f_eval = f_eval .* (e-s);
    end

    f_eval(isnan(f_eval)) = 0; % reset NaN vlaues to zero
  end

  function MVNPfunvalfinal = MVNPexact(t,b,sig)
    % MVNPexact calculates the true solution of multivariate normal probability
    % when the covariance matrix is in a special form: diagonal is 1 and off
    % diagonal elements are all the same.
    %
    % b   - the upper limits of the integral with size 1 x d
    % sig - the off diagonal element
    % dim - the dimension of the integral
    % t   - the variable
    MVNPfunval = (gail.stdnormcdf((b(1)+sqrt(sig)*t)/sqrt(1-sig)));
    dim =  length(b);
    for i =2:dim
      MVNPfunval= MVNPfunval.*(gail.stdnormcdf((b(i)+sqrt(sig)*t)/sqrt(1-sig)));
      %i=i+100;
    end
    MVNPfunvalfinal = MVNPfunval.*exp(-t.^2/2);
  end

  function report_integration_result(testId,algo,abstol,reltol,exactsol,approxsol,timeSec,nSample)
    fprintf('%s: %s\n', testId,algo)
    fprintf('  Estimated probability: %f \n', approxsol)
    if ~isnan(exactsol)
      fprintf('       True probability: %f \n', exactsol)
    end

    fprintf('  The algorithm took %1.3f seconds and %d points \n', timeSec,nSample)

    if ~isnan(exactsol)
      errTol = gail.tolfun(abstol,reltol,1,exactsol,'max');
      errReal = abs(exactsol-approxsol);
      if errReal > errTol
        ME = MException('cubBayesLattice_g_demo:errorExceeded', ...
          'Real error %1.2e exceeds given tolerance %1.2e',errReal,errTol);
        throw(ME)
      else
        fprintf('  Real error is %1.3e, which is less than the tolerance %1.3e\n',...
          errReal, errTol)
      end
    end
  end

end


%% References
%
% [1] Fred J. Hickernell, Lluis Antoni Jimenez Rugama "Reliable adaptive
%     cubature using digital sequences", Monte Carlo and Quasi-Monte Carlo
%     Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D. Nuyens,
%     eds.), Springer Proceedings in Mathematics and Statistics, vol. 163,
%     Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp.
%     367-383.
%
% [2] Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
%     "Guaranteed conservative fixed width confidence intervals via Monte
%     Carlo sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012
%     (J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
%     Springer-Verlag, Berlin, pp. 105-128, 2014.
%
% [3] Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
%     Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
%     Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
%     Integration Library (Version 2.3.1) [MATLAB Software], 2020. Available
%     from http://gailgithub.github.io/GAIL_Dev/
%
% [4] Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
%     Means of Random Variables, PhD Thesis, Illinois Institute of
%     Technology, 2016.
%
% [5] R. Jagadeeswaran and F. J. Hickernell, "Fast Automatic Bayesian
%     cubature using Lattice sampling", In review,
%     Proceedings of Prob Num 2018, Journal of Statistics and Computing,
%     arXiv:1809.09803 [math.NA]
%%% Estimation of the Expected Value Using |meanMC_CLT|
% For a random variable $Y = f(\mathbf{X})$, we will estimate its expected value:
% 
% $$ \mu=\mathbf{E}(Y) = \mathbf{E}[f(\mathbf{X})]=\int_{\mathbf{R}^d}
%  f(\mathbf{X}) \rho(\mathbf{X})\,\mathrm{d}\mathbf{X}
%  \approx \frac{1}{n}\sum_{i=1}^{n}{f(x_i)}, \qquad {x_i} \, \mathrm{IID} \sim
%  \, \rho $$
% 
%
% We will approximate $\mu$ using the meanMC_CLT GAIL algorithm. It is an
% IID Monte Carlo algorithm using the Central Limit Theorem to determine
% the sample size required to meet the desired error tolerance. In order to
% improve computation efficiency, we will also use control variates.
%%

function demo_meanMC_CLT
%% Initialize the workspace and set the display parameters
% This script cleans up the workspace and makes the display beautiful.

gail.InitializeWorkspaceDisplay
%initialize the workspace and the display parameters

%% Example 1: Estimate $\mathbf{E}[f(X)]$ where $f(X)=\exp(-X^2)$ and $X \sim \mathcal{U} (0,1)$ using $x \mapsto x$ as a control variate
% In this example a function that cannot be integrated analytically is
% integrated using our adaptive IID Monte Carlo method.  We also use the
% function  $x \mapsto x$ as a control variate.
%
% First we plot our function

f = @(x)[exp(-x.^2), x]; YXn = @(n)f(rand(n,1)); 
%set up the random variable 
figure %plot f(x)
x = (0:0.001:1)'; %x values for plot
fx = f(x);
plot(x,fx(:,1),'-'); 
ylabel('$\exp(-x^2)$')
xlabel('$x$')

%%
% Next we set up the parameters for Monte Carlo integration

absTol = 1e-3; %absolute tolerance 
relTol = 0; %relative tolerance
s = struct('Y',YXn,'nY',1,'trueMuCV',1/2); 
% create a structure containing random variables, 
% number of random variables and mean of the control variates 
[hmu,out] = meanMC_CLT(s,absTol,relTol);  % calculate the mean
exactsol = erf(1)*sqrt(pi)/2; %true mean
disp('Example 1')
disp(['Estimated mean is: ' num2str(hmu)])
disp(['     True mean is: ' num2str(exactsol)])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nSample) ' points.'])
disp(['Real error was ' ...
    num2str(abs(exactsol-hmu))...
    ' which is less than the user input tolerance '...
    num2str(absTol) '.'])


%% Example 2: Price European Call Option
% This time we price a European call option using the stock price as a
% control variate.
%
% Initialize option parameters for a European call option
inp.timeDim.timeVector = 0:0.25:1; %time increments
inp.assetParam.initPrice = 10; %initial stock price
inp.assetParam.interest = 0.01; %risk-free interest rate
inp.assetParam.volatility = 0.5; %volatility
inp.payoffParam.strike = 10; %strike price
inp.priceParam.absTol = 0.01; %absolute tolerance of a penny
inp.priceParam.relTol = 0; %relative tolerance
EuroCall = optPayoff(inp); %create a European call option payoff object
EuroCallPayoff = @(n) genOptPayoffs(EuroCall,n); %identify the payoff function

%%
% Plot an empirical distribution of the European call option
n = 1e4; %number of payoffs to plot
payoffs = EuroCallPayoff(n); %generate n payoffs
sortedpay = sort(payoffs); %sort them
figure
plot(sortedpay,((1:n)-1/2)/n,'-'); %plot the empirical distribution function scenarios
xlabel('Payoff in dollars')
ylabel('CDF')
axis([0 50 0 1])

%%
% Note that the option has a positive payoff only about 60% of the time.

%%
% Next we price the option using simple IID Monte Carlo.  We happen to have
% a formula for the exact price that can be used to check our error.
[hmu,out] = meanMC_CLT(EuroCallPayoff,inp.priceParam.absTol,inp.priceParam.relTol);
disp(['Estimated price is: ' num2str(hmu)])
disp(['    Exact price is: ' num2str(EuroCall.exactPrice)])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nSample) ' points.'])
disp(['Real error was ' ...
    num2str(abs(EuroCall.exactPrice-hmu))...
    ' which is less than the user input tolerance '...
    num2str(inp.priceParam.absTol) '.'])

%%
% Now we use the stock price as a control variate.  To do that, we need to
% modify the option payoff object by adding the stock price.

EuroCallCV = optPayoff(EuroCall); %make a copy of the European call option parameters
EuroCallCV.payoffParam = struct('optType', {{'euro','stockprice'}}, ...
   'putCallType', {{'call',''}}); %identify the option type
EuroCallCVPayoff = @(n) genOptPayoffs(EuroCallCV,n); %identify the payoff function
s=struct('Y',EuroCallCVPayoff,'nY',1,'trueMuCV',EuroCallCV.assetParam.initPrice);
[hmu,out] = meanMC_CLT(s,inp.priceParam.absTol,inp.priceParam.relTol);
disp(['Estimated price is: ' num2str(hmu)])
disp(['    Exact price is: ' num2str(EuroCallCV.exactPrice(1))])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nSample) ' points.'])
disp(['Real error was ' ...
    num2str(abs(EuroCallCV.exactPrice(1)-hmu))...
    ' which is less than the user input tolerance '...
    num2str(inp.priceParam.absTol) '.'])

%%
% Let's compare how meanMC_CLT works on pricing a European call option with
% and without a control variate. While the two methods both give answers
% within the error tolerance, adding a control variate uses only about 1/9
% of sample points and takes 1/5 of the time.
%
% Acutally, in GAIL we have a function genOptPrice that could direcelty compute a call
% option price.

%%
% First we compute the answer without control variates
inp.priceParam.cubMethod = 'IID_MC_CLT'; %set method
EuroCall = optPrice(inp); %create an object for computation of the price
[EuroCallPrice, out] = genOptPrice(EuroCall); %compute the option price
disp(['Estimated price is: ' num2str(EuroCallPrice)])
disp(['    Exact price is: ' num2str(EuroCall.exactPrice)])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nPaths) ' points.'])
disp(['Real error was ' ...
    num2str(abs(EuroCall.exactPrice-EuroCallPrice))...
    ' which is less than the user input tolerance '...
    num2str(inp.priceParam.absTol) '.'])
%%
% And then we compute the answer with control variates
EuroCallCV = optPrice(EuroCall); %make a copy of the European call option parameters
EuroCallCV.payoffParam = struct('optType', {{'euro','stockprice'}}, ...
   'putCallType', {{'call',''}}); %identify the option type
[EuroCallCVPrice, out] = genOptPrice(EuroCallCV); %compute the option price
disp(['Estimated price is: ' num2str(EuroCallCVPrice(1))])
disp(['    Exact price is: ' num2str(EuroCallCV.exactPrice(1))])
disp(['The algorithm took ' num2str(out.time) ' seconds and '...
    num2str(out.nPaths) ' points.'])
disp(['Real error was ' ...
    num2str(abs(EuroCallCV.exactPrice(1)-EuroCallCVPrice(1)))...
    ' which is less than the user input tolerance '...
    num2str(inp.priceParam.absTol) '.'])


%% Example 3: Keister's multidimensional integration
% We will evaluate the Keister's integral $I$ using meanMC_CLT.  Note
% that the we do a change of variable $\mathbf{t} = \mathbf{X}/a$
% and transform the integral:
%
% $$
% \begin{array}{cl}
% I &= \int_{\mathbf{R}^d} \cos(\| \mathbf{X} \|) 
% \exp(-\| \mathbf{X} \|^2) \, \mathrm{d} \mathbf{X} \\
% &= \int_{\mathbf{R}^d} a^d \cos(a \|  \mathbf{t}
% \| ) \exp(-a^2 \|  \mathbf{t} \| ^2) \, \mathrm{d}
% \mathbf{t}, \qquad a > 0, \\
% & = \int_{\mathbf{R}^d} \underbrace{(2\pi a^2)^{d/2} 
% \cos(a \| \mathbf{t} \|) \exp((1/2-a^2) \| \mathbf{t} \|^2)}_{f(\mathbf{t})}
% \times \underbrace{\frac{\exp(-\| \mathbf{t} \|^2/2)}
% {(2\pi)^{d/2}}}_{\varrho(\mathbf{t})} \, \mathrm{d} \mathbf{t} \\
% & = \mathbf{E}[f(\mathbf{T})], \qquad \mathrm{where} \quad \mathbf{T} 
% \sim \mathcal{N}(\mathbf{0},\mathbf{I}). 
% \end{array}
% $$
% 
% 
%%
% define an anonymous function $f$ as follows:

normsqd = @(t) sum(t.*t,2); %squared l_2 norm of t
f1 = @(normt,a,d) ((2*pi*a^2).^(d/2)) * cos(a*sqrt(normt)) ...
   .* exp((1/2-a^2)*normt);
f = @(t,a,d) f1(normsqd(t),a,d);
abstol = 0; %absolute tolerance
reltol = 0.01; %relative tolerance
dvec = 1:5; %vector of dimensions
avec = [1 1/sqrt(2) 1/sqrt(1.5) 1/sqrt(3)]; %trial values of a 
IMCvec(size(dvec)) = 0; %vector of answers
f2 = @(t,d) cell2mat(arrayfun(@(a) f(t,a,d),avec,'UniformOutput',false)); %a vector of functions for each value of a
outT = zeros(size(dvec)); %vector of time
outN = zeros(size(dvec)); %vector of points
for d = dvec
   f3 = @(t)f2(t,d);%integration in dimension d
   YXn = @(n)f3(randn(n,d));%random generator
   s = struct('Y',YXn,'nY',size(avec,2)); 
   [IMCvec(d),out] = meanMC_CLT(s,abstol,reltol);
   outT(d) = out.time;
   outN(d) = out.nSample;
end
Ivec = zeros(size(dvec)); %vector of true integration
for d = dvec
   Ivec(d) = Keistertrue(d); %true integration
end
error=abs(Ivec-IMCvec)./abs(Ivec); %calculate the relative error
%display a result table
title ='Dimension  Est. Mean  True Mean  Rel. Error     N      Time';
datasave=[dvec; IMCvec; Ivec; error; outN; outT];
disp(title)
fprintf('    %1.0f        %2.4f     %2.4f     %2.4f  %8.0f   %2.4f \n', datasave)


%%
% Author _Yueyi Li_
