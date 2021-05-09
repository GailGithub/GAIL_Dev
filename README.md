[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4018190.svg)](https://doi.org/10.5281/zenodo.4018190)


Guaranteed Automatic Integration Library (GAIL)
-----------------------------------------------
GAIL Version 2.3.2, 2021.

See LICENSE.m for copyright and disclaimer. Refer to ReleaseNotes.m for
what is new in this version.

GAIL is a suite of algorithms for integration problems in one and many
dimensions, and whose answers are guaranteed to be correct.


Developed by
-------------

Fred Hickernell, Sou-Cheng Choi, and their collaborators including Yuhan
Ding, Lan Jiang, Lluis Antoni Jimenez Rugama, Da Li, Jagadeeswaran
Rathinavel, Kan Zhang, Yizhi Zhang, and Xuan Zhou, Department of Applied
Mathematics, Illinois Institute of Technology (IIT) and Xin Tong,
Department of Mathematics, Statistics, and Computer Science, University
of Illinois at Chicago.

We thank the contributions of Aleksei Sorokin, Noah Grudowski, Francisco
Hernandez, Cu Hauw Hung, Yueyi Li, Xincheng Sheng, Xiaoyang Zhao, Tianci
Zhu, and the IIT classes of SCI 498 Adaptive Monte Carlo Algorithms with
Applications to Financial Risk Management, Summer 2016; MATH 491 Reading
& Research, Summer 2015; SCI 498/MATH 491 Computational Social Sciences,
Summer 2016; MATH 491-195 Solving Problems in the Social Sciences Using
Tools from Computational Mathematics and Statistics, Summer 2015; Math
573 Reliable Mathematical Software, Fall 2013.



Please cite the following software, papers, and materials:

Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
Antoni Jimenez Rugama, Da Li, Jagadeeswaran Rathinavel, Xin Tong, Kan
Zhang, Yizhi Zhang, and Xuan Zhou, GAIL: Guaranteed Automatic
Integration Library (Version 2.3.2) [MATLAB Software], 2021. Available from
http://gailgithub.github.io/GAIL_Dev/
(this software)

Sou-Cheng T. Choi, Fred J. Hickernell, Jagadeeswaran Rathinavel, Michael McCourt,
and Aleksei Sorokin, QMCPy: A quasi-Monte Carlo Python library, 2020. Available
from https://qmcsoftware.github.io/QMCSoftware/. Working.
(open-source Python package for Quasi-Monte Carlo methods with some GAIL
algorithms)

Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible Research via
Supportable Scientific Software," Journal of Open Research Software, Volume 2,
Number 1, e22, pp. 1-7, 2014.
(describes principles of Reliable Reproducible Research and Supportable
Scientific Software)

Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Xin Tong, "Local
Adaption for Approximation and Minimization of Univariate Functions,"
Journal of Complexity 40, pp. 17-33, 2017.
(describes funappx_g.m and funmin_g.m)

Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable Mathematical
Software" [Course Slides], Illinois Institute of Technology, Chicago, IL, 2013.
Available from http://gailgithub.github.io/GAIL_Dev/
(develops practices of Reliable Reproducible Research and Supportable
Scientific Software)

Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
(describes deprecated integralTrap_g.m and deprecated funappxtau_g.m)

Yuhan Ding, "Guaranteed Adaptive Univariate Function Approximation," PhD
thesis, Illinois Institute of Technology, 2015.
(describes deprecated funappxPenalty_g.m)

Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
"Guaranteed conservative fixed width confidence intervals via Monte
Carlo sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012
(J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
Springer-Verlag, Berlin, pp. 105-128, 2014.
(describes meanMC_g.m and cubMC_g.m)

Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable adaptive
cubature using digital sequences", Monte Carlo and Quasi-Monte Carlo
Methods: MCQMC, Leuven, Belgium, April 2014 (R. Cools and D. Nuyens, eds.),
Springer Proceedings in Mathematics and Statistics, vol. 163,
Springer-Verlag, Berlin, 2016, arXiv:1410.8615 [math.NA], pp. 367-383.
(describes cubSobol_g.m)

Fred J. Hickernell, Lluis Antoni Jimenez Rugama, and Da Li, "Adaptive
quasi-Monte Carlo methods, 2017+, submitted for publication,
arXiv:1702.01491 [math.NA].

Yizhi Zhang, Guaranteed Adaptive Automatic Algorithms for
Univariate Integration: Methods, Costs and Implementations, PhD
Thesis, Illinois Institute of Technology, 2018.
(describes integral_g.m)

Lan Jiang, Guaranteed Adaptive Monte Carlo Methods for Estimating
Means of Random Variables, PhD Thesis, Illinois Institute of
Technology, 2016.
(describes meanMC_g.m and cubMC_g.m)

Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive
multidimensional integration based on rank-1 lattices," Monte Carlo
and Quasi-Monte Carlo  Methods: MCQMC, Leuven, Belgium, April 2014
(R. Cools and D. Nuyens, eds.), Springer Proceedings in Mathematics
and Statistics, vol. 163, Springer-Verlag, Berlin, 2016, arXiv:1411.1966,
pp. 407-422.
(describes cubLattice_g.m)

Fred J. Hickernell, Sou-Cheng T. Choi, Lan Jiang, and Lluis Antoni Jimenez
Rugama, "Monte Carlo simulation, automatic stopping criteria for", Wiley
StatsRef: Statistics Reference Online, pp. 1-7, 2018.
(review of cubMC_g, cubLattice_g, and cubSobol_g)

Jagadeeswaran Rathinavel and Fred J. Hickernell, "Fast automatic Bayesian
cubature using lattice sampling," Statistics and Computingt 29,
pp. 1215-1229, 2019. Available from https://doi.org/10.1007/s11222-019-09895-9
(describes cubBayesLattice_g)

Jagadeeswaran Rathinavel, "Fast automatic Bayesian cubature using matching
kernels and designs," PhD thesis, Illinois Institute of Technology, 2019.
(describes cubBayesLattice_g and cubBayesNet_g)

Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
Workshop On Sustainable Software for Science: Practice and Experiences
(WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
pp. 1-21, 2014.
(discusses practice and challenges in Sustainable Software for Science)

Daniel S. Katz, Sou-Cheng T. Choi, Nancy Wilkins-Diehr, Neil Chue Hong,
Colin C. Venters, James Howison, Frank J. Seinstra, Matthew Jones, Karen
Cranston, Thomas L. Clune,  Miguel de Val-Borro, and Richard Littauer,
"Report on the Second Workshop on Sustainable Software for Science:
Practice and Experiences (WSSSPE2)," Journal of Open Research Software,
Volume 4, Number 1, e7, 2016.
(discusses practice and challenges in Sustainable Software for Science)

Daniel S. Katz, Sou-Cheng T. Choi, Kyle E. Niemeyer, James
Hetherington, Frank Loffler, Dan Gunter, Ray Idaszak, Steven R. Brandt,
Mark A. Miller, Sandra Gesing, Nick D. Jones, Nic Weber, Suresh Marru,
Gabrielle Allen, Birgit Penzenstadler, Colin C. Venters, Ethan Davis,
Lorraine Hwang, Ilian Todorov, Abani Patra, and Miguel de Val-Borro,
"Report on the Third Workshop on Sustainable Software for Science:
Practice and Experiences (WSSSPE3)," Journal of Open Research Software,
Volume 4, Number 1, e37, 2016.
(discusses practice and challenges in Sustainable Software for Science)

Da Li, "Reliable Quasi-Monte Carlo with Control Variates," Master's thesis,
Illinois Institute of Technology, 2016.
(describes cubSobol_g.m for control variates)

Arfon M. Smith, Daniel S. Katz, Kyle E. Niemeyer, and FORCE11 Software
Citation Working Group, "Software citation principles," PeerJ Computer
Science, Volume 2, e86, 2016.
(motivates research software citation)

Xin Tong, "A Guaranteed, Adaptive, Automatic Algorithm for Univariate
Function Minimization," MS thesis, Illinois Institute of Technology, 2014.
(describes deprecated funmin01_g.m)


Downloads
------------
GAIL can be downloaded from http://gailgithub.github.io/GAIL_Dev/.

Alternatively, you can get a local copy of the GAIL repository with
this command:

  git clone https://github.com/GailGithub/GAIL_Dev.git


Requirements
------------

You will need to install MATLAB 7 or a later version.

GAIL is developed in MATLAB versions R2016a to R2020a. In particular, three
of our core algorithms, cubSobol_g, cubBayesNet_g, and cubBayesLattice_g
require the following MATLAB add-on toolboxes: Signal Processing Toolbox,
Optimization Toolbox, Statistics and Machine Learning Toolbox. In MATLAB,
we could use the following command to find out toolbox dependencies of an
algorithm:

  names = dependencies.toolboxDependencyAnalysis({'cubBayesNet_g'})

For development and testing purposes, we use the third-party toolboxes,
Chebfun version 5.7.0 and Doctest for MATLAB, version 2010.

Documentation
-------------

Detailed documentation is available in the folder,
GAIL_Matlab/Documentation/html/GAIL.html.

You can also go to MATLAB's Help. Under the section of Supplemental Software,
you will find GAIL Toolbox's searchable HTML documentation.

A PDF version of GAIL's documentation with selected examples is available at
https://github.com/GailGithub/GAIL_Dev/blob/master/Documentation/gail_ug_2_3_1.pdf.

General Usage Notes
-------------------

GAIL Version 2.3.2 includes the following ten algorithms:

1.  funappx_g: One-dimensional function approximation on bounded interval

2.  funmin_g: global minimum value of univariate function on a closed interval

3.  integral_g: One-dimensional integration on bounded interval

4.  meanMC_g: Monte Carlo method for estimating mean of a random variable

5.  cubMC_g: Monte Carlo method for numerical multiple integration

6.  cubSobol_g: Quasi-Monte Carlo method using Sobol' cubature for
d-dimensional integration

7.  cubLattice_g: Quasi-Monte Carlo method using rank-1 Lattices cubature
for d-dimensional integration

8. cubBayesLattice_g: Bayesian cubature method for d-dimensional integration
   using lattice points

9. cubBayesNet_g: Bayesian cubature method for d-dimensional integration
   using Sobol points

10.  meanMC_CLT: Monte Carlo method with Central Limit Theorem (CLT)
confidence intervals for estimating mean of a random variable

Each one of our key GAIL algorithms, with the exception of
cubBayesLattice_g and cubBayesNet_g, can parse inputs with the following
three patterns of APIs, where f is a real-valued MATLAB function or
function handle; in_param and out_param are MATLAB structure arrays; and
x is an estimated output:

1. Ordered input values: [x, out_param] = algo(f, inputVal1, inputVal2, inputVal3,...)

2. Input structure array: [x, out_param] = algo(f, in_param)

3. Ordered input values, followed by optional name-value pairs:
	[x, out_param] = algo(f, 'input2', inputVal2, 'input3', inputVal3,...)

For object classes cubBayesLattice_g and cubBayesNet_g, the output pattern
is [out, x], where out is an instance of the corresponding object class.

Installation Instruction
------------------------

1.  Unzip the contents of the zip file to a directory and maintain the
    existing directory and subdirectory structure. (Please note: If you
    install into the "toolbox" subdirectory of the MATLAB program
    hierarchy, you will need to click the button "Update toolbox path
    cache" from the File/Preferences... dialog in MATLAB.)

2.  In MATLAB, add the GAIL directory to your path. This can be done
    by running "GAIL_Install.m".  Alternatively, this can be done by
    selecting "File/Set Path..." from the main or Command window
    menus, or with the command "pathtool". We recommend that you
    select the "Save" button on this dialog so that GAIL is on the
    path automatically in future MATLAB sessions.

3.  To check if you have installed GAIL successfully, type "help
    funappx_g" to see if its documentation shows up.

Alternatively, you could do this:

1.  Download DownloadInstallGail_2_3_2.m and put it where you want
    GAIL to be installed.

2.  Execute it in MATLAB.

To uninstall GAIL, execute "GAIL_Uninstall".

To reinstall GAIL, execute "GAIL_Install".

Tests
-----

We provide quick doctests for each of the functions above. To run
doctests in funappx_g, for example, issue the command doctest
funappx_g.

We also provide unit tests for MATLAB version 8 or later. To run unit
tests for funmin_g, for instance, execute run(ut_funmin_g).

We execute automated nightly fast tests and weekly long tests
on our server. Moreover, these tests are now conducted for all MATLAB
versions from R2016a to R2020a. The test reports are available on Mega
cloud storage at https://mega.nz/. More specifically, fast and long test
reports are archived in text files, gail_daily_tests*
and gail_weekly_tests* at
https://mega.nz/folder/FlMEjI5a#jVixXyAoI05ppbCstz8yEg
respectively. Output files such as images of test scripts are archived
at https://mega.nz/folder/I0cAEKJD#AyQ_8tmxkknfIsuEW0_jnA
respectively.

Known Issues
-------------

During our documentation development with MATLAB releases 2019a and 2020a, the
software's internal HTML viewer is found to display LaTeX expression in
larger font size than it is intended to be. This is an aesthetic issue with
no impact on the content accuracy.  Users may use a web browser to view our
HTML documentation instead. The main page to GAIL's HTML documentation is
GAIL.html, located in the subfolder, Documentation/html/.


Contact Information
--------------------

Please send any queries, questions, or comments to
gail-users@googlegroups.com or visit our project website:
http://gailgithub.github.io/GAIL_Dev/


Acknowledgements
----------------

Our work was supported in part by the National Science Foundation under
grants NSF-DMS-1522687 and NSF-DMS-1115392; and by the Office of Advanced
Scientific Computing Research, Office of Science, U.S. Department of
Energy, under contract DE-AC02-06CH11357.


Organizations
--------------

<img src="./Documentation/logo/illinois-institute-of-technology-vector-logo.jpg" alt="IIT logo"/>

<img src="./Documentation/logo/kamakura-corporation-vector-logo.png" alt="Kamakura logo"/>

