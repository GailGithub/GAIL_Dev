Guaranteed Automatic Integration Library (GAIL)
GAIL Version 2.1, Mar 14, 2015.
See LICENSE.m for copyright and disclaimer.

GAIL is a suite of algorithms for integration problems in one and many
dimensions, and whose answers are guaranteed to be correct.


Developed by
-------------

Fred Hickernell, Sou-Cheng Choi, and their collaborators including
Yuhan Ding, Lan Jiang, Lluis Antoni Jimenez Rugama, Yizhi 
Zhang and Xuan Zhou, Department of Applied Mathematics, Illinois 
Institute of Technology (IIT) and  Xin Tong, Department of Mathematics, 
Statistics, and Computer Science, University of Illinois at Chicago. 
We thank the contributions of Xincheng Sheng, and the IIT class of Math 
573 Reliable Mathematical Software, Fall 2013.


Please cite the following software, papers, and materials:


Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang, and Xuan Zhou, "GAIL: 
Guaranteed Automatic Integration Library (Version 2.1)" [MATLAB Software],
2015. Available from http://code.google.com/p/gail/

Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible Research via 
Supportable Scientific Software", Journal of Open Research Software, Volume 2, 
Number 1 (2014), e22, pp. 1-7, DOI: http://dx.doi.org/10.5334/jors.bb
(describes principles of Reliable Reproducible Research and 
Supportable Scientific Software)

Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable Mathematical
Software" [Course Slides], Illinois Institute of Technology, Chicago, IL, 2013.
(develops practices of Reliable Reproducible Research and 
Supportable Scientific Software)

Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
Cones, Not Balls", Journal of Complexity 30 (2014), pp. 21-45.
(describes integral_g.m and funappx_g.m)

Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
"Guaranteed conservative fixed width confidence intervals via Monte
Carlo sampling", Monte Carlo and Quasi-Monte Carlo Methods 2012
(J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
Springer-Verlag, Berlin, 2014, pp. 105-128. 
(describes meanMC_g.m and cubMC_g.m)

Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable Adaptive 
Cubature Using Digital Sequences", submitted for publication, 2014.
(describes cubSobol_g.m)

Lan Jiang and Fred J. Hickernell, "Guaranteed Conservative Confidence Intervals 
for Means of Bernoulli Random Variables", submitted for publication, 2014.
(describes meanMCBer_g)

Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive Multidimensional 
Integration Based on Rank-1 Lattices", submitted for publication, 2014.
(describes cubLattice_g.m)

Xin Tong, "A Guaranteed, Adaptive, Automatic Algorithm for Univariate 
Function Minimization". MS thesis, Illinois Institute of Technology, 2014. 
(describes funmin_g.m)

Downloads
------------
GAIL can be downloaded from http://code.google.com/p/gail/.

Alternatively, you can get a local copy of the GAIL repository with
this command:

  git clone https://github.com/GailGithub/GAIL_Dev.git


Requirements
------------

You will need to install MATLAB 7 or a later version.


Documentation
-------------

Detailed documentation is available at GAIL_Matlab/Documentation.


General Usage Notes
-------------------

GAIL Version 2.1 includes the following eight algorithms:

1.  funappx_g: One-dimensional function approximation on bounded interval
2.  integral_g:  One-dimensional integration on bounded interval
3.  meanMC_g:  Monte Carlo method for estimating mean of a random variable
4.  cubMC_g: Monte Carlo method for numerical multiple integration
5.  meanMCBer_g:  Monte Carlo method to estimate the mean of a Bernoulli random variable
6.  funmin_g: global minimum value of univariate function on a closed interval
7.  cubSobol_g: Quasi-Monte Carlo method using Sobol' cubature for a d-dimensional integration
8.  cubLattice_g: Quasi-Monte Carlo method using rank-1 Lattices cubature for a d-dimensional integration

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

1.  Download DownloadInstallGail_2_1.m and put it where you want
    GAIL to be installed.

2.  Execute it in MATLAB.

To uninstall GAIL, execute "GAIL_Uninstall".

To reinstall GAIL, execute "GAIL_Install".

Contact Information
--------------------

Please send any queries, questions, or comments to
gail-users@googlegroups.com or visit our project website:
http://code.google.com/p/gail/


Acknowledgements
----------------

Our work was supported in part by grants from the National Science
Foundation under grant NSF-DMS-1115392, and the Office of Advanced
Scientific Computing Research, Office of Science, U.S. Department of
Energy, under contract DE-AC02-06CH11357.
