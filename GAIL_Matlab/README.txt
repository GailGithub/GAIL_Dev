Guaranteed Automatic Integration Library (GAIL)
GAIL Version 1.0.0, Sept 3, 2013.
See LICENSE.m for copyright and disclaimer.

GAIL is a suite of algorithms for integration problems in one and many
dimensions, and whose answers are guaranteed to be correct.


Developed by
-------------

Fred Hickernell (Illinois Institute of Technology), Sou-Cheng Choi
(University of Chicago, Argonne National Laboratory, and IIT), and
their collaborators including Yuhan Ding (IIT), Lan Jiang (IIT), and
Yizhi Zhang (IIT), Department of Applied Mathematics, Illinois
Institute of Technology.


Please cite the following software and papers:


Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, and
Yizhi Zhang. "GAIL: Guaranteed Automatic Integration Library (Version
1)" [MATLAB Software], 2013.  Available from
http://code.google.com/p/gail/


Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
"Guaranteed conservative fixed width confidence intervals via Monte
Carlo sampling", Monte Carlo and Quasi-Monte Carlo Methods 2012
(J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
Springer-Verlag, Berlin, 2014, to appear, arXiv:1208.4318 [math.ST]


Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
Yizhi Zhang, "The complexity of guaranteed automatic algorithms:
Cones, not balls", preprint, 2013 arXiv:1303.2412 [math.ST].


Downloads
------------
GAIL can be downloaded from http://code.google.com/p/gail/.

Alternatively, you can Get a local copy of the gail repository with this command:

git clone https://code.google.com/p/gail/


Requirements
------------

You will need to install MATLAB 7 or a later version.


Documentation
-------------

Detailed documentation is available at GAIL_Matlab/Documentation.


General Usage Notes
-------------------

GAIL Version 1.0.0 includes the following three algorithms:

1.  funappx_g: One-dimensional function approximation on unit interval
2.  integral_g:  One-dimensional integration on unit interval
3.  meanMC_g:  Monte Carlo method for estimating mean of a random variable


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

3.  To check if you have installed GAIL successfully, type: help
    funappx_g to see if it shows up.


To uninstall GAIL, execute "GAIL_Uninstall".


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