%% Guaranteed Automatic Integration Library (GAIL) 2.1 Documentation
%
% GAIL (Guaranteed Automatic Integration Library) is created, developed, 
% and maintained by Fred Hickernell (Illinois Institute of Technology), 
% Sou-Cheng Choi (NORC at the University of Chicago and IIT), 
% Yuhan Ding (IIT), Lan Jiang (IIT), Lluis Antoni Jimenez Rugama (IIT), Xin
% Tong (University of Illinois at Chicago), Yizhi Zhang (IIT), and Xuan 
% Zhou (IIT). 
%
% GAIL is a suite of algorithms for integration problems in one, many, and 
% infinite dimensions, and whose answers are guaranteed to be correct.
%
%% Introduction
%
% <html>
% <a href="help_license.html">GAIL License</a>
% <a href="help_readme.html">README</a>
% </html>
%
%% Functions
%
% <html>
% <a href="help_funappx_g.html">funappx_g</a>
% <a href="help_integral_g.html">integral_g</a>
% <a href="help_funmin_g.html">funmin_g</a>
% <a href="help_meanMC_g.html">meanMC_g</a>
% <a href="help_cubMC_g.html">cubMC_g</a>
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
% <a href="help_cubSobol_g.html">cubSobol_g</a>
% <a href="help_cubLattice_g.html">cubLattice_g</a>
% </html>
%
%
%% Website
% For more information about GAIL, visit
% <http://code.google.com/p/gail/ Gailteam>
%
%% GAIL License
%Copyright © 2015, Illinois Institute of Technology. All rights reserved.
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
% ANY USE OF THE SOFTWARE SHALL BE AT USER?S OWN RISK. THE COPYRIGHT 
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
% AGREED THAT THE MAXIMUM LIABILITY OF THE COPYRITGHT HOLDER AND 
% CONTRIBUTORS SHALL NOT EXCEED US$100.00.
%  
% THE USER OF THE SOFTWARE ACKNOWLEDGES THAT THE SOFTWARE IS BEING PROVIDED 
% WITHOUT CHARGE, AND AS A RESULT, THE USER, ACKNOWLEDGING THAT HE OR SHE 
% HAS READ THE SAME, AGREES THAT THE FOREGOING LIMITATIONS AND RESTRICTIONS 
% REPRESENT A REASONABLE ALLOCATION OF RISK.%% Guaranteed Automatic Integration Library (GAIL)
% GAIL Version 2.1, Mar 14, 2015.
% See LICENSE.m for copyright and disclaimer.
% 
% GAIL is a suite of algorithms for integration problems in one and many
% dimensions, and whose answers are guaranteed to be correct.
% 
% 
%% Developed by
% 
% Fred Hickernell, Sou-Cheng Choi, and their collaborators including
% Yuhan Ding, Lan Jiang, Lluis Antoni Jimenez Rugama, Yizhi 
% Zhang and Xuan Zhou, Department of Applied Mathematics, Illinois 
% Institute of Technology (IIT) and  Xin Tong, Department of Mathematics, 
% Statistics, and Computer Science, University of Illinois at Chicago. 
% We thank the contributions of Xincheng Sheng, and the IIT class of Math 
% 573 Reliable Mathematical Software, Fall 2013.
% 
% 
% Please cite the following software, papers, and materials:
% 
% 
% Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang, Lluis
% Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang, and Xuan Zhou, GAIL: 
% Guaranteed Automatic Integration Library (Version 2.1) [MATLAB Software],
% 2015. Available from http://code.google.com/p/gail/
% (this software)
% 
% Sou-Cheng T. Choi, "MINRES-QLP Pack and Reliable Reproducible Research via 
% Supportable Scientific Software," Journal of Open Research Software, Volume 2, 
% Number 1, e22, pp. 1-7, 2014.
% (describes principles of Reliable Reproducible Research and Supportable 
% Scientific Software)
% 
% Sou-Cheng T. Choi and Fred J. Hickernell, "IIT MATH-573 Reliable Mathematical
% Software" [Course Slides], Illinois Institute of Technology, Chicago, IL, 2013.
% Available from http://code.google.com/p/gail/
% (develops practices of Reliable Reproducible Research and Supportable 
% Scientific Software)
% 
% Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, "Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1)," Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
% (discusses practice and challenges in Sustainable Software for Science)
% 
% Nicholas Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, "The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls," Journal of Complexity 30, pp. 21-45, 2014.
% (describes integral_g.m and funappx_g.m)
% 
% Yuhan Ding, Fred J. Hickernell, and Sou-Cheng T. Choi, "Locally
% Adaptive Method for Approximating Univariate Functions in Cones with a
% Guarantee for Accuracy," working, 2015.
% (describes funappx_g.m)
% 
% Fred J. Hickernell, Lan Jiang, Yuewei Liu, and Art B. Owen,
% "Guaranteed conservative fixed width confidence intervals via Monte
% Carlo sampling," Monte Carlo and Quasi-Monte Carlo Methods 2012
% (J. Dick, F. Y. Kuo, G. W. Peters, and I. H. Sloan, eds.),
% Springer-Verlag, Berlin, pp. 105-128, 2014.
% (describes meanMC_g.m and cubMC_g.m)
% 
% Fred J. Hickernell and Lluis Antoni Jimenez Rugama, "Reliable Adaptive 
% Cubature Using Digital Sequences," submitted for publication, 2014.
% (describes cubSobol_g.m)
% 
% Lan Jiang and Fred J. Hickernell, "Guaranteed Conservative Confidence Intervals 
% for Means of Bernoulli Random Variables," submitted for publication, 2014.
% (describes meanMCBer_g)
% 
% Lluis Antoni Jimenez Rugama and Fred J. Hickernell, "Adaptive Multidimensional 
% Integration Based on Rank-1 Lattices," submitted for publication, 2014.
% (describes cubLattice_g.m)
% 
% Xin Tong, "A Guaranteed, Adaptive, Automatic Algorithm for Univariate 
% Function Minimization," MS thesis, Illinois Institute of Technology, 2014. 
% (describes funmin_g.m)
% 
%% Downloads
% 
% GAIL can be downloaded from http://code.google.com/p/gail/.
% 
% Alternatively, you can get a local copy of the GAIL repository with
% this command:
% 
%   git clone https://github.com/GailGithub/GAIL_Dev.git
% 
% 
%% Requirements
% 
% You will need to install MATLAB 7 or a later version.
% 
% 
%% Documentation
% 
% Detailed documentation is available at GAIL_Matlab/Documentation.
% 
% 
%% General Usage Notes
% 
% GAIL Version 2.1 includes the following eight algorithms:
% 
% 1.  funappx_g: One-dimensional function approximation on bounded interval
% 2.  integral_g:  One-dimensional integration on bounded interval
% 3.  meanMC_g:  Monte Carlo method for estimating mean of a random variable
% 4.  cubMC_g: Monte Carlo method for numerical multiple integration
% 5.  meanMCBer_g:  Monte Carlo method to estimate the mean of a Bernoulli random variable
% 6.  funmin_g: global minimum value of univariate function on a closed interval
% 7.  cubSobol_g: Quasi-Monte Carlo method using Sobol' cubature for a d-dimensional integration
% 8.  cubLattice_g: Quasi-Monte Carlo method using rank-1 Lattices cubature for a d-dimensional integration
% 
%% Installation Instruction
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
% 1.  Download DownloadInstallGail_2_1.m and put it where you want
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
% We provide quick doctests for each of the functions above. To run
% doctests in *funappx_g*, for example, issue the command *doctest
% funappx_g*.
%
% We also provide unit tests for MATLAB version 8 or later. To run unit
% tests for *funmin_g*, for instance, execute *run(ut_funmin_g)*.
%
%% Contact Information
% 
% Please send any queries, questions, or comments to
% gail-users@googlegroups.com or visit our project website:
% http://code.google.com/p/gail/
% 
% 
%% Acknowledgements
% 
% Our work was supported in part by grants from the National Science
% Foundation under grant NSF-DMS-1115392, and the Office of Advanced
% Scientific Computing Research, Office of Science, U.S. Department of
% Energy, under contract DE-AC02-06CH11357.
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
%% High dimension integration
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
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
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
%
%% funappx_g
% 1-D guaranteed locally adaptive function approximation (or function recovery) on [a,b]
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
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end point of interval, default value is 0
%
% * in_param.b --- right end point of interval, default value is 1
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6
%
% *Optional Input Arguments*
%
% * in_param.nlo --- lower bound of initial number of points we used,
%  default value is 10
%
% * in_param.nhi --- upper bound of initial number of points we used,
%  default value is 1000
%
% * in_param.nmax --- when number of points hits the value, iteration
%  will stop, default value is 1e7
%
% * in_param.maxiter --- max number of iterations, default value is 1000
%
% *Output Arguments*
%
% * fappx --- approximated function handle (Note: When Matlab version is
%  higher or equal to 8.3, fappx is an interpolant generated by
%  griddedInterpolant. When Matlab version is lower than 8.3, fappx is a
%  function handle generated by ppval and interp1.)
%
% * out_param.f --- input function
%
% * out_param.a --- left end point of interval
%
% * out_param.b --- right end point of interval
%
% * out_param.abstol --- guaranteed absolute error tolerance
%
% * out_param.nlo --- a lower bound of initial number of points we use
%
% * out_param.nhi --- an upper bound of initial number of points we use
%
% * out_param.nmax --- when number of points hits the value, iteration
%  will stop
%
% * out_param.maxiter --- max number of iterations
%
% * out_param.ninit --- initial number of points we use for each sub
%  interval
%
% <html>
% <ul type="square">
%  <li>out_param.exit --- this is a number defining the conditions of
%  success or failure satisfied when finishing the algorithm. The 
%  algorithm is considered successful (with out_param.exit == 0) if no 
%  other flags arise warning that the results are certainly not 
%  guaranteed. The initial value is 0 and the final value of this
%  parameter is encoded as follows:</li>
%   <ul type="circle">
%    <li>1   If reaching overbudget. It states whether
%                 the max budget is attained without reaching the
%                    guaranteed error tolerance.</li>
%    <li>2   If reaching overiteration. It states whether
%                    the max iterations is attained without reaching the
%                    guaranteed error tolerance.</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.iter --- number of iterations
%
% * out_param.npoints --- number of points we need to reach the
%  guaranteed absolute error tolerance
%
% * out_param.errest --- an estimation of the absolute error for the
%  approximation
%
% * out_param.nstar --- final value of the parameter defining the cone of
%  functions for which this algorithm is guaranteed for each
%  subinterval; nstar = ninit-2 initially
%
%% Guarantee
%
% For $[a,b]$, there exists a partition
%
% $$ P=\{[t_0,t_1], [t_1,t_2],  \ldots, [t_{L-1},t_L]\},  a=t_0 < t_1 < \cdots < t_L=b.$$
% 
% If the function to be approximated,  $f$ satisfies the cone condition
%
% $$\|f''\|_\infty \le \frac { 2\mathrm{nstar} }{t_l-t_{l-1} } \left\|f'-\frac{f(t_l)-f(t_{l-1})}{t_l-t_{l-1}}\right\|_\infty$$
% 
% for each sub interval $[t_{l-1},t_l]$, where $1 \le l \le L$, then the
% $fappx$ |output by this algorithm is guaranteed to satisfy
%
% $$\| f-fappx \|_{\infty} \le \mathrm{abstol}.$$
%
%
%% Examples
% *Example 1*

f = @(x) x.^2; [fappx, out_param] = funappx_g(f)

% Approximate function x^2 with default input parameter to make the error
% less than 1e-6.
%%
% *Example 2*

[fappx, out_param] = funappx_g(@(x) x.^2,0,100,1e-7,10,1000,1e8)

% Approximate function x^2 on [0,100] with error tolerance 1e-7, cost
% budget 10000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100

%%
% *Example 3*

clear in_param; in_param.a = -20; in_param.b = 20; in_param.nlo = 10;
in_param.nhi = 100; in_param.nmax = 1e8; in_param.abstol = 1e-7; 
[fappx, out_param] = funappx_g(@(x) x.^2, in_param)

% Approximate function x^2 on [-20,20] with error tolerance 1e-7, cost
% budget 1000000, lower bound of initial number of points 10 and upper
% bound of initial number of points 100
%%
% *Example 4*

clear in_param; f = @(x) x.^2;
[fappx, out_param] = funappx_g(f,'a',-10,'b',50,'nmax',1e6,'abstol',1e-7)

% Approximate function x^2 with error tolerance 1e-7, cost budget 1000000,
% lower bound of initial number of points 10 and upper
% bound of initial number of points 100
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
%% References
%
% [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, _The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls,_ Journal of Complexity 30, pp. 21-45, 2014.
%  
% [2]  Yuhan Ding, Fred J. Hickernell, and Sou-Cheng T. Choi, _Locally
% Adaptive Method for Approximating Univariate Functions in Cones with a
% Guarantee for Accuracy,_ working, 2015.
%          
% [3]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1) [MATLAB
% Software], 2015. Available from http://code.google.com/p/gail/
%
% [4] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [5] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [6] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% funmin_g
% 1-D guaranteed global minimum value on [a,b] and the subset
% containing optimal solutions
%% Syntax
% fmin = *funmin_g*(f)
%
% fmin = *funmin_g*(f,a,b,abstol,TolX)
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%
% fmin = *funmin_g*(f,in_param)
%
% [fmin, out_param] = *funmin_g*(f,...)
%% Description
%
% fmin = *funmin_g*(f) finds minimum value of function f on the default
%  interval [0,1] within the guaranteed absolute error tolerance of 1e-6
%  and the X tolerance of 1e-3. Default initial number of points is 100
%  and default cost budget is 1e7. Input f is a function handle.
%
% fmin = *funmin_g*(f,a,b,abstol,TolX) finds minimum value of
%  function f with ordered input parameters that define the finite
%  interval [a,b], a guaranteed absolute error tolerance abstol and a
%  guaranteed X tolerance TolX.
%
% fmin = *funmin_g*(f,'a',a,'b',b,'abstol',abstol,'TolX',TolX)
%  finds minimum value of function f on the interval [a,b] with a 
%  guaranteed absolute error tolerance abstol and a guaranteed X tolerance 
%  TolX. All five
%  field-value pairs are optional and can be supplied in different order.
%
% fmin = *funmin_g*(f,in_param) finds minimum value of function f on the
%  interval [in_param.a,in_param.b] with a guaranteed absolute error
%  tolerance in_param.abstol and a guaranteed X tolerance in_param.TolX.
%  If a field is not specified, the default value is used.
%
% [fmin, out_param] = *funmin_g*(f,...) returns minimum value fmin of
%  function f and an output structure out_param.
%
% *Input Arguments*
%
% * f --- input function
%
% * in_param.a --- left end point of interval, default value is 0
%
% * in_param.b --- right end point of interval, default value is 1
%
% * in_param.abstol --- guaranteed absolute error tolerance, default
%  value is 1e-6.
%
% * in_param.TolX --- guaranteed X tolerance, default value is 1e-3.
%
% *Optional Input Arguments*
%
% * in_param.nlo --- lower bound of initial number of points we used,
%  default value is 10
%
% * in_param.nhi --- upper bound of initial number of points we used,
%  default value is 1000
%
% * in_param.nmax --- cost budget, default value is 1e7.
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
% * out_param.TolX --- guaranteed X tolerance
%
% * out_param.nlo --- a lower bound of initial number of points we use
%
% * out_param.nhi --- an upper bound of initial number of points we use
%
% * out_param.nmax --- cost budget
%
% * out_param.ninit --- initial number of points we use
%
% * out_param.tau --- latest value of tau
%
% * out_param.npoints --- number of points needed to reach the guaranteed
%  absolute error tolerance or the guaranteed X tolerance
%
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0  Success</li>
%    <li>1  Number of points used is greater than out_param.nmax</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.errest --- estimation of the absolute error bound
%
% * out_param.volumeX --- the volume of intervals containing the point(s)
%  where the minimum occurs
%
% * out_param.tauchange --- it is 1 if out_param.tau changes, otherwise
%  it is 0
%
% * out_param.intervals --- the intervals containing point(s) where the
%  minimum occurs. Each column indicates one interval where the first
%  row is the left point and the second row is the right point.
%
%% Guarantee
%    
% If the function to be minimized, $f$ satisfies the cone condition
%
% $$\|f''\|_\infty \le  \frac {\tau}{b-a}\left\|f'-\frac{f(b)-f(a)}{b-a}
% \right\|_\infty,$$
%      
% then the $\mathrm{fmin}$ output by this algorithm is guaranteed to
% satisfy
%
% $$| \min f-\mathrm{fmin}| \le \mathrm{abstol},$$
%
% or
%
% $$\mathrm{volumeX} \le \mathrm{TolX},$$
%
% provided the flag $\mathrm{exitflag} = 0.$
%
%
%% Examples
% *Example 1*

f=@(x) (x-0.3).^2+1; [fmin,out_param] = funmin_g(f)

% Minimize function (x-0.3)^2+1 with default input parameter.
%%
% *Example 2*

f=@(x) (x-0.3).^2+1;
[fmin,out_param] = funmin_g(f,-2,2,1e-7,1e-4,10,10,1000000)

% Minimize function (x-0.3)^2+1 on [-2,2] with error tolerance 1e-4, X
% tolerance 1e-2, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 10

%%
% *Example 3*

clear in_param; in_param.a = -13; in_param.b = 8;
in_param.abstol = 1e-7; in_param.TolX = 1e-4;
in_param.nlo = 10; in_param.nhi = 100;
in_param.nmax = 10^6;
[fmin,out_param] = funmin_g(f,in_param)

% Minimize function (x-0.3)^2+1 on [-13,8] with error tolerance 1e-7, X
% tolerance 1e-4, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 100

%%
% *Example 4*

f=@(x) (x-0.3).^2+1;
[fmin,out_param] = funmin_g(f,'a',-2,'b',2,'nhi',100,'nlo',10,...
    'nmax',1e6,'abstol',1e-4,'TolX',1e-2)

% Minimize function (x-0.3)^2+1 on [-2,2] with error tolerance 1e-4, X
% tolerance 1e-2, cost budget 1000000, lower bound of initial number of
% points 10 and upper bound of initial number of points 100
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
% [1]  Xin Tong. A Guaranteed, _Adaptive, Automatic Algorithm for
% Univariate Function Minimization,_ MS thesis, Illinois Institute of 
% Technology, 2014.
% 
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1)
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% integral_g
% 1-D guaranteed function integration using trapezoidal rule
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
%  the interval [a,b] by trapezoidal rule with in a guaranteed absolute
%  error of 1e-6. Default starting number of sample points taken is 100
%  and default cost budget is 1e7. Input f is a function handle. The
%  function y = f(x) should accept a vector argument x and return a vector
%  result y, the integrand evaluated at each element of x.
%
% q = *integral_g*(f,a,b,abstol) computes q, the definite integral of
%  function f on the finite interval [a,b] by trapezoidal rule with the
%  ordered input parameters, and guaranteed absolute error tolerance
%  abstol.
%
% q = *integral_g*(f,'a',a,'b',b,'abstol',abstol) computes q, the definite
%  integral of function f on the finite interval [a,b] by trapezoidal rule
%  within a guaranteed absolute error tolerance abstol. All four
%  field-value pairs are optional and can be supplied.
%
% q = *integral_g*(f,in_param) computes q, the definite integral of
%  function f by trapezoidal rule within a guaranteed absolute error
%  in_param.abstol. If a field is not specified, the default value is
%  used.
%
% [q, out_param] = *integral_g*(f,...) returns the approximated 
%  integration q and output structure out_param.
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
% * out_param.tauchange --- it is true if the cone constant has been
%  changed, false otherwise. See [1] for details. If true, you may wish to
%  change the input in_param.ninit to a larger number.
% 
% * out_param.tauchange --- it is true if the cone constant has been
%  changed, false otherwise. See [1] for details. If true, you may wish to
%  change the input in_param.ninit to a larger number.
% 
% * out_param.iter --- number of iterations
%
% * out_param.npoints --- number of points we need to 
%  reach the guaranteed absolute error tolerance abstol.
%
% * out_param.errest --- approximation error defined as the differences
%  between the true value and the approximated value of the integral.
%
% * out_param.nstar --- final value of the parameter defining the cone of
%  functions for which this algorithm is guaranteed; nstar = ninit-2
%  initially and is increased as necessary
%
% <html>
% <ul type="square">
%  <li>out_param.exit --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0   Success</li> 
%    <li>1   Number of points used is greater than out_param.nmax</li> 
%    <li>2   Number of iterations is greater than out_param.maxiter</li> 
%   </ul>
% </ul>
% </html>
%
%% Guarantee
%    
% If the function to be integrated, $f$ satisfies the cone condition
%
% $$\|f''\|_1 \le \frac { \mathrm{nstar} }{2(b-a)}
% \left\|f'-\frac{f(b)-f(a)}{b-a}\right\|_1,$$
% 
% then the $q$ output by this algorithm is guaranteed to satisfy
%
% $$\left\| \int_{a}^{b} f(x) dx - q \right\|_{1} \le \mathrm{abstol},$$
%
% provided the flag $\mathrm{exceedbudget} = 0.$
%
% And the upper bound of the cost is
%
% $$\sqrt{ \frac{\mathrm{nstar}* (b-a)^2 \mathrm{Var}(f')}{2 \times \mathrm{abstol}}}
% + 2 \times \mathrm{nstar} +4.$$
%
%
%% Examples
% *Example 1*

f = @(x) x.^2; [q, out_param] = integral_g(f)

% Integrate function x with default input parameter to make the error less
% than 1e-7.

%%
% *Example 2*

[q, out_param] = integral_g(@(x) exp(-x.^2),'a',1,'b',2,...
   'nlo',100,'nhi',10000,'abstol',1e-5,'nmax',1e7)

% Integrate function x^2 with starting number of points 52, cost budget
% 10000000 and error tolerance 1e-8
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
% <a href="help_funappx_g.html">funappx_g</a>
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
%% References
%
% [1]  Nick Clancy, Yuhan Ding, Caleb Hamilton, Fred J. Hickernell, and
% Yizhi Zhang, _The Cost of Deterministic, Adaptive, Automatic Algorithms:
% Cones, Not Balls,_ Journal of Complexity 30, pp. 21-45, 2014.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1)
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% meanMC_g
% Monte Carlo method to estimate the mean of a random variable
%% Syntax
% tmu = *meanMC_g*(Yrand)
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha)
%
% tmu = *meanMC_g*(Yrand,'abstol',abstol,'reltol',reltol,'alpha',alpha)
%
% [tmu, out_param] = *meanMC_g*(Yrand,in_param)
%% Description
%
% tmu = *meanMC_g*(Yrand) estimates the mean, mu, of a random variable Y to
%  within a specified generalized error tolerance, 
%  tolfun:=max(abstol,reltol*| mu |), i.e., | mu - tmu | <= tolfun with
%  probability at least 1-alpha, where abstol is the absolute error
%  tolerance, and reltol is the relative error tolerance. Usually the
%  reltol determines the accuracy of the estimation, however, if the | mu |
%  is rather small, the abstol determines the accuracy of the estimation.
%  The default values are abstol=1e-2, reltol=1e-1, and alpha=1%. Input
%  Yrand is a function handle that accepts a positive integer input n and
%  returns an n x 1 vector of IID instances of the random variable Y.
%
% tmu = *meanMC_g*(Yrand,abstol,reltol,alpha) estimates the mean of a
%  random variable Y to within a specified generalized error tolerance
%  tolfun with guaranteed confidence level 1-alpha using all ordered
%  parsing inputs abstol, reltol, alpha.
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
% * Yrand --- the function for generating n IID instances of a random
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
%  percentage. default value is 1%.
%
% *Optional Input Arguments*
%
% * in_param.fudge --- standard deviation inflation factor, which should
%  be larger than 1, default value is 1.2.
%
% * in_param.nSig --- initial sample size for estimating the sample
%  variance, which should be a moderate large integer at least 30, the
%  default value is 1e4.
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
%  <li>out_param.exit --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0   Success</li>
%    <li>1   Not enough samples to estimate the mean</li>
%   </ul>
% </ul>
% </html>
%
% * out_param.kurtmax --- the upper bound on modified kurtosis.
%
% * out_param.time --- the time elapsed in seconds.
%
% <html>
% <ul type="square">
%  <li>out_param.flag --- parameter checking status:</li>
%   <ul type="circle">
%    <li>1   checked by meanMC_g</li>
%   </ul>
% </ul>
% </html>
%
%%  Guarantee
% This algorithm attempts to calculate the mean, mu, of a random variable
% to a prescribed error tolerance, tolfun:= max(abstol,reltol*| mu |), with
% guaranteed confidence level 1-alpha. If the algorithm terminated without
% showing any warning messages and provide an answer tmu, then the follow
% inequality would be satisfied:
% 
% Pr(| mu - tmu | <= tolfun) >= 1-alpha
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

% Calculate the mean of x^2 when x is uniformly distributed in
% [0 1], with the absolute error tolerance = 1e-3 and uncertainty 5%.

  in_param.reltol=0; in_param.abstol = 1e-3;in_param.reltol = 0;
  in_param.alpha = 0.05; Yrand=@(n) rand(n,1).^2;
  tmu = meanMC_g(Yrand,in_param)

%%
% *Example 2*

% Calculate the mean of exp(x) when x is uniformly distributed in
% [0 1], with the absolute error tolerance 1e-3.

  tmu = meanMC_g(@(n)exp(rand(n,1)),1e-3,0)

%%
% *Example 3*

% Calculate the mean of cos(x) when x is uniformly distributed in
% [0 1], with the relative error tolerance 1e-2 and uncertainty 0.05.

  tmu = meanMC_g(@(n)cos(rand(n,1)),'reltol',1e-2,'abstol',0,...
      'alpha',0.05)
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
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
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
%% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, _Guaranteed
% conservative fixed width confidence intervals via Monte Carlo
% sampling,_ Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
% Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin,
% 2014. arXiv:1208.4318 [math.ST]
%          
% [2]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1) [MATLAB
% Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% meanMCBer_g
% Monte Carlo method to estimate the mean of a Bernoulli random
% variable to within a specified absolute error tolerance with guaranteed
% confidence level 1-alpha.
%% Syntax
% pHat = *meanMCBer_g*(Yrand)
%
% pHat = *meanMCBer_g*(Yrand,abstol,alpha,nmax)
%
% pHat = *meanMCBer_g*(Yrand,'abstol',abstol,'alpha',alpha,'nmax',nmax)
%
% [pHat, out_param] = *meanMCBer_g*(Yrand,in_param)
%% Description
%
% pHat = *meanMCBer_g*(Yrand) estimates the mean of a Bernoulli random
%  variable Y to within a specified absolute error tolerance with
%  guaranteed confidence level 99%. Input Yrand is a function handle that
%  accepts a positive integer input n and returns a n x 1 vector of IID
%  instances of the Bernoulli random variable Y.
% 
% pHat = *meanMCBer_g*(Yrand,abstol,alpha,nmax) estimates the mean
%  of a Bernoulli random variable Y to within a specified absolute error
%  tolerance with guaranteed confidence level 1-alpha using all ordered
%  parsing inputs abstol, alpha and nmax.
% 
% pHat = *meanMCBer_g*(Yrand,'abstol',abstol,'alpha',alpha,'nmax',nmax)
%  estimates the mean of a Bernoulli random variable Y to within a
%  specified absolute error tolerance with guaranteed confidence level
%  1-alpha. All the field-value pairs are optional and can be supplied in
%  different order.
% 
% [pHat, out_param] = *meanMCBer_g*(Yrand,in_param) estimates the mean
%  of a Bernoulli random variable Y to within a specified absolute error
%  tolerance with the given parameters in_param and produce the estimated
%  mean pHat and output parameters out_param.
% 
% *Input Arguments*
%
% * Yrand --- the function for generating IID instances of a Bernoulli
%            random variable Y whose mean we want to estimate.
%
% * pHat --- the estimated mean of Y.
%
% * in_param.abstol --- the absolute error tolerance, the default value is 1e-2.
% 
% * in_param.alpha --- the uncertainty, the default value is 1%.
% 
% * in_param.nmax --- the sample budget, the default value is 1e9.
% 
% *Output Arguments*
%
% * out_param.n --- the total sample used.
%
% * out_param.time --- the time elapsed in seconds.
% 
% <html>
% <ul type="square">
%  <li>out_param.exit --- the state of program when exiting:</li>
%   <ul type="circle">
%    <li>0   success</li>
%    <li>1   Not enough samples to estimate p with guarantee</li>
%   </ul>
% </ul>
% </html>                    
%
%%  Guarantee
%
% If the sample size is calculated according Hoeffding's inequality, which
% equals to ceil(log(2/out_param.alpha)/(2*out_param.abstol^2)), then the
% following inequality must be satisfied:
%
% Pr(| p - pHat | <= abstol) >= 1-alpha.
% 
% Here p is the true mean of Yrand, and pHat is the output of MEANMCBER_G.
%
% Also, the cost is deterministic.
%
%%   Examples
%%
% *Example 1*

% Calculate the mean of a Bernoulli random variable with true p=1/90,
% absolute error tolerance 1e-3 and uncertainty 0.01.
 
    in_param.abstol=1e-3; in_param.alpha = 0.01; in_param.nmax = 1e9; 
    p=1/9; Yrand=@(n) rand(n,1)<p;
    pHat = meanMCBer_g(Yrand,in_param)
 
%%
% *Example 2*

% Using the same function as example 1, with the absolute error tolerance
% 1e-4.

    pHat = meanMCBer_g(Yrand,1e-4)
    
%%
% *Example 3*

% Using the same function as example 1, with the absolute error tolerance
% 1e-2 and uncertainty 0.05.

    pHat = meanMCBer_g(Yrand,'abstol',1e-2,'alpha',0.05)
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
%% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, _Guaranteed
% conservative fixed width confidence intervals via Monte Carlo
% sampling,_ Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
% Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), Springer-Verlag, Berlin,
% 2014. arXiv:1208.4318 [math.ST]
%
% [2]  Lan Jiang and Fred J. Hickernell, _Guaranteed Conservative
% Confidence Intervals for Means of Bernoulli Random Variables,_
% submitted for publication, 2014.
%          
% [3]  Sou-Cheng T. Choi, Yuhan Ding, Fred J. Hickernell, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1) [MATLAB
% Software], 2015. Available from http://code.google.com/p/gail/
%
% [4] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [5] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [6] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
%% cubMC_g
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
%  max(abstol, reltol*| I |), i.e., | I - Q | <= tolfun with probability at
%  least 1-alpha, where abstol is the absolute error tolerance, and reltol
%  is the relative error tolerance. Usually the reltol determines the
%  accuracy of the estimation, however, if the | I | is rather small, the
%  abstol determines the accuracy of the estimation. The default values
%  are abstol=1e-2, reltol=1e-1, and alpha=1%. Input f is a function
%  handle that accepts an n x d matrix input, where d is the dimension of
%  the hyperbox, and n is the number of points being evaluated
%  simultaneously. The input hyperbox is a 2 x d matrix, where the first
%  row corresponds to the lower limits and the second row corresponds to
%  the upper limits.
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
%  the default is 'uniform'. The other measure could be handled is
%  'normal'/'Gaussian'. The input should be a string type, hence with
%  quotes.
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
% * out_param.ntot --- total sample used.
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
%    <li>11  hyperbox is not 2 x d</li>
%    <li>12  hyperbox is only a point in one direction</li>
%    <li>13  hyperbox is infinite when measure is 'uniform'</li>
%    <li>14  hyperbox is not doubly infinite when measure
%                        is 'normal'</li>
%   </ul>
% </ul>
% </html>
% 
%%  Guarantee
% This algorithm attempts to calculate the integral of function f over a
% hyperbox to a prescribed error tolerance tolfun:= max(abstol,reltol*| I |)
% with guaranteed confidence level 1-alpha. If the algorithm terminated
% without showing any warning messages and provide an answer Q, then the
% follow inequality would be satisfied:
% 
% Pr(| Q - I | <= tolfun) >= 1-alpha
%
% The cost of the algorithm, N_tot, is also bounded above by N_up, which is
% a function in terms of abstol, reltol, nSig, n1, fudge, kurtmax, beta.
% And the following inequality holds:
% 
% Pr (N_tot <= N_up) >= 1-beta
%
% Please refer to our paper for detailed arguments and proofs.
%
%%  Examples
% *Example 1*

% Estimate the integral with integrand f(x) = sin(x) over the interval
% [1;2]
% 

 f = @(x) sin(x); interval = [1;2];
 Q = cubMC_g(f,interval,'uniform',1e-3,1e-2)
 
%% 
% *Example 2*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) over the
% hyperbox [0 0;1 1], where x is a vector x = [x1 x2].
% 

 f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [0 0;1 1];
 Q = cubMC_g(f,hyperbox,'measure','uniform','abstol',1e-3,...
     'reltol',1e-13)

%%
% *Example 3*

% Estimate the integral with integrand f(x) = 2^d*prod(x1*x2*...*xd) +
% 0.555 over the hyperbox [zeros(1,d);ones(1,d)], where x is a vector x =
% [x1 x2... xd].


  d = 3;f = @(x) 2^d*prod(x,2)+0.555; hyperbox = [zeros(1,d);ones(1,d)];
  in_param.abstol = 1e-3; in_param.reltol=1e-3;
  Q = cubMC_g(f,hyperbox,in_param)

%%
% *Example 4* 

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% hyperbox [-inf -inf;inf inf], where x is a vector x = [x1 x2].


 f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-inf -inf;inf inf];
 Q = cubMC_g(f,hyperbox,'normal',0,1e-2)
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
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
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
%% References
%
% [1]  F. J. Hickernell, L. Jiang, Y. Liu, and A. B. Owen, _Guaranteed
% conservative fixed width confidence intervals via Monte Carlo
% sampling,_ Monte Carlo and Quasi-Monte Carlo Methods 2012 (J. Dick, F.
% Y. Kuo, G. W. Peters, and I. H. Sloan, eds.), pp. 105-128,
% Springer-Verlag, Berlin, 2014. DOI: 10.1007/978-3-642-41095-6_5
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1)
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
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
%  guaranteed not to be greater than a specific generalized error tolerance,
%  tolfun:=max(abstol,reltol*| integral(f) |). Input f is a function handle. f should
%  accept an n x d matrix input, where d is the dimension and n is the 
%  number of points being evaluated simultaneously. The input hyperbox is
%  a 2 x d matrix, where the first row corresponds to the lower limits 
%  and the second row corresponds to the upper limits of the integral.
%  Given the construction of our Lattices, d must be a positive integer
%  with 1<=d<=250.
% 
% q = *cubLattice_g*(f,hyperbox,measure,abstol,reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either. Inputs f and hyperbox 
%  are required. The other optional inputs are in the correct order:
%  measure,abstol,reltol,shift,mmin,mmax,fudge,transform,toltype and
%  theta.
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
%  greater than 250. By default f is f=@ x.^2.
%
% * hyperbox --- the integration region defined by its bounds. It must be
%  a 2 x d matrix, where the first row corresponds to the lower limits 
%  and the second row corresponds to the upper limits of the integral.
%  The default value is [0;1].
%
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in the hyperbox
%  or normally distributed with covariance matrix I_d. The only possible
%  values are 'uniform' or 'normal'. For 'uniform', the hyperbox must be
%  a finite volume while for 'normal', the hyperbox can only be defined as 
%  (-Inf,Inf)^d. By default it is 'uniform'.
%
% * in_param.abstol --- the absolute error tolerance, abstol>=0. By 
%  default it is 1e-4.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  in [0,1]. Default value is 1e-2.
% 
% *Optional Input Arguments*
%
% * in_param.shift --- the Rank-1 lattices can be shifted to avoid the
%  origin or other particular points. By default we consider a uniformly
%  [0,1) random shift.
% 
% * in_param.mmin --- the minimum number of points to start is 2^mmin.
%  The cone condition on the Fourier coefficients decay requires a
%  minimum number of points to start. The advice is to consider at least
%  mmin=10. mmin needs to be a positive integer with mmin<=mmax. By
%  default it is 10.
% 
% * in_param.mmax --- the maximum budget is 2^mmax. By construction of
%  our Lattices generator, mmax is a positive integer such that
%  mmin<=mmax<=26. The default value is 24.
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
%  periodic functions. If the input function f is not, there are 5
%  types of transform to periodize it without modifying the result. 
%  By default it is the Baker's transform. The options are:</li>
%   <ul type="circle">
%    <li>id : no transformation.</li>
%    <li>Baker : Baker's transform or tent map in each coordinate. Preserving
%              only continuity but simple to compute. Chosen by default.</li>
%    <li>C0 : polynomial transformation only preserving continuity.</li>
%    <li>C1 : polynomial transformation preserving the first derivative.</li>
%    <li>C1sin : Sidi's transform with sine, preserving the first derivative.
%              This is in general a better option than 'C1'.</li>
%   </ul>
%  </ul>
% </html>
%
% * in_param.toltype --- this is the generalized tolerance function.
%  There are two choices, 'max' which takes
%  max(abstol,reltol*| integral(f) | ) and 'comb' which is the linear combination
%  theta*abstol+(1-theta)*reltol*| integral(f) | . Theta is another 
%  parameter to be specified with 'comb'(see below). For pure absolute
%  error, either choose 'max' and set reltol = 0 or choose 'comb' and set
%  theta = 1. For pure relative error, either choose 'max' and set 
%  abstol = 0 or choose 'comb' and set theta = 0. Note that with 'max',
%  the user can not input abstol = reltol = 0 and with 'comb', if theta = 1
%  abstol con not be 0 while if theta = 0, reltol can not be 0.
%  By default toltype is 'max'.
% 
% * in_param.theta --- this input is parametrizing the toltype 
%  'comb'. Thus, it is only active when the toltype
%  chosen is 'comb'. It establishes the linear combination weight
%  between the absolute and relative tolerances
%  theta*abstol+(1-theta)*reltol*| integral(f) |. Note that for theta = 1, 
%  we have pure absolute tolerance while for theta = 0, we have pure 
%  relative tolerance. By default, theta=1.
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
% * out_param.bound_err --- predicted bound on the error based on the cone
%  condition. If the function lies in the cone, the real error will be
%  smaller than generalized tolerance.
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
%                    <li>1 : If reaching overbudget. It states whether
%                    the max budget is attained without reaching the
%                    guaranteed error tolerance.</li> 
%                    <li>2 : If the function lies outside the cone. In
%                    this case, results are not guaranteed. Note that
%                    this parameter is computed on the transformed
%                    function, not the input function. For more
%                    information on the transforms, check the input
%                    parameter in_param.transform; for information about
%                    the cone definition, check the article mentioned
%                    below.</li>
%   </ul>
%  </ul>
% </html>
%
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in
% dimension d with a prescribed generalized error tolerance. The Fourier
% coefficients of the integrand are assumed to be absolutely convergent. If
% the algorithm terminates without warning messages, the output is given
% with guarantees under the assumption that the integrand lies inside a
% cone of functions. The guarantee is based on the decay rate of the
% Fourier coefficients. For more details on how the cone is defined, please
% refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval
% [0,1)^2:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)];
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','C1sin')

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
  q = cubLattice_g(f,hyperbox,'normal',1e-3,1e-3,'transform','C1sin')

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2)^2:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
  q = cubLattice_g(f,hyperbox,'uniform',1e-3,1e-2,'transform','C1')

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
  q = cubLattice_g(f,hyperbox,'normal',1e-4,1e-2,'transform','C1sin')

%%
% *Example 5*

% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the
% interval [0,1)^5 with pure absolute error 1e-5.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0)

%%
% *Example 6*

% Estimate the integral with integrand f(x) = 3./(5-4*(cos(2*pi*x))) in the interval
% [0,1) with pure absolute error 1e-5.

  f = @(x) 3./(5-4*(cos(2*pi*x))); hyperbox = [0;1];
  q = cubLattice_g(f,hyperbox,'uniform',1e-5,0,'transform','id')
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
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%% References
%
% [1] Lluis Antoni Jimenez Rugama and Fred J. Hickernell, _Adaptive
% Multidimensional Integration Based on Rank-1 Lattices,_ 2014. Submitted
% for publication: arXiv:1411.1966.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1)
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
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
%  number of points being evaluated simultaneously. The input hyperbox is
%  a 2 x d matrix, where the first row corresponds to the lower limits 
%  and the second row corresponds to the upper limits of the integral.
%  Given the construction of Sobol' sequences, d must be a positive 
%  integer with 1<=d<=1111.
%
% q = *cubSobol_g*(f,hyperbox,measure,abstol,reltol)
%  estimates the integral of f over the hyperbox. The answer
%  is given within the generalized error tolerance tolfun. All parameters
%  should be input in the order specified above. If an input is not specified,
%  the default value is used. Note that if an input is not specified,
%  the remaining tail cannot be specified either. Inputs f and hyperbox 
%  are required. The other optional inputs are in the correct order:
%  measure,abstol,reltol,mmin,mmax,fudge,toltype and
%  theta.
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
% * f --- the integrand whose input should be a matrix n x d where n is
%  the number of data points and d the dimension, which cannot be
%  greater than 1111. By default f is f=@ x.^2.
%
% * hyperbox --- the integration region defined by its bounds. It must be
%  a 2 x d matrix, where the first row corresponds to the lower limits 
%  and the second row corresponds to the upper limits of the integral.
%  The default value is [0;1].
%
% * in_param.measure --- for f(x)*mu(dx), we can define mu(dx) to be the
%  measure of a uniformly distributed random variable in the hyperbox
%  or normally distributed with covariance matrix I_d. The only possible
%  values are 'uniform' or 'normal'. For 'uniform', the hyperbox must be
%  a finite volume while for 'normal', the hyperbox can only be defined as 
%  (-Inf,Inf)^d. By default it is 'uniform'.
%
% * in_param.abstol --- the absolute error tolerance, abstol>=0. By 
%  default it is 1e-4.
%
% * in_param.reltol --- the relative error tolerance, which should be
%  in [0,1]. Default value is 1e-2.
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
% * in_param.toltype --- this is the generalized tolerance function.
%  There are two choices, 'max' which takes
%  max(abstol,reltol*| integral(f) | ) and 'comb' which is the linear combination
%  theta*abstol+(1-theta)*reltol*| integral(f) | . Theta is another 
%  parameter to be specified with 'comb'(see below). For pure absolute
%  error, either choose 'max' and set reltol = 0 or choose 'comb' and set
%  theta = 1. For pure relative error, either choose 'max' and set 
%  abstol = 0 or choose 'comb' and set theta = 0. Note that with 'max',
%  the user can not input abstol = reltol = 0 and with 'comb', if theta = 1
%  abstol con not be 0 while if theta = 0, reltol can not be 0.
%  By default toltype is 'max'.
% 
% * in_param.theta --- this input is parametrizing the toltype 
%  'comb'. Thus, it is only active when the toltype
%  chosen is 'comb'. It establishes the linear combination weight
%  between the absolute and relative tolerances
%  theta*abstol+(1-theta)*reltol*| integral(f) |. Note that for theta = 1, 
%  we have pure absolute tolerance while for theta = 0, we have pure 
%  relative tolerance. By default, theta=1.
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
% <html>
% <ul type="square">
%  <li>out_param.exitflag --- this is a binary vector stating whether
%  warning flags arise. These flags tell about which conditions make the
%  final result certainly not guaranteed. One flag is considered arisen
%  when its value is 1. The following list explains the flags in the
%  respective vector order:</li>
%   <ul type="circle">
%                    <li>1    If reaching overbudget. It states whether
%                    the max budget is attained without reaching the
%                    guaranteed error tolerance.</li>
%                    <li>2   If the function lies outside the cone. In
%                    this case, results are not guaranteed. For more
%                    information about the cone definition, check the
%                    article mentioned below.</li>
%   </ul>
%  </ul>
% </html>
%
%%  Guarantee
%
% This algorithm computes the integral of real valued functions in
% dimension d with a prescribed generalized error tolerance. The
% Walsh-Fourier coefficients of the integrand are assumed to be absolutely
% convergent. If the algorithm terminates without warning messages, the
% output is given with guarantees under the assumption that the integrand
% lies inside a cone of functions. The guarantee is based on the decay rate
% of the Walsh-Fourier coefficients. For more details on how the cone is
% defined, please refer to the references below.
%
%% Examples
%
%%
% *Example 1*

% Estimate the integral with integrand f(x) = x1.*x2 in the interval
% [0,1)^2:

  f = @(x) prod(x,2); hyperbox = [zeros(1,2);ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0)

%%
% *Example 2*

% Estimate the integral with integrand f(x) = x1.^2.*x2.^2.*x3.^2
% in the interval R^3 where x1, x2 and x3 are normally distributed:

  f = @(x) x(:,1).^2.*x(:,2).^2.*x(:,3).^2; hyperbox = [-inf(1,3);inf(1,3)];
  q = cubSobol_g(f,hyperbox,'normal',1e-3,1e-3)

%%
% *Example 3*

% Estimate the integral with integrand f(x) = exp(-x1^2-x2^2) in the
% interval [-1,2)^2:

  f = @(x) exp(-x(:,1).^2-x(:,2).^2); hyperbox = [-ones(1,2);2*ones(1,2)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-3,1e-2)

%%
% *Example 4*

% Estimate the price of an European call with S0=100, K=100, r=sigma^2/2,
% sigma=0.05 and T=1.

  f = @(x) exp(-0.05^2/2)*max(100*exp(0.05*x)-100,0); hyperbox = [-inf(1,1);inf(1,1)];
  q = cubSobol_g(f,hyperbox,'normal',1e-4,1e-2)

%%
% *Example 5*

% Estimate the integral with integrand f(x) = 8*x1.*x2.*x3.*x4.*x5 in the
% interval [0,1)^5 with pure absolute error 1e-5.

  f = @(x) 8*prod(x,2); hyperbox = [zeros(1,5);ones(1,5)];
  q = cubSobol_g(f,hyperbox,'uniform',1e-5,0)
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
% <a href="help_meanMCBer_g.html">meanMCBer_g</a>
% </html>
%
% <html>
% <a href="help_integral_g.html">integral_g</a>
% </html>
%
%% References
%
% [1] Fred J. Hickernell and Lluis Antoni Jimenez Rugama, _Reliable adaptive
% cubature using digital sequences,_ 2014. Submitted for publication:
% arXiv:1410.8615.
%
% [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
% Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
% GAIL: Guaranteed Automatic Integration Library (Version 2.1)
% [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
% [3] Sou-Cheng T. Choi, _MINRES-QLP Pack and Reliable Reproducible
% Research via Supportable Scientific Software,_ Journal of Open Research
% Software, Volume 2, Number 1, e22, pp. 1-7, 2014.
%
% [4] Sou-Cheng T. Choi and Fred J. Hickernell, _IIT MATH-573 Reliable
% Mathematical Software_ [Course Slides], Illinois Institute of
% Technology, Chicago, IL, 2013. Available from
% http://code.google.com/p/gail/ 
%
% [5] Daniel S. Katz, Sou-Cheng T. Choi, Hilmar Lapp, Ketan Maheshwari,
% Frank Loffler, Matthew Turk, Marcus D. Hanwell, Nancy Wilkins-Diehr,
% James Hetherington, James Howison, Shel Swenson, Gabrielle D. Allen,
% Anne C. Elster, Bruce Berriman, Colin Venters, _Summary of the First
% Workshop On Sustainable Software for Science: Practice And Experiences
% (WSSSPE1),_ Journal of Open Research Software, Volume 2, Number 1, e6,
% pp. 1-21, 2014.
%
% If you find GAIL helpful in your work, please support us by citing the
% above papers, software, and materials.
%
