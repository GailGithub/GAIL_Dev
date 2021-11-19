---

<!-- 
1. A list of the authors of the software and their affiliations, using the correct format (see the example below).
2. A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
3. A Statement of Need section that clearly illustrates the research purpose of the software.
4. A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
5. Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
6. Acknowledgement of any financial support. 
-->

title: 'Guaranteed Automatic Integration Library (GAIL)'
tags:
  - Univariate Function Approximation
  - Multivariate Integration
  - Univariate Optimization
  - Adaptive
  - (Quasi-)Monte Carlo
  - Bayesian Cubature
authors:
  - name: Xin Tong^[corresponding author, co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-4718-1198
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Sou-Cheng T. Choi^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 
    affiliation: "2, 3"
  - name: Yuhan Ding^[co-first author]
    orcid: 
    affiliation: 2
  - name: Fred J. Hickernell^[co-first author]
    orcid: 
    affiliation: 2
  - name: Lan Jiang^[co-first author]
    orcid: 
    affiliation: 4
  - name: Lluı`s Antoni Jime`nez Rugama^[co-first author]
    orcid: 
    affiliation: 5
  - name: Jagadeeswaran Rathinavel^[co-first author]
    orcid: 
    affiliation: 6
  - name: Kan Zhang^[co-first author]
    orcid: 
    affiliation: 2
  - name: Yizhi Zhang^[co-first author]
    orcid: 
    affiliation: 7
  - name: Xuan Zhou^[co-first author]
    orcid: 
    affiliation: 8


affiliations:
 - name: Department of Scientific Computing, Florida State University
   index: 1
 - name: Department of Applied Mathematics, Illinois Institute of Technology
   index: 2
 - name: Kamakura Corporation
   index: 3
 - name: Compass Inc.
   index: 4
 - name: UBS Financial Services Inc.
   index: 5
 - name: Wi-Tronix LLC
   index: 6
 - name: Jamran International Inc.
   index: 7
 - name: J.P. Morgan.
   index: 8
date: xx xxx 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 
aas-journal: 
---

# Summary

Function approximation, integration, and optimization are three fundamental mathematical problems. They are especially challenging when the functions involved fluctuate wildly in certain parts of the domain, or if the domain is high dimensional. Ideally, these algorithms should possess a rigorous mathematical framework, databased (probabilistic) error bounds, and advanced sampling strategies for efficiency.

The Guaranteed Automatic Integration Library (GAIL) is our multi-year research effort directed towards addressing the aforementioned challenge. Briefly, GAIL is a free, open-source MATLAB software library with nine main algorithms undergirded by over a dozen peer-reviewed publications. GAIL solves problems in univariate and multivariate integration, and in univariate function approximation and optimization. GAIL algorithms sample data points of a given function in an adaptive manner and automatically stop when the error tolerance has been reached. In some cases, GAIL algorithms are proven to have asymptotically optimal computational cost. We consistently employ good software development practices for GAIL such as unit tests, searchable online documentation, and git version control. 


# Statement of need
Function approximation, integration, and optimization are fundamental problems requiring numerical solutions. A fundamental question is how and when to stop the computation. Theoretical error bounds typically contain unknown quantities, such as the norm of the input function. This makes them impractical as stopping criteria.

Therefore, practical algorithms that adapt the computation to the error requirement are often based on heuristics. These include a popular adaptive quadrature algorithm of Shampine `[@Sha08a]` that is a part of MATLAB [28] and the Chebfun library [8]. Heuristics based on function data tend to lack theoretical support; one does not know when they work and when they do not. A warning against commonly used adaptive quadrature stopping criteria is given by Lyness [25].



# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

Hickernell, Ding, and Choi wish to thank students from the following IIT courses for discussion: SCI 498 Adaptive Monte Carlo Algorithms with Applications to Fi- nancial Risk Management, Summer 2016; MATH 491 Reading & Research, Summer 2015; SCI 498/MATH 491 Computational Social Sciences, Summer 2016; MATH 491-195 Solving Problems in the Social Sciences Using Tools from Computational Mathematics and Statistics, Summer 2015; Math 573/SCI 498 Reliable Mathemat- ical Software, Fall 2013 and Fall 2018.

# References

bibliography: paper.bib