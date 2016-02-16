function dt_funappxNoPenalty_g
%DT_funappxNoPenalty_g small doctest for funappxNoPenalty_g
%
%   Example 1:
%
%   >> funappxNoPenalty_g
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 2:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2, in_param)
%
%   out_param = 
% 
%              a: -10
%         abstol: 1.0000e-***8
%              b: 10
%              f: @(x)x.^2
%        maxiter: 1000
%            nhi: 1000
%            nlo: 10
%           nmax: 10000000
%          ninit: 804
%           exit: [2x1 logical]
%           iter: 9
%        npoints: 205569
%         errest: 2.8480e-***9
%              x: [1x205569 double]
%   >> out_param.bytes <= 25288738
%      1
%
% 
%   Example 3: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-8); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2, in_param)
%
%   out_param = 
% 
%           a: -5
%      abstol: 1.0000e-08
%           b: 5
%           f: @(x)x.^2
%     maxiter: 1000
%         nhi: 500
%         nlo: 10
%        nmax: 10000000
%       ninit: 351
%        exit: [2x1 logical]
%        iter: 9
%     npoints: 89601
%      errest: 3.7478e-***9
%           x: [1x89601 double]
%   >> out_param.bytes <= 11025042
%      1
%      
%
%   Example 4:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
%   out_param = 
% 
%                a: -2
%           abstol: 1.0000e-***7
%                b: 2
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            ninit: 42
%             exit: [2x1 logical]
%             iter: 9
%          npoints: 10497
%           errest: 4.3699e-***8
%                x: [1x10497 double]
%   >> out_param.bytes <= 1295406
%      1
%
%
%   Example 5:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
%   out_param = 
% 
%                a: -3
%           abstol: 1.0000e-***6
%                b: 0
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 40
%              nlo: 20
%             nmax: 10000000
%            ninit: 34
%             exit: [2x1 logical]
%             iter: 7
%          npoints: 2113
%           errest: 6.1248e-***7
%                x: [1x2113 double]
%   >> out_param.bytes <= 263930
%      1
%
% 
%   Example 6:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-2,5,1e-6,20,50)
%   
%   out_param = 
% 
%           a: -2
%      abstol: 1.0000e-***6
%           b: 5
%           f: @(x)x.^2
%     maxiter: 1000
%         nhi: 50
%         nlo: 20
%        nmax: 10000000
%       ninit: 45
%        exit: [2x1 logical]
%        iter: 8
%     npoints: 5633
%      errest: 4.6617e-***7
%           x: [1x5633 double]
%   >> out_param.bytes <= 696546
%      1
%
%
%   Example 7:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-3,3,1e-7,20,50)
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-***7
%                b: 3
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            ninit: 44
%             exit: [2x1 logical]
%             iter: 9
%          npoints: 11009
%           errest: 8.9388e-***8
%                x: [1x11009 double]
%   >> out_param.bytes <= 1357794
%      1
%
%
%   Example 8:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-5,10,1e-7,10,20)
%
%   out_param = 
% 
%                a: -5
%           abstol: 1.0000e-07
%                b: 10
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 20
%              nlo: 10
%             nmax: 10000000
%            ninit: 20
%             exit: [2x1 logical]
%             iter: 12
%          npoints: 38913
%           errest: 4.4596e-***8
%                x: [1x38913 double]
%   >> out_param.bytes <= 4789986
%      1
