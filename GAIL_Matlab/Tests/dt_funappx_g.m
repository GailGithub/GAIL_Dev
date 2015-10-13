function dt_funappx_g
%DT_FUNAPPX_G small doctest for funappx_g
%
%   Example 1:
%
%   >> funappx_g
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 2:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [~, out_param] = funappx_g(@(x) x.^2, in_param)
%
% out_param = 
% 
%              a: -10
%         abstol: 1.0000e-***8
%              b: 10
%              f: @(x)x.^2
%        maxiter: 1000
%            nhi: 1000
%            nlo: 10
%           nmax: 10000000
%          nstar: [1x256 double]
%          ninit: 1609
%           exit: [2x1 logical]
%           iter: 9
%        npoints: 411649
%         errest: 8.3292e-***9
%              x: [1x411649 double]
%          bytes: 34605930
%
% 
%   Example 3: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-8); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappx_g(@(x) x.^2, in_param)
%
% out_param = 
% 
%           a: -5
%      abstol: 1.0000e-08
%           b: 5
%           f: @(x)x.^2
%     maxiter: 1000
%         nhi: 500
%         nlo: 10
%        nmax: 10000000
%       nstar: [1x256 double]
%       ninit: 703
%        exit: [2x1 logical]
%        iter: 9
%     npoints: 179713
%      errest: 9.4370e-***9
%           x: [1x179713 double]
%       bytes: 15123674
%
%
%   Example 4:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
% out_param = 
% 
%                a: -2
%           abstol: 1.0000e-***7
%                b: 2
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            nstar: [1x512 double]
%            ninit: 85
%             exit: [2x1 logical]
%             iter: 10
%          npoints: 43009
%           errest: 5.9830e-***8
%                x: [1x43009 double]
%            bytes: 3664502
%
%
%   Example 5:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
% out_param = 
% 
%                a: -3
%           abstol: 1.0000e-***6
%                b: 0
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 40
%              nlo: 20
%             nmax: 10000000
%            nstar: [1x128 double]
%            ninit: 69
%             exit: [2x1 logical]
%             iter: 8
%          npoints: 8705
%           errest: 8.8908e-***7
%                x: [1x8705 double] 
%            bytes: 747010
%
%
%   Example 6:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-2,5,1e-6,20,50)
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
%       nstar: [1x256 double]
%       ninit: 91
%        exit: [2x1 logical]
%        iter: 9
%     npoints: 23041
%      errest: 6.2507e-***7
%           x: [1x23041 double]
%       bytes: 1962794
% 
%
%   Example 7:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-3,3,1e-7,20,50)
%
% out_param = 
%                a: -3
%           abstol: 1.0000e-***7
%                b: 3
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            nstar: [1x1024 double]
%            ninit: 89
%             exit: [2x1 logical]
%             iter: 11
%          npoints: 90113
%           errest: 3.0223e-***8
%                x: [1x90113 double]
%            bytes: 7668266
%
%
%   Example 8:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-5,10,1e-7,10,20)
%
% out_param = 
% 
%                a: -5
%           abstol: 1.0000e-07
%                b: 10
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 20
%              nlo: 10
%             nmax: 10000000
%            nstar: [1x4096 double]
%            ninit: 41
%             exit: [2x1 logical]
%             iter: 13
%          npoints: 163841
%           errest: 3.0999e-***8
%                x: [1x163841 double]
%            bytes: 14147114
