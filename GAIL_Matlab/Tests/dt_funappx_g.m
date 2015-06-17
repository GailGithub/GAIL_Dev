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
%         abstol: 1.0000e-08
%              b: 10
%              f: @(x)x.^2
%        maxiter: 1000
%            nhi: 1000
%            nlo: 10
%           nmax: 10000000
%          ninit: 804
%           exit: [2x1 logical]
%           iter: 10
%        npoints: 411137
%         errest: 3.5719e-09
%          nstar: [1x512 double]
%          bytes: 31303594
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
%       ninit: 351
%        exit: [2x1 logical]
%        iter: 10
%     npoints: 179201
%      errest: 3.9377e-09
%       nstar: [1x512 double]
%       bytes: 13676826
%
%
%   Example 4:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
% out_param = 
% 
%                a: -2
%           abstol: 1.0000e-07
%                b: 2
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            ninit: 42
%             exit: [2x1 logical]
%             iter: 11
%          npoints: 41985
%           errest: 2.9281e-08
%            nstar: [1x1024 double]
%            bytes: 3302326
%
%
%   Example 5:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
% out_param = 
% 
%                a: -3
%           abstol: 1.0000e-06
%                b: 0
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 40
%              nlo: 20
%             nmax: 10000000
%            ninit: 34
%             exit: [2x1 logical]
%             iter: 9
%          npoints: 8449
%           errest: 4.3863e-07
%            nstar: [1x256 double]
%            bytes: 672706
%
%
%   Example 6:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-2,5,1e-6,20,50)
%   
%   out_param = 
% 
%           a: -2
%      abstol: 1.0000e-06
%           b: 5
%           f: @(x)x.^2
%     maxiter: 1000
%         nhi: 50
%         nlo: 20
%        nmax: 10000000
%       ninit: 45
%        exit: [2x1 logical]
%        iter: 10
%     npoints: 22529
%      errest: 3.0527e-07
%       nstar: [1x512 double]
%       bytes: 1769322
% 
%
%   Example 7:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-3,3,1e-7,20,50)
%
% out_param = 
%                a: -3
%           abstol: 1.0000e-07
%                b: 3
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            ninit: 44
%             exit: [2x1 logical]
%             iter: 11
%          npoints: 44033
%           errest: 5.9080e-08
%            nstar: [1x1024 double]
%            bytes: 3457386
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
%            ninit: 20
%             exit: [2x1 logical]
%             iter: 13
%          npoints: 77825
%           errest: 5.9705e-08
%            nstar: [1x4096 double]
%            bytes: 6348138
