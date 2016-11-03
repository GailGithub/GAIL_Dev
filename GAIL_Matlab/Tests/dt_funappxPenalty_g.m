function dt_funappxPenalty_g
%DT_FUNAPPX_G small doctest for funappxPenalty_g
%
%   Example 1:
%
%   >> funappxPenalty_g
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 2:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2, in_param)
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
%           exit: [0 0]
%           iter: 9
%        npoints: 411649
%         errest: 8.3292e-***9
%
% 
%   Example 3: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-8); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2, in_param)
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
%        exit: [0 0]
%        iter: 9
%     npoints: 179713
%      errest: 9.4370e-***9
%
%
%   Example 4:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
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
%            nstar: [1x512 double]
%            ninit: 85
%             exit: [0 0]
%             iter: 10
%          npoints: 43009
%           errest: 5.9830e-***8
%
%
%   Example 5:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
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
%             exit: [0 0]
%             iter: 8
%          npoints: 8705
%           errest: 8.8908e-***7
%
%
%   Example 6:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,-2,5,1e-6,20,50)
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
%        exit: [0 0]
%        iter: 9
%     npoints: 23041
%      errest: 6.2507e-***7
%
%
%   Example 7:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,-3,3,1e-7,20,50)
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
%             exit: [0 0]
%             iter: 11
%          npoints: 90113
%           errest: 3.0223e-***8
%
%
%   Example 8:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,-5,10,1e-7,10,20)
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
%             exit: [0 0]
%             iter: 13
%          npoints: 163841
%           errest: 3.0999e-***8
%
%
%   Example 9:
%
%   >> [~, out_param] = funappxPenalty_g(@(x) x.^2,'memorytest',true,'output_x',true);
%   >> out_param.bytes <= 390000
%   1
%   >> length(out_param.x) == out_param.npoints
%   1
