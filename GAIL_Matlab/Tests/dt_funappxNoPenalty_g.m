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
%          nstar: 804
%          ninit: 1609
%           exit: [2x1 logical]
%           iter: 8
%        npoints: 205825
%         errest: ***.***e-***9
%              x: [1x205825 double]
%          bytes: 23673834
%
% 
%   Example 3: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-8); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2, in_param)
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
%       nstar: 351
%       ninit: 703
%        exit: [2x1 logical]
%        iter: 8
%     npoints: 89857
%      errest: ***.***e-***9
%           x: [1x89857 double]
%       bytes: 10337882
%
%
%   Example 4:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
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
%            nstar: 42
%            ninit: 85
%             exit: [2x1 logical]
%             iter: 8
%          npoints: 10753
%           errest: ***.***e-***8
%                x: [1x10753 double]
%            bytes: 1241078
%
%
%   Example 5:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
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
%            nstar: 34
%            ninit: 69
%             exit: [2x1 logical]
%             iter: 6
%          npoints: 2177
%           errest: ***.***e-***7
%                x: [1x2177 double] 
%            bytes: 254594
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
%       nstar: 45
%       ninit: 91
%        exit: [2x1 logical]
%        iter: 7
%     npoints: 5761
%      errest: ***.***e-***7
%           x: [1x5761 double]
%       bytes: 666410
% 
%
%   Example 7:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-3,3,1e-7,20,50)
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
%            nstar: 44
%            ninit: 89
%             exit: [2x1 logical]
%             iter: 8
%          npoints: 11265
%           errest: ***.***e-***8
%                x: [1x11265 double]
%            bytes: 1299370
%
%
%   Example 8:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-5,10,1e-7,10,20)
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
%            nstar: 20
%            ninit: 41
%             exit: [2x1 logical]
%             iter: 11
%          npoints: 40961
%           errest: ***.***e-***8
%                x: [1x40961 double]
%            bytes: 4714410
