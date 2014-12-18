%DT_FUNAPPX_G small doctest for funappx_g
%
%   >> funappx_g
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 1:
%
%   >> f = @(x) exp(-100*(x-sqrt(2)/2).^2); [~, out_param] = funappx_g(f) 
%
% out_param = 
% 
%                f: @(x)exp(-100*(x-sqrt(2)/2).^2)
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%             nmax: 10000000
%          maxiter: 1000
%            ninit: 100
%             exit: 0
%             iter: 9
%          npoints: 6733
%           errest: 9.4644e-07
%            nstar: [1x68 double]
% 
%
%   Example 3:
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
%           exit: 0
%           iter: 10
%        npoints: 411137
%         errest: 5.9832e-09
%          nstar: [1x512 double]
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappx_g(@(x) x.^2, in_param)
%
% out_param = 
% 
%                a: -5
%           abstol: 1.0000e-06
%                b: 5
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 500
%              nlo: 10
%             nmax: 10000000
%            ninit: 351
%             exit: 0
%             iter: 7
%          npoints: 22401
%           errest: 7.7860e-07
%            nstar: [1x64 double]
%   
%
%   Example 5:
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
%             exit: 0
%             iter: 11
%          npoints: 41985
%           errest: 7.8394e-08
%            nstar: [1x1024 double]
%
%
%   Example 6:
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
%             exit: 0
%             iter: 10
%          npoints: 16897
%           errest: 3.4229e-07
%            nstar: [1x512 double]
%
%
%   Example 7:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-2,5,1e-6,20,50)
%   
% out_param = 
% 
%                a: -2
%           abstol: 1.0000e-06
%                b: 5
%                f: @(x)x.^2
%          maxiter: 1000
%              nhi: 50
%              nlo: 20
%             nmax: 10000000
%            ninit: 45
%             exit: 0
%             iter: 10
%          npoints: 22529
%           errest: 7.8881e-07
%            nstar: [1x512 double]
%
%
%   Example 8:
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
%             exit: 0
%             iter: 12
%          npoints: 88065
%           errest: 3.8587e-08
%            nstar: [1x2048 double]
%
%
%   Example 9:
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
%             exit: 0
%             iter: 14
%          npoints: 155649
%           errest: 3.7614e-08
%            nstar: [1x8192 double]
%
