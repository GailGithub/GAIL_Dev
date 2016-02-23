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
%              f: @(x)x.^2
%              a: -10
%              b: 10
%         abstol: 1.0000e-***8
%            nlo: 10
%            nhi: 1000
%          ninit: 804
%           nmax: 10000000
%        maxiter: 1000
%           exit: [2x1 logical]
%           iter: 8
%        npoints: 205569
%         errest: 2.8480e-***9
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
%           f: @(x)x.^2
%           a: -5
%           b: 5
%      abstol: 1.0000e-08
%         nlo: 10
%         nhi: 500
%       ninit: 351
%        nmax: 10000000
%     maxiter: 1000
%        exit: [2x1 logical]
%        iter: 8
%     npoints: 89601
%      errest: 3.7478e-***9
%      
%
%   Example 4:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
%   out_param = 
%
%                f: @(x)x.^2
%                a: -2
%                b: 2
%           abstol: 1.0000e-***7
%              nlo: 20
%              nhi: 50
%            ninit: 42
%             nmax: 10000000
%          maxiter: 1000
%             exit: [2x1 logical]
%             iter: 8
%          npoints: 10497
%           errest: 4.3699e-***8
%
%
%   Example 5:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
%   out_param = 
% 
%                f: @(x)x.^2
%                a: -3
%                b: 0
%           abstol: 1.0000e-***6
%              nlo: 20
%              nhi: 40
%            ninit: 34
%             nmax: 10000000
%          maxiter: 1000
%             exit: [2x1 logical]
%             iter: 6
%          npoints: 2113
%           errest: 6.1248e-***7
%
% 
%   Example 6:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-2,5,1e-6,20,50)
%   
%   out_param = 
%
%           f: @(x)x.^2 
%           a: -2
%           b: 5
%      abstol: 1.0000e-***6
%         nlo: 20
%         nhi: 50
%       ninit: 45
%        nmax: 10000000
%     maxiter: 1000
%        exit: [2x1 logical]
%        iter: 7
%     npoints: 5633
%      errest: 4.6617e-***7
%
%
%   Example 7:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,-5,10,1e-7,10,20)
%
%   out_param = 
% 
%                f: @(x)x.^2
%                a: -5
%                b: 10
%           abstol: 1.0000e-07
%              nlo: 10
%              nhi: 20
%            ninit: 20
%             nmax: 10000000
%          maxiter: 1000
%             exit: [2x1 logical]
%             iter: 11
%          npoints: 38913
%           errest: 4.4596e-***8
%
%
%   Example 8:
%
%   >> [~, out_param] = funappxNoPenalty_g(@(x) x.^2,'memorytest',true,'output_x',1);
%   >> out_param.bytes <= 73444
%      1
%   >> length(out_param.x) == out_param.npoints
%      1
