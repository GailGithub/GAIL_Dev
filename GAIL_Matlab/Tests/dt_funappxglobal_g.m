function dt_funappxglobal_g
%DT_FUNAPPXGLOBAL_G small doctest for funappxglobal_g
%
%   >> funappxglobal_g
%
%   Function f must be specified. Now GAIL is using f(x)=x^2 and unit interval [0,1].
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; pp = funappxglobal_g(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x5051 double]
%      coefs: [5050x2 double]
%     pieces: 5050
%      order: 2
%        dim: 1
%     orient: 'first'
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [pp out_param] = funappxglobal_g(f)
%   
%  pp = 
%       form: 'pp'
%     breaks: [1x5051 double]
%      coefs: [5050x2 double]
%     pieces: 5050
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-***6
%              nlo: 10
%              nhi: 1000
%             nmax: 10000000
%            nstar: 100
%            ninit: 102
%     exceedbudget: 0
%          npoints: 5051
%           errest: 9.9990e-***7
%
%
%   Example 3:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x5668811 double]
%      coefs: [5668810x2 double]
%     pieces: 5668810
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -10
%           abstol: 1.0000e-***8
%                b: 10
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            nstar: 804
%            ninit: 806
%     exceedbudget: 0
%          npoints: 5668811
%           errest: 2.5023e-***9
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 100; in_param.nhi = 500;
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x207841 double]
%      coefs: [207840x2 double]
%     pieces: 207840
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
% 
%                a: -5
%           abstol: 1.0000e-***6
%                b: 5
%                f: @(x)x.^2
%              nhi: 500
%              nlo: 100
%             nmax: 10000000
%            nstar: 432
%            ninit: 434
%     exceedbudget: 0
%          npoints: 207841
%           errest: 2.5053e-***7
%   
%
%   Example 5:
%
%   >> clear in_param; in_param.a = -1; in_param.b = 1; in_param.Nmax = 10^6;
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 500;  
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x23291 double]
%      coefs: [23290x2 double]
%     pieces: 23290
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
% 
%                a: -1
%           abstol: 1.0000e-***6
%                b: 1
%                f: @(x)x.^2
%              nhi: 500
%              nlo: 10
%             nmax: 1000000
%            nstar: 136
%            ninit: 138
%     exceedbudget: 0
%          npoints: 23291
%           errest: 2.5219e-***7
%
%
%   Example 6: 
%
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x252801 double]
%      coefs: [252800x2 double]
%     pieces: 252800
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-***7
%                b: 2
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            nstar: 399
%            ninit: 401
%     exceedbudget: 0
%          npoints: 252801
%           errest: 2.5013e-***8
%
%
%   Example 7:
%
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x17291 double]
%      coefs: [17290x2 double]
%     pieces: 17290
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-***6
%                b: 0
%                f: @(x)x.^2
%              nhi: 40
%              nlo: 20
%             nmax: 10000000
%            nstar: 34
%            ninit: 36
%     exceedbudget: 0
%          npoints: 17291
%           errest: 2.5639e-***7
%
%
%   Example 8:
%
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2,'a',-30,'b',30,'nmax',1e7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x1828273 double]
%      coefs: [1828272x2 double]
%     pieces: 1828272
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -30
%           abstol: 1.0000e-***6
%                b: 30
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            nstar: 928
%            ninit: 930
%     exceedbudget: 0
%          npoints: 1828273
%           errest: 2.4999e-***7
% 
%
%   Example 9:
%
%   >> [pp, out_param] = funappxglobal_g(@(x) x.^2,-2,5)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x83473 double]
%      coefs: [83472x2 double]
%     pieces: 83472
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-***6
%                b: 5
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            nstar: 563
%            ninit: 565
%     exceedbudget: 0
%          npoints: 83473
%           errest: 9.9654e-***7
%
%
%   Example 10:
%
%   >> f = @(x) x.^2; 
%   >> [pp, out_param] = funappxglobal_g(f,-3,3,1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x431809 double]
%      coefs: [431808x2 double]
%     pieces: 431808
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-***7
%                b: 3
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            nstar: 518
%            ninit: 520
%     exceedbudget: 0
%          npoints: 431809
%           errest: 2.5033e-***8
%
%
%   Example 11:
%
%   >> f = @(x) x.^2; 
%   >> [pp, out_param] = funappxglobal_g(f,-5,10,1e-7,10,20)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x207061 double]
%      coefs: [207060x2 double]
%     pieces: 207060
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -5
%           abstol: 1.0000e-***7
%                b: 10
%                f: @(x)x.^2
%              nhi: 20
%              nlo: 10
%             nmax: 10000000
%            nstar: 20
%            ninit: 22
%     exceedbudget: 0
%          npoints: 207061
%           errest: 2.6242e-***8
%
