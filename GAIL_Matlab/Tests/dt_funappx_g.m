%DT_FUNAPPX_G small doctest for funappx_g
%
%   >> funappx_g
%
%   Function f must be specified. Now GAIL is using f(x)=x^2 and unit interval [0,1].
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; pp = funappx_g(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x9901 double]
%      coefs: [9900x2 double]
%     pieces: 9900
%      order: 2
%        dim: 1
%     orient: 'first'
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [pp out_param] = funappx_g(f)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x9901 double]
%      coefs: [9900x2 double]
%     pieces: 9900
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
% 
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 10
%              nhi: 1000
%             nmax: 10000000
%            ninit: 100
%            nstar: 98
%     exceedbudget: 0
%          npoints: 9901
%       errorbound: 2.5508e-***9 
%
%
%   Example 3:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [pp, out_param] = funappx_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x5661151 double]
%      coefs: [5661150x2 double]
%     pieces: 5661150
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -10
%           abstol: 1.0000e-08
%                b: 10
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            ninit: 804
%            nstar: 802
%     exceedbudget: 0
%          npoints: 5661151
%       errorbound: 3.1299e-***12
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 100; in_param.nhi = 500;
%   >> [pp, out_param] = funappx_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x207743 double]
%      coefs: [207742x2 double]
%     pieces: 207742
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
% 
%                a: -5
%           abstol: 1.0000e-06
%                b: 5
%                f: @(x)x.^2
%              nhi: 500
%              nlo: 100
%             nmax: 10000000
%            ninit: 432
%            nstar: 430
%     exceedbudget: 0
%          npoints: 207743
%       errorbound: 5.7929e-***10
%   
%
%   Example 5:
%
%   >> clear in_param; in_param.a = -1; in_param.b = 1; in_param.Nmax = 10^6;
%   >> in_param.abstol = 10^(-6); in_param.nlo = 10; in_param.nhi = 500;  
%   >> [pp, out_param] = funappx_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x23221 double]
%      coefs: [23220x2 double]
%     pieces: 23220
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
% 
%                a: -1
%           abstol: 1.0000e-06
%                b: 1
%                f: @(x)x.^2
%              nhi: 500
%              nlo: 10
%             nmax: 1000000
%            ninit: 136
%            nstar: 134
%     exceedbudget: 0
%          npoints: 23221
%       errorbound: 1.8547e-***9
%
%
%   Example 6: 
%
%   >> [pp, out_param] = funappx_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x252333 double]
%      coefs: [252332x2 double]
%     pieces: 252332
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-07
%                b: 2
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            ninit: 399
%            nstar: 397
%     exceedbudget: 0
%          npoints: 252333
%       errorbound: 6.2823e-***11
%
%
%   Example 7:
%
%   >> [pp, out_param] = funappx_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x16765 double]
%      coefs: [16764x2 double]
%     pieces: 16764
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-06
%                b: 0
%                f: @(x)x.^2
%              nhi: 40
%              nlo: 20
%             nmax: 10000000
%            ninit: 34
%            nstar: 32
%     exceedbudget: 0
%          npoints: 16765
%       errorbound: 8.0062e-***9
%
%
%   Example 8:
%
%   >> [pp, out_param] = funappx_g(@(x) x.^2,'a',-30,'b',30,'nmax',1e7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x1826191 double]
%      coefs: [1826190x2 double]
%     pieces: 1826190
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -30
%           abstol: 1.0000e-06
%                b: 30
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            ninit: 928
%            nstar: 926
%     exceedbudget: 0
%          npoints: 1826191
%       errorbound: 2.6994e-***10
% 
%
%   Example 9:
%
%   >> [pp, out_param] = funappx_g(@(x) x.^2,-2,5)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x166353 double]
%      coefs: [166352x2 double]
%     pieces: 166352
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-06
%                b: 5
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            ninit: 563
%            nstar: 561
%     exceedbudget: 0
%          npoints: 166353
%       errorbound: 4.4267e-***10
%
%
%   Example 10:
%
%   >> f = @(x) x.^2; 
%   >> [pp, out_param] = funappx_g(f,-3,3,1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x431179 double]
%      coefs: [431178x2 double]
%     pieces: 431178
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-07
%                b: 3
%                f: @(x)x.^2
%              nhi: 1000
%              nlo: 10
%             nmax: 10000000
%            ninit: 518
%            nstar: 516
%     exceedbudget: 0
%          npoints: 431179
%       errorbound: 4.8410e-***11
%
%
%   Example 11:
%
%   >> f = @(x) x.^2; 
%   >> [pp, out_param] = funappx_g(f,-5,10,1e-7,10,20)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x195891 double]
%      coefs: [195890x2 double]
%     pieces: 195890
%      order: 2
%        dim: 1
%     orient: 'first'
%
%   out_param = 
%                a: -5
%           abstol: 1.0000e-07
%                b: 10
%                f: @(x)x.^2
%              nhi: 20
%              nlo: 10
%             nmax: 10000000
%            ninit: 20
%            nstar: 18
%     exceedbudget: 0
%          npoints: 195891
%       errorbound: 1.4659e-***9
%
