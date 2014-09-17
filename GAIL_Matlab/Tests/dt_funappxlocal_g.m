%DT_FUNAPPXLOCAL_G small doctest for funappxlocal_g
%
%   >> funappxlocal_g
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; pp = funappxlocal_g(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x3969 double]
%      coefs: [3968x2 double]
%     pieces: 3968
%      order: 2
%        dim: 1
%     orient: 'first'
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [pp out_param] = funappxlocal_g(f)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x3969 double]
%      coefs: [3968x2 double]
%     pieces: 3968
%      order: 2
%        dim: 1
%     orient: 'first'
%
% out_param = 
% 
%              f: @(x)x.^2
%              a: 0
%              b: 1
%         abstol: 1.0000e-06
%            nlo: 20
%            nhi: 50
%           nmax: 10000000
%          ninit: 32
%        npoints: 3969
%     errorbound: 7.5421e-07
%          nstar: [1x128 double]
%           iter: 8
%
%
%   Example 3:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x770049 double]
%      coefs: [770048x2 double]
%     pieces: 770048
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -10
%         abstol: 1.0000e-08
%              b: 10
%              f: @(x)x.^2
%            nhi: 50
%            nlo: 20
%           nmax: 10000000
%          ninit: 48
%        npoints: 770049
%     errorbound: 5.2640e-09
%          nstar: [1x16384 double]
%           iter: 15
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 100; in_param.nhi = 500;
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x110337 double]
%      coefs: [110336x2 double]
%     pieces: 110336
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -5
%         abstol: 1.0000e-06
%              b: 5
%              f: @(x)x.^2
%            nhi: 500
%            nlo: 100
%           nmax: 10000000
%          ninit: 432
%        npoints: 110337
%     errorbound: 2.8441e-07
%          nstar: [1x256 double]
%           iter: 9
%   
%
%   Example 5:
%
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x41985 double]
%      coefs: [41984x2 double]
%     pieces: 41984
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-07
%              b: 2
%              f: @(x)x.^2
%            nhi: 50
%            nlo: 20
%           nmax: 10000000
%          ninit: 42
%        npoints: 41985
%     errorbound: 7.8394e-08
%          nstar: [1x1024 double]
%           iter: 11
%
%
%   Example 6:
%
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x16897 double]
%      coefs: [16896x2 double]
%     pieces: 16896
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -3
%         abstol: 1.0000e-06
%              b: 0
%              f: @(x)x.^2
%            nhi: 40
%            nlo: 20
%           nmax: 10000000
%          ninit: 34
%        npoints: 16897
%     errorbound: 3.4229e-07
%          nstar: [1x512 double]
%           iter: 10
%
%
%   Example 7:
%
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2,-2,5)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x22529 double]
%      coefs: [22528x2 double]
%     pieces: 22528
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-06
%              b: 5
%              f: @(x)x.^2
%            nhi: 50
%            nlo: 20
%           nmax: 10000000
%          ninit: 45
%        npoints: 22529
%     errorbound: 7.8881e-07
%          nstar: [1x512 double]
%           iter: 10
%
%
%   Example 8:
%
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2,-3,3,1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x88065 double]
%      coefs: [88064x2 double]
%     pieces: 88064
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% 
% out_param = 
% 
%              a: -3
%         abstol: 1.0000e-07
%              b: 3
%              f: @(x)x.^2
%            nhi: 50
%            nlo: 20
%           nmax: 10000000
%          ninit: 44
%        npoints: 88065
%     errorbound: 3.8587e-08
%          nstar: [1x2048 double]
%           iter: 12
%
%
%   Example 9:
%
%   >> [pp, out_param] = funappxlocal_g(@(x) x.^2,-5,10,1e-7,10,20)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x155649 double]
%      coefs: [155648x2 double]
%     pieces: 155648
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -5
%         abstol: 1.0000e-07
%              b: 10
%              f: @(x)x.^2
%            nhi: 20
%            nlo: 10
%           nmax: 10000000
%          ninit: 20
%        npoints: 155649
%     errorbound: 3.7614e-08
%          nstar: [1x8192 double]
%           iter: 14
%
