%DT_FUNAPPXLOCAL_G small doctest for funappxlocal_g_CSC
%
%   >> funappxlocal_g_CSC
%
%   Function f must be specified. Now GAIL is using f(x)=exp(-100*(x-0.5)^2) and unit interval [0,1].
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; pp = funappxlocal_g_CSC(f)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x1857 double]
%      coefs: [1856x2 double]
%     pieces: 1856
%      order: 2
%        dim: 1
%     orient: 'first'
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [pp out_param] = funappxlocal_g_CSC(f)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x1857 double]
%      coefs: [1856x2 double]
%     pieces: 1856
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
%            nlo: 9
%            nhi: 100
%          ninit: 30
%        npoints: 1857
%     errorbound: 7.7413e-07
%          nstar: [1x64 double]
%
%
%   Example 3:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2, in_param)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x364545 double]
%      coefs: [364544x2 double]
%     pieces: 364544
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
%            nhi: 100
%            nlo: 9
%          ninit: 90
%        npoints: 364545
%     errorbound: 6.5402e-09
%          nstar: [1x4096 double]
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 100; in_param.nhi = 500;
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2, in_param)
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
%          ninit: 432
%        npoints: 110337
%     errorbound: 2.6331e-07
%          nstar: [1x256 double]
%   
%
%   Example 5:
%
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x31233 double]
%      coefs: [31232x2 double]
%     pieces: 31232
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-07
%              b: 2
%              f: @(x)x.^2
%            nhi: 100
%            nlo: 9
%          ninit: 62
%        npoints: 31233
%     errorbound: 3.7139e-08
%          nstar: [1x512 double]
%
%
%   Example 6:
%
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
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
%          ninit: 34
%        npoints: 16897
%     errorbound: 3.4229e-07
%          nstar: [1x512 double]
%
%
%   Example 7:
%
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2,-2,5)
%   
% pp = 
% 
%       form: 'pp'
%     breaks: [1x18945 double]
%      coefs: [18944x2 double]
%     pieces: 18944
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -2
%         abstol: 1.0000e-06
%              b: 5
%              f: @(x)x.^2
%            nhi: 100
%            nlo: 9
%          ninit: 75
%        npoints: 18945
%     errorbound: 3.0204e-07
%          nstar: [1x256 double]
%
%
%   Example 8:
%
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2,-3,3,1e-7)
%
% pp = 
% 
%       form: 'pp'
%     breaks: [1x35841 double]
%      coefs: [35840x2 double]
%     pieces: 35840
%      order: 2
%        dim: 1
%     orient: 'first'
% 
% out_param = 
% 
%              a: -3
%         abstol: 1.0000e-07
%              b: 3
%              f: @(x)x.^2
%            nhi: 100
%            nlo: 9
%          ninit: 71
%        npoints: 35841
%     errorbound: 6.2381e-08
%          nstar: [1x512 double]
%
%
%   Example 9:
%
%   >> [pp, out_param] = funappxlocal_g_CSC(@(x) x.^2,-5,10,1e-7,10,20)
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
%          ninit: 20
%        npoints: 155649
%     errorbound: 3.7614e-08
%          nstar: [1x8192 double]
%
