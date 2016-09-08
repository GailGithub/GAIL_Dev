function dt_funappx_g
%DT_funappx_g small doctest for funappx_g
%
%   Example 1:
%
%   >> funappx_g
%
%   Warning: Function f must be a function handle. ***
%
%
%   Example 2:
%
%   >> clear in_param; in_param.a = -5; in_param.b =5; in_param.abstol = 10^(-7); 
%   >> [~, out_param] = funappx_g(@(x) x.^2, in_param)
%
% out_param = 
% 
%            f: @(x)x.^2
%            a: -5
%            b: 5
%       abstol: 1.0000e-07
%          nlo: 10
%          nhi: 1000
%        ninit: 658
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 24
%      npoints: 42049
%       errest: 4.2866e-***8
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
%            f: @(x)x.^2
%            a: -5
%            b: 5
%       abstol: 1.0000e-08
%          nlo: 10
%          nhi: 500
%        ninit: 351
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 32
%      npoints: 89601
%       errest: 9.3664e-***9
%      
%
%   Example 4:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
% out_param = 
% 
%            f: @(x)x.^2
%            a: -2
%            b: 2
%       abstol: 1.0000e-07
%          nlo: 20
%          nhi: 50
%        ninit: 42
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 19
%      npoints: 20993
%       errest: 2.7266e-***8
%
%
%   Example 5:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
%   out_param = 
% 
%            f: @(x)x.^2
%            a: -3
%            b: 0
%       abstol: 1.0000e-06
%          nlo: 20
%          nhi: 40
%        ninit: 34
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 22
%      npoints: 4225
%       errest: 3.8024e-***7
%
% 
%   Example 6:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-2,5,1e-6,20,50)
%   
%   out_param = 
%
%            f: @(x)x.^2
%            a: -2
%            b: 5
%       abstol: 1.0000e-06
%          nlo: 20
%          nhi: 50
%        ninit: 45
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 34
%      npoints: 11265
%       errest: 2.9039e-***7
%
%
%   Example 7:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,-5,10,1e-7,10,20)
%
% out_param = 
% 
%            f: @(x)x.^2
%            a: -5
%            b: 10
%       abstol: 1.0000e-07
%          nlo: 10
%          nhi: 20
%        ninit: 20
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 59
%      npoints: 77825
%       errest: 2.7867e-***8
%
%
%   Example 8:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'memorytest',1,'output_x',1);
%   >> out_param.bytes <= 168199
%      1
%   >> length(out_param.x) == out_param.npoints
%      1
