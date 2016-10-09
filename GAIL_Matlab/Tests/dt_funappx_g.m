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
%         iter: 7
%      npoints: 42049
%       errest: ***e-***8
%
% 
%   Example 3: 
%
%   >> clear in_param; f = @(x) sin(x); in_param.a = -1; in_param.b = 1; 
%   >> in_param.abstol = 10^(-8); in_param.nlo = 10; in_param.nhi = 500;
%   >> [~, out_param] = funappx_g(f, in_param)
%
% out_param = 
% 
%            f: @(x)sin(x)
%            a: -1
%            b: 1
%       abstol: 1.0000e-08
%          nlo: 10
%          nhi: 500
%        ninit: 136
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 8
%      npoints: 14361
%       errest: ***e-***9
%      
%
%   Example 4:
%
%   >> [~, out_param] = funappx_g(@(x) x.^3,'a',-2,'b',2,'abstol',1e-7,'nlo',20,'nhi',50)
%
% out_param = 
% 
%            f: @(x)x.^3
%            a: -2
%            b: 2
%       abstol: 1.0000e-07
%          nlo: 20
%          nhi: 50
%        ninit: 42
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 11
%      npoints: 33211
%       errest: ***e-***8
%
%
%   Example 5:
%
%   >> [~, out_param] = funappx_g(@(x) exp(-100*(x-0.7).^2),'a',0,'b',1,'nlo',20,'nhi',40)
%
%   out_param = 
% 
%            f: @(x)exp(-100*(x-0.7).^2)
%            a: 0
%            b: 1
%       abstol: 1.0000e-06
%          nlo: 20
%          nhi: 40
%        ninit: 29
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0 0 0 0]
%         iter: 10
%      npoints: 4573
%       errest: ***e-***7
%
%
%   Example 6:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'memorytest',1,'output_x',1);
%   >> out_param.bytes <= 189591
%      1
%   >> length(out_param.x) == out_param.npoints
%      1
