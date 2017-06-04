function dt_funappx_g
%DT_FUNAPPX_G fast doctest for funappx_g
%
%   Example 1:
%
%   >> funappx_g
%
%   Warning: Function f must be a function handle. 
%
%
%   Example 2:
%
%   >> clear in_param; in_param.a = -5; in_param.b =5; in_param.abstol = 10^(-7); 
%   >> [~, out_param] = funappx_g(@(x) x.^2, in_param)
% 
% out_param =***
% 
%            a: -5
%       abstol: 1.0000e-07
%            b: 5
%            f: @(x)x.^2
%      maxiter: 1000
%        ninit: 20
%         nmax: 10000000
%     exitflag: [0 0]
%         iter: 13
%      npoints: 81921
%       errest: 3.7262e-***8
% 
%
%   Example 3: 
%
%   >> clear in_param; f = @(x) sin(x); in_param.a = -1; in_param.b = 1; 
%   >> in_param.abstol = 10^(-8); in_param.ninit = 20; 
%   >> [~, out_param] = funappx_g(f, in_param)
% 
% out_param =***
% 
%            a: -1
%       abstol: 1.0000e-08
%            b: 1
%            f: @(x)sin(x)
%      maxiter: 1000
%        ninit: 20
%         nmax: 10000000
%     exitflag: [0 0]
%         iter: 12
%      npoints: 18159
%       errest: 2.5089e-***9
%      
%
%   Example 4:
%
%   >> [~, out_param] = funappx_g(@(x) x.^3,'a',-2,'b',2,'abstol',1e-7,'ninit',41)
%
% out_param =***
% 
%            a: -2
%       abstol: 1.0000e-07
%            b: 2
%            f: @(x)x.^3
%      maxiter: 1000
%        ninit: 41
%         nmax: 10000000
%     exitflag: [0 0]
%         iter: 12
%      npoints: 48787
%       errest: 3.4055e-***08
%
%
%   Example 5:
%
%   >> [~, out_param] = funappx_g(@(x) exp(-100*(x-0.7).^2),'a',0,'b',1,'ninit',28)
%
% out_param =***
% 
%            a: 0
%       abstol: 1.0000e-06
%            b: 1
%            f: @(x)exp(-100*(x-0.7).^2)
%      maxiter: 1000
%        ninit: 28
%         nmax: 10000000
%     exitflag: [0 0]
%         iter: 11
%      npoints: 6659
%       errest: 3.0439e-***7
%
%
%   Example 6:
%
%   >> [~, out_param] = funappx_g(@(x) x.^2,'memorytest',1,'output_x',1);
%   >> out_param.bytes <= 280674
%      1
%   >> length(out_param.x) == out_param.npoints
%      1
%
%
