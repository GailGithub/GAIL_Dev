%DT_FUNAPPX_G small doctest for funappx_g
%
%   >> funappx_g
%
%   Function f must be specified. Now GAIL is using f(x)=x^2 and unit interval [0,1].
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; fappx = funappx_g(f)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [fappx out_param] = funappx_g(f)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%                f: @(x)x.^2
%                a: 0
%                b: 1
%           abstol: 1.0000e-06
%              nlo: 52
%              nhi: 52
%             nmax: 10000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 7039
%       errorbound: 5.0471e-***9
% 
%
%   Example 3:
%
%   >> clear in_param; in_param.a = -10; in_param.b =10; in_param.abstol = 10^(-8); 
%   >> [fappx, out_param] = funappx_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%                a: -10
%           abstol: 1.0000e-08
%                b: 10
%                f: @(x)x.^2
%              nhi: 52
%              nlo: 52
%             nmax: 10000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 1400359
%       errorbound: 5.1001e-***11
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.a = -5; in_param.b = 5; 
%   >> in_param.abstol = 10^(-6); in_param.nlo = 100; in_param.nhi = 500;
%   >> [fappx, out_param] = funappx_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
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
%   >> [fappx, out_param] = funappx_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
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
%   >> [fappx, out_param] = funappx_g(@(x) x.^2,'a',-2,'b',2,'abstol',1e-7)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-07
%                b: 2
%                f: @(x)x.^2
%              nhi: 52
%              nlo: 52
%             nmax: 10000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 88639
%       errorbound: 5.0912e-***10
%
%
%   Example 7:
%
%   >> [fappx, out_param] = funappx_g(@(x) x.^2,'a',-3,'b',0,'nlo',20,'nhi',40)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
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
%   >> [fappx, out_param] = funappx_g(@(x) x.^2,'a',-30,'b',30,'nmax',1e6)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%                a: -30
%           abstol: 1.0000e-06
%                b: 30
%                f: @(x)x.^2
%              nhi: 52
%              nlo: 52
%             nmax: 1000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 420139
%       errorbound: 5.0987e-***9
% 
%
%   Example 9:
%
%   >> [fappx, out_param] = funappx_g(@(x) x.^2,-2,5)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%                a: -2
%           abstol: 1.0000e-06
%                b: 5
%                f: @(x)x.^2
%              nhi: 52
%              nlo: 52
%             nmax: 10000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 49063
%       errorbound: 5.0892e-***9
%
%
%   Example 10:
%
%   >> f = @(x) x.^2; 
%   >> [fappx, out_param] = funappx_g(f,-3,3,1e-7)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%                a: -3
%           abstol: 1.0000e-07
%                b: 3
%                f: @(x)x.^2
%              nhi: 52
%              nlo: 52
%             nmax: 10000000
%            ninit: 52
%            nstar: 50
%     exceedbudget: 0
%          npoints: 132907
%       errorbound: 5.0951e-***10
%
%
%   Example 11:
%
%   >> f = @(x) x.^2; 
%   >> [fappx, out_param] = funappx_g(f,-5,10,1e-7,10,20)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
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
