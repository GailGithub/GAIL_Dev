%DT_FUNAPPX01_G small doctest for funappx01_g
%
%   >> funappx01_g
%
%   Function f must be specified. Now GAIL is using f(x)=x^2. 
%
%
%   Example 1:
%   
%   >> f = @(x) x.^2; fappx = funappx01_g(f)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%
%   Example 2:
%
%   >> f = @(x) x.^2; [fappx out_param] = funappx01_g(f)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%            ninit: 52
%             nmax: 10000000
%            nstar: 50
%     exceedbudget: 0
%          npoints: 7039
%         errbound: 5.0471e-***9
% 
%
%   Example 3:
%
%   >> clear in_param; in_param.abstol = 10^(-8); 
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 52
%             nmax: 10000000
%            nstar: 50    
%     exceedbudget: 0
%          npoints: 70075
%         errbound: 5.0913e-***11
%
% 
%   Example 4: 
%
%   >> clear in_param; in_param.ninit = 10; 
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 10000000
%            nstar: 8
%     exceedbudget: 0
%          npoints: 2683
%         errbound: 3.4755e-***8
%   
%
%   Example 5:
%
%   >> clear in_param; in_param.Nmax = 10^6; 
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%                f: @(x)x.^2
%            ninit: 52
%             nmax: 1000000
%            nstar: 50
%     exceedbudget: 0
%          npoints: 7039
%         errbound: 5.0471e-***9
%
%
%   Example 6: 
%
%   >> clear in_param; in_param.abstol = 10^(-8); 
%   >> in_param.ninit = 10; in_param.Nmax = 10^6; 
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2, in_param)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 1000000
%            nstar: 8
%     exceedbudget: 0
%          npoints: 26677
%         errbound: 3.5132e-***10
%
%
%   Example 7:
%
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2,'abstol',1e-8)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 52
%             nmax: 10000000
%            nstar: 50    
%     exceedbudget: 0
%          npoints: 70075
%         errbound: 5.0913e-***11
%
%
%   Example 8:
%
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2,'ninit',10)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 10000000
%            nstar: 8
%     exceedbudget: 0
%          npoints: 2683
%         errbound: 3.4755e-***8
%
%
%   Example 9:
%
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2,'nmax',1e6)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
% 
%           abstol: 1.0000e-***6
%                f: @(x)x.^2
%            ninit: 52
%             nmax: 1000000
%            nstar: 50
%     exceedbudget: 0
%          npoints: 7039
%         errbound: 5.0471e-***9
%
%
%   Example 10:
%
%   >> f = @(x) x.^2; 
%   >> [fappx, out_param] = funappx01_g(f,'ninit',10,'nmax',1e6,'abstol',1e-8)
%
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 1000000
%            nstar: 8
%     exceedbudget: 0
%          npoints: 26677
%         errbound: 3.5132e-***10
%
%
%   Example 11:
%
%   >> [fappx, out_param] = funappx01_g(@(x) x.^2,1e-8,10,1e6)
%   
%   fappx = 
%
%       @(x)interp1(x1,y1,x,'linear')
%
%   out_param = 
%           abstol: 1.0000e-***8
%                f: @(x)x.^2
%            ninit: 10
%             nmax: 1000000
%            nstar: 8
%     exceedbudget: 0
%          npoints: 26677
%         errbound: 3.5132e-***10
%
