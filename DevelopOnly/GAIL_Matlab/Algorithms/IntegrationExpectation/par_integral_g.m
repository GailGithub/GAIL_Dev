function [q, out_param] = par_integral_g(workers,varargin)
%   Example 1:
%   Create a pool of 8 workers:
%      myCluster = parcluster('local')
%      myCluster.NumWorkers = 8;
%      saveProfile(myCluster);
%      if isempty(gcp), parpool; end
%
%   >> f = @(x) x.^2;
%   >> [q, out_param] = par_integral_g(4,f,-2,2,1e-7,20,20)
%   q = 
%      5.3333
%   out_param = 
% 
%    a: -2
%             abstol: 1.0000e-07
%                  b: 2
%                  f: @(x)x.^2
%            maxiter: 1000
%                nhi: 20
%                nlo: 20
%               nmax: 10000000
%              ninit: 20
%                tau: ***
%       exceedbudget: 0
%          tauchange: 0
%               iter: 3
%                  q: 5.3333
%            npoints: ***
%             errest: ***e-08
%
%
%  To release workers:
%   delete(gcp); 

in_param = gail.integral_g_in_param(varargin{:});
out_param = in_param.toStruct();
f = in_param.f;
q = 0;

a = out_param.a;
b = out_param.b;
h = (b-a)/workers;
%aa = zeros(1,workers); 
qa = zeros(1,workers);
ou = cell(1,workers);
%% parallel tasks
parfor i=1:workers,
    aa=a+(i-1)*h;
    [q1,out] = integral_g(f, aa, aa+h, out_param.abstol/workers, ...
        max(5,out_param.nlo/workers), max(5,out_param.nhi/workers), ...
        ceil(out_param.nmax/workers));
    qa(i) = q1;
    ou{i} = out;
end

%% postprocessing
q = sum(qa);
out_all = ou{1};
out_all.abstol = out_param.abstol;
out_all.q = q;
out_all.b = out_param.b;
out_all.nlo = ou{1}.nlo * workers;
out_all.nhi = out_all.nlo;
out_all.ninit = ou{1}.ninit * workers;
out_all.nmax = ou{1}.nmax * workers;
out_all.npoints = ou{1}.npoints;
for i = 2:workers,
    out_all.iter =  max(out_all.iter, ou{i}.iter);
    out_all.npoints = out_all.npoints + ou{i}.npoints;
    out_all.errest = out_all.errest +  ou{i}.errest;
    %out_all.exit = max(out_all.exit, ou{1}.exit);
end 
out_param = out_all;