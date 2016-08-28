function [fmin, out_param] = par_funmin_g(workers,varargin)
% par_funmin_g Executes GAIL's funmin_g using MATLAB Parallel Toolbox
%
%   [fmin, out_param] = par_funmin_g(w, f) minimizes a function f by applying the
%   GAIL algorithm funmin_g on w even subintervals of [0,1] using a
%   parallel pool of w MATLAB workers. If w equals 1, then no parallel
%   pool is created and it is equivalent to calling fmin = funmin_g(f).
%
%   Input Arguments
%      
%     w --- number of MATLAB parallel pool workers, positive integer
% 
%     Refer to documentation of funappx_g for additional optional inputs
%     
%
%   Output Arguments
%
%     Same as the outputs of funmin_g. Refer to documentation of funmin_g.
%
%
%   Example 1:
%   Create a pool of 8 workers:
%      myCluster = parcluster('local')
%      myCluster.NumWorkers = 8;
%      saveProfile(myCluster);
%      if isempty(gcp), parpool; end
%
%   >> f = @(x) x.^2;
%   >> [fmin, out_param] = par_funmin_g(4,f,-2,2,1e-7,100,100); out_param
%   out_param = 
%             f: @(x)x.^2
%             a: -2
%             b: 2
%        abstol: 1.0000e-07
%           nlo: 100
%           nhi: 100
%         ninit: 100
%          nmax: 10000000
%       maxiter: 1000
%      exitflag: [0 0]
%          iter: 12
%       npoints: 190
%        errest: ***e-08
%     intervals: ***
%      output_x: 0
% 
%  >> out_param.intervals
%     ans = 
%       1.0e-03 *
%        -0.3***
%         0.3***
%
%  >> abs(fmin - 0) < 1e-7
%    1
%
%
%  To release workers:
%   delete(gcp); 
%
%
%  See also FUNMIN_G, FMINBND, PAR_FUNAPPX_G 
%
%
%  References
%
%   [1]  Sou-Cheng T. Choi, Yuhan Ding, Fred J.Hickernell, Xin Tong, "Local
%   Adaption for Approximation and Minimization of Univariate Functions,"
%   working, 2016.
%
%   [2]  Xin Tong. A Guaranteed, "Adaptive, Automatic Algorithm for
%   Univariate Function Minimization," MS thesis, Illinois Institute of
%   Technology, 2014.
%
%   [2] Sou-Cheng T. Choi, Fred J. Hickernell, Yuhan Ding, Lan Jiang,
%   Lluis Antoni Jimenez Rugama, Xin Tong, Yizhi Zhang and Xuan Zhou,
%   GAIL: Guaranteed Automatic Integration Library (Version 2.1)
%   [MATLAB Software], 2015. Available from http://code.google.com/p/gail/
%
%   If you find GAIL helpful in your work, please support us by citing the
%   above papers, software, and materials.
%


in_param = gail.funmin_g_in_param(varargin{:});
out_param = in_param.toStruct();
f = in_param.f;

a = out_param.a;
b = out_param.b;
h = (b-a)/workers;
fa = zeros(1,workers);
ou = cell(1,workers);
abstol = out_param.abstol;
nlo = out_param.nlo;
nhi = out_param.nhi;
nmax = out_param.nmax;
parfor (i=1:workers, double(workers>1))
    aa=a+(i-1)*h;
    [fm,out] = funmin_g(f, aa, aa+h, abstol, ...
         max(5,ceil(nlo/workers)), max(5,ceil(nhi/workers)), ...
        ceil(nmax/workers));
    fa(i) = fm;   
    ou{i} = out;
end


% postprocessing
[fmin,index] = min(fa);
out_param = ou{index};
[~,indices] = find(fa == fmin);
indices_to_comine = find(diff(indices)==1);
if length(indices_to_comine) >= 1
  ind = indices(indices_to_comine(1));
  out_param = ou{ind};
  for i = indices_to_comine
    ind = indices(indices_to_comine(i));
    out_param.intervals = [min(out_param.intervals), max(ou{ind+1}.intervals)];
  end
end


out_param.a = a;
out_param.b = b;
out_param.abstol = abstol;
out_param.iter = ou{1}.iter;
for i = 2:workers,
  out_param.nlo = out_param.nlo + ou{i}.nlo;
  out_param.nhi = out_param.nhi + ou{i}.nhi;
  out_param.ninit = out_param.ninit + ou{i}.ninit;
  out_param.nmax = out_param.nmax + ou{i}.nmax;
  out_param.iter =  max(out_param.iter, ou{i}.iter);
  out_param.npoints = out_param.npoints + ou{i}.npoints - 1;
  out_param.errest =  max(out_param.errest, ou{i}.errest);
  out_param.exitflag = max(out_param.exitflag, ou{1}.exitflag);
end




