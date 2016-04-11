function [fappx, out_param] = par_funappx_g(varargin)
%   Example 1:
%
%   >> f = @(x) x.^2;
%   >> [fappx, out_param] = par_funappx_g(f,-2,2,1e-7,20,20); out_param
% 
%   out_param = 
% 
%            f: @(x)x.^2
%            a: -2
%            b: 2
%       abstol: 1.0000e-07
%          nlo: 20
%          nhi: 20
%        ninit: 20
%         nmax: 10000000
%      maxiter: 1000
%     exitflag: [0 0]
%         iter: 10
%      npoints: ***
%       errest: ***e-08
% 
%  >> x=-2:1e-7:2; f = @(x) x.^2; err = max(abs(f(x)-fappx(x))) 
%     err =
%       ***e-08

in_param = gail.funappx_g_in_param(varargin{:});
out_param = in_param.toStruct();
f = in_param.f;
MATLABVERSION = gail.matlab_version;

workers = 4;
a = out_param.a;
b = out_param.b;
h = (b-a)/workers;
aa = zeros(1,workers); 
fa = cell(1,workers);
ou = cell(1,workers);
parfor i=1:workers,
    aa=a+(i-1)*h;
    [fappx,out] = funappxNoPenalty_g(f, aa, aa+h, out_param.abstol, out_param.nlo/workers, out_param.nhi/workers, out_param.nmax/workers, out_param.maxiter, 'output_x', 1);
    fa{i} = fappx;
    ou{i} = out;
end
for i = 1:workers,
    if i == 1,
        out_all.f = ou{1}.f;
        out_all.a = ou{1}.a;
        out_all.b = ou{workers}.b;
        out_all.abstol = ou{1}.abstol;
        out_all.nlo = ou{1}.nlo * workers;
        out_all.nhi = out_all.nlo;
        out_all.ninit = ou{1}.ninit * workers;
        out_all.nmax = ou{1}.nmax * workers;
        out_all.maxiter = ou{1}.maxiter;
        out_all.exitflag = ou{1}.exitflag;
        out_all.iter = ou{1}.iter;
        out_all.npoints = ou{1}.npoints;
        out_all.errest = ou{1}.errest;
        x = ou{1}.x;
        y = ou{1}.y;
    else
        out_all.iter =  max(out_all.iter, ou{i}.iter);
        out_all.npoints = out_all.npoints + ou{i}.npoints - 1;
        out_all.errest =  max(out_all.errest, ou{i}.errest);
        out_all.exitflag = max(out_all.exitflag, ou{1}.exitflag);
        x = [x ou{i}.x(2:end)];
        y = [y ou{i}.y(2:end)];
    end
end

if MATLABVERSION >= 8.3
    fappx = griddedInterpolant(x,y,'linear');
else
    fappx = @(t) ppval(interp1(x,y,'linear','pp'), t);
end;
out_param = out_all;


