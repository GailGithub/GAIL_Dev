classdef multivarGauss < handle
   %MULTIVARGAUSS is a class that computes the probability of a
   %mulitivariate Gaussian distribution over a hyperbox
   
   properties
      a = -1.96 %left/lower endpoint vector of the hyperbox
      b = 1.96 %right/upper endpoint vector of the hyperbox
      mu = 0 %mean vector of the distribution
      Cov = 1 %covariance vector of the distribution
      absTol = 0.01 %absolute tolerance of the calculation
      relTol = 0 %relative error tolerance of the calculation
      n = 1e4 %number of samples
      intMeth = 'aff' %method for forming integrand
      cubMeth = 'Sobol' %cubature method for calculating the integral
      errMeth = 'n' %method for determining that the error is met
   end
   
   properties (SetAccess = private)
      CovProp %square root and determinant of covariance matrix
      f
   end
   
   methods
      function obj = multivarGauss(varargin)
         if nargin > 0
            iStart = 1;
            if isa(varargin{1},'multivarGauss')
               obj = copy(varargin{1});
               iStart = 2;
            end
           if nargin >= iStart
              wh = find(strcmp(varargin(iStart:end),'a'));
              if ~isempty(wh), obj.a = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'b'));
               if ~isempty(wh), obj.b = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'mu'));
               if ~isempty(wh), obj.mu = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'Cov'));
               if ~isempty(wh), obj.Cov = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'absTol'));
               if ~isempty(wh), obj.absTol = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'relTol'));
               if ~isempty(wh), obj.relTol = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'n'));
               if ~isempty(wh), obj.n = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'intMeth'));
               if ~isempty(wh), obj.intMeth = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'cubMeth'));
               if ~isempty(wh), obj.cubMeth = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'errMeth'));
               if ~isempty(wh), obj.errMeth = varargin{wh+iStart}; end
           end            
         end
         updateCovProp(obj)
      end
      
      function set.Cov(obj,val)
         obj.Cov = val;
         updateCovProp(obj);
         updateInteg(obj);
      end
      
       function set.a(obj,val)
         obj.a = val;
         updateInteg(obj);
       end
       
       function set.b(obj,val)
         obj.b = val;
         updateInteg(obj);
       end
       
       function set.mu(obj,val)
         obj.mu = val;
         updateInteg(obj);
       end
       
       function set.intMeth(obj,val)
         obj.intMeth = val;
         updateInteg(obj);
       end
             
       function updateCovProp(obj)
         obj.CovProp.C = chol(obj.Cov)';
         obj.CovProp.detSig = det(obj.Cov);
         obj.CovProp.invSig = inv(obj.Cov);
      end
      
      function updateInteg(obj)
         if strcmp(obj.intMeth,'aff')
            stretch = @(t) bsxfun(@plus,obj.a,bsxfun(@times,obj.b-obj.a,t));
            multinorm = ...
               @(x) exp(-0.5*sum((x - obj.mu).*((x - obj.mu)*obj.CovProp.invSig),2)) ...
               ./sqrt(((2*pi)^numel(obj.a))*obj.CovProp.detSig);
            obj.f = @(t) prod(obj.b-obj.a)*multinorm(stretch(t));
         elseif strcmp(obj.intMeth,'Genz')
            obj.f = @(t) Genz(t,obj);
         else
            error('intMeth not recognized')
         end
      end
      
      function fval = Genz(w,obj)
         dim = numel(obj.a);
         nn = size(w,1);
         am = obj.a - obj.mu;
         bm = obj.b - obj.mu;
         a1 = am(1)/obj.CovProp.C(1,1); 
         b1 = bm(1)/obj.CovProp.C(1,1); 
         d = gail.stdnormcdf(a1); 
         e = gail.stdnormcdf(b1);
         fval = (e-d)*ones(nn,1);
         y = zeros(nn,dim-1);
         for i = 2:dim
            y(:,i-1) = gail.stdnorminv(d+w(:,i-1).*(e-d));
            aux = sum(bsxfun(@times,obj.CovProp.C(i,1:i-1),y(:,1:i-1)),2);
            a1 = (am(i)-aux)/obj.CovProp.C(i,i);
            b1 = (bm(i)-aux)/obj.CovProp.C(i,i);
            d = gail.stdnormcdf(a1);
            e = gail.stdnormcdf(b1);
            fval = fval .* (e-d);
         end
      end
      
      function [prob,out] = compProb(obj)
         redDim = strcmp(obj.intMeth,'Genz');
         if strcmp(obj.errMeth,'n')
            out = [];
            dim = numel(obj.a);
            nmax = max(obj.n);
            if strcmp(obj.cubMeth,'IID')
               if dim - redDim >= 1
                  x = rand(nmax,dim-redDim);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'Sobol')
               if dim - redDim >= 1
                  x = net(scramble(sobolset(dim-redDim), ...
                  'MatousekAffineOwen'),nmax);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'uSobol')
               if dim - redDim >= 1
                  x = net(sobolset(dim-redDim),nmax);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'SobolOpt')
               if dim - redDim >= 1
                  x = net(scramble(sobolset(dim-redDim), ...
                     'MatousekAffineOwen'),nmax);
                  %acosx = acos(1-2*x)/pi;
                  nn = numel(obj.n);
                  prob(nn,1) = 0;
                  temp = obj.f(x);
%                  temp = obj.f(acosx);
                  for ii = 1:nn
                     nii = obj.n(ii);
                     [K,kvec] = kernelFun(x(1:nii,:),'Mat1');
%                     [K,kvec] = kernelFun(acosx(1:nii,:));
                     w = pinv(K)*kvec;
                     prob(ii) = w'*temp(1:nii);
                  end               
               else
                  prob = obj.f(0)*ones(size(obj.n));
               end
            end
         elseif strcmp(obj.errMeth,'g')
            dim = numel(obj.a)-1;
            if strcmp(obj.cubMeth,'IID')
               [prob, out] = meanMC_g(@(m) obj.f(rand(m,dim)), ...
                  obj.absTol,obj.relTol);
            elseif strcmp(obj.cubMeth,'Sobol')
               [prob, out] = cubSobol_g(obj.f,[zeros(1,dim); ones(1,dim)], ...
                  'uniform',obj.absTol,obj.relTol);
            elseif strcmp(obj.cubMeth,'lattice')
               [prob, out] = cubLattice_g(obj.f,[zeros(1,dim); ones(1,dim)], ...
                  'uniform',obj.absTol,obj.relTol);
            end
          end
      end
               
   end
   
end
