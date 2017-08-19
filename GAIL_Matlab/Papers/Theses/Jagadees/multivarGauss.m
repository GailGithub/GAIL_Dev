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
      transMeth = 'none' %transformation method
      errMeth = 'n' %method for determining that the error is met
      bernPolyOrder = 2 %Bernoulli polynomial order for Gaussian cubature
      ptransform = 'C1sin'; %Periodization transform for Bayesian cubature
      fName = '' %function name for plot title
      figSavePath  = '' % path to save the figure
      arbMean = true  % Non zero mean m for Bayesian cubature
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
               wh = find(strcmp(varargin(iStart:end),'transMeth'));
               if ~isempty(wh), obj.transMeth = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'errMeth'));
               if ~isempty(wh), obj.errMeth = varargin{wh+iStart}; end

               wh = find(strcmp(varargin(iStart:end),'BernPolyOrder'));
               if ~isempty(wh), obj.bernPolyOrder = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'ptransform'));
               if ~isempty(wh), obj.ptransform = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'fName'));
               if ~isempty(wh), obj.fName = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'figSavePath'));
               if ~isempty(wh), obj.figSavePath = varargin{wh+iStart}; end
               wh = find(strcmp(varargin(iStart:end),'arbMean'));
               if ~isempty(wh), obj.arbMean = varargin{wh+iStart}; end
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

        function set.transMeth(obj,val)
         obj.transMeth = val;
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
            ff = @(t) prod(obj.b-obj.a)*multinorm(stretch(t));
         elseif strcmp(obj.intMeth,'Genz')
            ff = @(t) Genz(t,obj);
         else
            error('intMeth not recognized')
         end
         if strcmp(obj.transMeth,'none')
            obj.f = ff;
         elseif strcmp(obj.transMeth,'tent')
            obj.f = @(x) ff(1 - abs(2*x -1));
         elseif strcmp(obj.transMeth,'C0')
            obj.f = @(x) ff((x.^2).*(3-2*x)).*prod(6*x.*(1-x),2);
         else
            error ('transMeth not recognized')
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
         for i = 2:dim;
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
         dim = numel(obj.a);
         realDim = dim - redDim;
         out = [];
         if strcmp(obj.errMeth,'n')
            nmax = max(obj.n);
            if strcmp(obj.cubMeth,'IID')
               if realDim >= 1
                  x = rand(nmax,realDim);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'Sobol')
               if realDim >= 1
                  x = net(scramble(sobolset(realDim), ...
                  'MatousekAffineOwen'),nmax);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'uSobol')
               if realDim >= 1
                  x = net(sobolset(realDim),nmax);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            elseif strcmp(obj.cubMeth,'SobolOpt')
               if realDim >= 1
                  x = net(scramble(sobolset(realDim), ...
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
            elseif strcmp(obj.cubMeth,'SobolMLE')
               if realDim >= 1
                  [prob,out] = cubMLE(obj.f,obj.n,[zeros(1,realDim); ones(1,realDim)]);
               else
                  prob = obj.f(0)*ones(size(obj.n));
               end
            elseif strcmp(obj.cubMeth,'MLELattice')
               if realDim >= 1
                  %[prob,out] = cubMLE(obj.f,obj.n,[zeros(1,realDim); ones(1,realDim)],...
                  %    'Lattice1','Fourier','Thompson',obj.bernPolyOrder,obj.ptransform,...
                  %    obj.fName,obj.figSavePath);
                  testAll=true;
                  [prob, out] = cubMLELattice(obj.f, ...
                  realDim,obj.absTol,obj.relTol,obj.bernPolyOrder,obj.ptransform, ...
                  testAll,obj.figSavePath,obj.fName,obj.arbMean);
               else
                  prob = obj.f(0)*ones(size(obj.n));
               end
            elseif strcmp(obj.cubMeth,'lattice')
               if realDim >= 1
                  x = mod(bsxfun(@plus,gail.lattice_gen(1,nmax,realDim), ...
                     + rand(1,realDim)),1);
               else
                  x(nmax,1) = 0;
               end
               temp = cumsum(obj.f(x),1);
               prob = temp(obj.n)./obj.n(:);
            end
         elseif strcmp(obj.errMeth,'g')
            if strcmp(obj.cubMeth,'IID')
               [prob, out] = meanMC_g(@(m) obj.f(rand(m,realDim)), ...
                  obj.absTol,obj.relTol);
            elseif strcmp(obj.cubMeth,'Sobol')
               [prob, out] = cubSobol_g(obj.f, ...
                  [zeros(1,realDim); ones(1,realDim)], ...
                  'uniform',obj.absTol,obj.relTol);
            elseif strcmp(obj.cubMeth,'MLELattice')
               testAll=false;
               [prob, out] = cubMLELattice(obj.f, ...
                  realDim,obj.absTol,obj.relTol,obj.bernPolyOrder,obj.ptransform, ...
                  testAll,obj.figSavePath,obj.fName);
            end
          end
      end

      function val = sameProblem(obj1,obj2)
         val = all(obj1.a == obj2.a) && ...
            all(obj1.b == obj2.b) && ...
            all(obj1.mu == obj2.mu) && ...
            all(all(obj1.Cov == obj2.Cov)) && ...
            obj1.absTol == obj2.absTol && ...
            obj1.relTol == obj2.relTol && ...
            all(obj1.n == obj2.n) && ...
            strcmp(obj1.intMeth,obj2.intMeth) && ...
            strcmp(obj1.cubMeth,obj2.cubMeth) && ...
            strcmp(obj1.transMeth,obj2.transMeth) && ...
            strcmp(obj1.errMeth,obj2.errMeth);
   end


   end

end
