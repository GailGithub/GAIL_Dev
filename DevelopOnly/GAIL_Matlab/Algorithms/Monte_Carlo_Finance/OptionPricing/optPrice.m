classdef optPrice < optPayoff

%% optPrice
% is a object of that computes the price of an object via (quasi-)Monte
% Carlo methods.
% 
%
% Example 1
% >> obj=whiteNoise
% obj = 
%   whiteNoise with properties:
% 
%       distribName: 'Uniform'
%        sampleName: 'IID'
%          xDistrib: 'Uniform'
%        qrandState: []
%              name: 'WhiteNoise'
%        timeVector: [1 2 3]
%         startTime: 1
%           endTime: 3
%            nSteps: 3
%     timeIncrement: [1 1]
%               dim: 1
%             nCols: 3
%         inputType: 'n'

%% Properties
% This process inherits properties from the |stochProcess| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class

   properties (SetAccess = public) %so they can only be set by the constructor
      priceParam = struct('cubMethod', 'IID_MC', ... %type of pricing scheme
         'absTol', 1, ... %absolute tolerance
         'relTol', 0, ... %relative tolerance
         'alpha', 0.01) %alpha = uncertainty
      
   end
   
   properties (Constant, Hidden) %do not change & not seen
      allowCubMethod = {'IID_MC','Sobol','lattice','IID_MC_new', 'IID_MC_newtwo', ...
         'IID_MC_abs'} 
   end
   

%% Methods
% The constructor for |assetPath| uses the |brownianMotion| constructor
% and then parses the other properties. The function |genStockPaths| generates
% the asset paths based on |whiteNoise| paths.

   methods
        
      % Creating an asset path process
      function obj = optPrice(varargin)         
         obj@optPayoff(varargin{:}) %parse basic input
         if isfield(obj.restInput,'priceParam')
            val = obj.restInput.priceParam;
            obj.priceParam = val;
            obj.restInput = rmfield(obj.restInput,'priceParam');
         end
      end
      
      function set.priceParam(obj,val)
         if isfield(val,'cubMethod') %cubature method
            assert(any(strcmp(val.cubMethod,obj.allowCubMethod)))
            obj.priceParam.cubMethod=val.cubMethod;
            if any(strcmp(obj.priceParam.cubMethod,{'Sobol','lattice'}))
               obj.inputType = 'x';
               obj.wnParam.sampleKind = obj.priceParam.cubMethod;
               obj.wnParam.xDistrib = 'Uniform';
            end
         end 
         if isfield(val,'absTol') %absolute error tolerance
            assert(val.absTol >= 0)
            obj.priceParam.absTol=val.absTol;
         end 
         if isfield(val,'relTol') %relative error tolerance
            assert(val.relTol >= 0 && val.relTol <= 1)
            obj.priceParam.relTol=val.relTol;
         end   
         if isfield(val,'exactOptType') %exact option type
            obj.priceParam.exactOptType = val.exactOptType;
         end 
        
      end
           
      % Generate option prices
      function [price, out] = genOptPrice(obj)
         if strcmp(obj.priceParam.cubMethod,'IID_MC')
            [price, outtemp] = meanMC_g(@(n) genOptPayoffs(obj,n), ...
               obj.priceParam.absTol, obj.priceParam.relTol, ...
               obj.priceParam.alpha);
            out.nPaths=outtemp.ntot;
         elseif strcmp(obj.priceParam.cubMethod,'IID_MC_abs')
            [price, outtemp] = meanMCabs_g(@(n) genOptPayoffs(obj,n), ...
               obj.priceParam.absTol, obj.priceParam.alpha);
            out.nPaths=outtemp.n;
         elseif strcmp(obj.priceParam.cubMethod,'IID_MC_new')
            [price, outtemp] = meanMCnew_g(@(n) genOptPayoffs(obj,n), ...
               obj.priceParam.absTol, obj.priceParam.relTol, ...
               obj.priceParam.alpha);
            out.nPaths=outtemp.ntot;
         elseif strcmp(obj.priceParam.cubMethod,'IID_MC_newtwo')
            [price, outtemp] = meanMCnew2_g(@(n) genOptPayoffs(obj,n), ...
               obj.priceParam.absTol, obj.priceParam.relTol, ...
               obj.priceParam.alpha);
            out.nPaths=outtemp.ntot;
         elseif strcmp(obj.priceParam.cubMethod,'Sobol')
            [price, outtemp] = cubSobol_g(@(x) genOptPayoffs(obj,x), ...
               [zeros(1,obj.timeDim.nCols); ones(1,obj.timeDim.nCols)], ...
               'uniform', obj.priceParam.absTol, obj.priceParam.relTol);
            out.nPaths=outtemp.n;
         elseif strcmp(obj.priceParam.cubMethod,'lattice')
            [price, outtemp] = cubLattice_g(@(x) genOptPayoffs(obj,x), ...
               [zeros(1,obj.timeDim.nSteps); ones(1,obj.timeDim.nSteps)], ...
               'uniform', obj.priceParam.absTol, obj.priceParam.relTol);
            out.nPaths=outtemp.n;
         end
         out.time=outtemp.time;
      end
       
   end
    
   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@optPayoff(obj);
         propList.priceParam_cubMethod = obj.priceParam.cubMethod;
         propList.priceParam_absTol = obj.priceParam.absTol;
         propList.priceParam_relTol = obj.priceParam.relTol;
         propList.priceParam_alpha = obj.priceParam.alpha;
      end

   end
end

