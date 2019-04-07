classdef optPrice < optPayoff

    %% OPTPRICE computes option price using (quasi-)Monte Carlo methods
    % 
    %
    % Example 1
    % >> obj = optPrice
    % obj =***
    %    optPrice with properties:
    % 
    %                   inputType: 'n'
    %          timeDim_timeVector: [1 2 3]
    %           timeDim_startTime: 1
    %             timeDim_endTime: 3
    %            timeDim_initTime: 0
    %           timeDim_initValue: 10
    %                 timeDim_dim: 1
    %          wnParam_sampleKind: 'IID'
    %         wnParam_distribName: {'Gaussian'}
    %            wnParam_xDistrib: 'Gaussian'
    %        bmParam_assembleType: 'diff'
    %                bmParam_whBM: 1
    %         assetParam_pathType: 'GBM'
    %        assetParam_initPrice: 10
    %         assetParam_interest: 0.0100
    %        assetParam_meanShift: 0
    %       assetParam_volatility: 0.5000
    %           assetParam_nAsset: 1
    %         payoffParam_optType: {'euro'}
    %     payoffParam_putCallType: {'call'}
    %          payoffParam_strike: 10
    %                  exactPrice: 3.4501
    %        priceParam_cubMethod: 'IID_MC'
    %           priceParam_absTol: 1
    %           priceParam_relTol: 0
    %            priceParam_alpha: 0.0100
    % 
    %   ***


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
      allowCubMethod = {'IID_MC','Sobol','SobolCV','lattice','IID_MC_new', 'IID_MC_newtwo', ...
         'IID_MC_abs','IID_MC_CLT'} 
   end
   

%% Methods
% The constructor for |assetPath| uses the |brownianMotion| constructor
% and then parses the other properties. The function |genStockPaths| generates
% the asset paths based on |whiteNoise| paths.

   methods
        
      % Creating an asset path process
      function obj = optPrice(varargin)         
         obj@optPayoff(varargin{:}) %parse basic input
         if nargin>0
            val=varargin{1};
            if isa(val,'optPrice')
               obj.priceParam = val.priceParam;
               if nargin == 1
                  return
               end
            end
            if isfield(obj.restInput,'priceParam')
               val = obj.restInput.priceParam;
               obj.priceParam = val;
               obj.restInput = rmfield(obj.restInput,'priceParam');
            end       
         end
      end
      
      function set.priceParam(obj,val)
         if isfield(val,'cubMethod') %cubature method
            assert(any(strcmp(val.cubMethod,obj.allowCubMethod)))
            obj.priceParam.cubMethod=val.cubMethod;
            if any(strcmp(obj.priceParam.cubMethod, ...
                  {'IID_MC','IID_MC_new', 'IID_MC_newtwo','IID_MC_abs'}))
               obj.wnParam.sampleKind = 'IID';
            end
            if any(strcmp(obj.priceParam.cubMethod,{'Sobol','lattice'}))
               obj.inputType = 'x';
               obj.wnParam.sampleKind = obj.priceParam.cubMethod;
               obj.wnParam.xDistrib = 'Uniform';
            end           
            if any(strcmp(obj.priceParam.cubMethod,{'SobolCV'}))
               obj.inputType = 'x';
               obj.wnParam.sampleKind = 'Sobol';
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
         elseif strcmp(obj.priceParam.cubMethod,'IID_MC_CLT')
            [price, outtemp] = meanMC_CLT('Y', @(n) genOptPayoffs(obj,n), ...
               'absTol', obj.priceParam.absTol, 'relTol', obj.priceParam.relTol, ...
               'alpha', obj.priceParam.alpha, 'trueMuCV', obj.exactPrice(2:end));
            out.nPaths=outtemp.nSample;
         elseif strcmp(obj.priceParam.cubMethod,'Sobol')
            if strcmp(obj.payoffParam.optType,'american')
                [price, outtemp] = cubSobol_american_g(@(x) genOptPayoffs(obj,x), ...
                    [zeros(1,obj.timeDim.nCols); ones(1,obj.timeDim.nCols)], ...
                    obj.priceParam.absTol, obj.priceParam.relTol);
                 out.nPaths=outtemp.n;
            else
                [price, outtemp] = cubSobol_g(@(x) genOptPayoffs(obj,x), ...
                    [zeros(1,obj.timeDim.nCols); ones(1,obj.timeDim.nCols)], ...
                    'uniform', obj.priceParam.absTol, obj.priceParam.relTol);
                 out.nPaths=outtemp.n;
            end
         elseif strcmp(obj.priceParam.cubMethod,'SobolCV') %control variates
            f.func = @(x) genOptPayoffs(obj,x);
            f.cv = obj.exactPrice(2:end);
            if strcmp(obj.payoffParam.optType,'american')
                [price, outtemp] = cubSobol_american_g(f, ...
                    [zeros(1,obj.timeDim.nCols); ones(1,obj.timeDim.nCols)], ...
                    obj.priceParam.absTol, obj.priceParam.relTol);
                 out.nPaths=outtemp.n;
            else
                [price, outtemp] = cubSobol_g(f, ...
                    [zeros(1,obj.timeDim.nCols); ones(1,obj.timeDim.nCols)], ...
                    'uniform', obj.priceParam.absTol, obj.priceParam.relTol);
                 out.nPaths=outtemp.n;
            end
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

