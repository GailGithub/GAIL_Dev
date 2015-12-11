classdef assetPath < brownianMotion

%% assetPath
% is a class of discretized stochastic processes that model the values of
% an asset with respect to time. Browniam motions are used to build these
% asset paths.
% 
%
% Example 1
% >> obj = assetPath
% obj = 
%   assetPath with properties:
% 
%                 inputType: 'n'
%        timeDim_timeVector: [1 2 3]
%         timeDim_startTime: 1
%           timeDim_endTime: 3
%          timeDim_initTime: 0
%         timeDim_initValue: 10
%        wnParam_sampleKind: 'IID'
%       wnParam_distribName: 'Gaussian'
%          wnParam_xDistrib: 'Uniform'
%      bmParam_assembleType: 'diff'
%       assetParam_pathType: 'GBM'
%      assetParam_initPrice: 10
%       assetParam_interest: 0.0100
%     assetParam_volatility: 0.5000

% Authors: Fred J. Hickernell

%% Properties
% This process inherits properties from the |stochProcess| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class
   properties (SetAccess=public) %so they can only be set by the constructor
      assetParam = struct('pathType','GBM', ... %type of asset path
         'initPrice', 10, ... %initial asset price
         'interest', 0.01, ... %interest rate
         'volatility', 0.5,... %volatility      
         'drift', 0,... %drift
         'nAsset', 1,... %number of assets 
         'corrMat', 1,...%A transpose  
         'dividend',0,...% dividend
         'Vinst',0.04,...%instantaneous variance
         'Vlong',0.04,...%long term variance
         'kappa',0.5,...%mean reversion speed
         'epsilon',1,...%volatility of variance
         'rho',-0.9)%correlation of the Brownian motions
   end
   
   properties (Constant, Hidden) %do not change & not seen
      allowPathType = {'GBM','QE_m'} 
         %kinds of asset paths that we can generate
   end
   
   properties (Constant)
       psiC=1.5,% Between 1 and two. Threshold for different estimation distributions
       gamma1=0.5%For PredictorCorrector
       gamma2=0.5%For PredictorCorrector
   end
   
   properties (Dependent = true)
       sqCorr
       dT %time step
       k1
       k2
       k3
       K1
       K2
       K3
       K4
       c1 %adjustment due to drift
       A %futher adjustment
       
   end
 


%% Methods
% The constructor for |assetPath| uses the |brownianMotion| constructor
% and then parses the other properties. The function |genStockPaths| generates
% the asset paths based on |whiteNoise| paths.

   methods
        
      % Creating an asset path process
      function obj = assetPath(varargin)         
         obj@brownianMotion(varargin{:}) %parse basic input
         if nargin>0
            val=varargin{1};
            if isa(val,'assetPath')
               obj.assetParam = val.assetParam;
               if nargin == 1
                  return
               end
            end
            if isfield(obj.restInput,'assetParam')
               val = obj.restInput.assetParam;
               obj.assetParam = val;
               obj.restInput = rmfield(obj.restInput,'assetParam');
            end
         end
         obj.timeDim = struct('initTime',0, ...
            'initValue',obj.assetParam.initPrice,...
            'dim',obj.assetParam.nAsset);
      end
           
      % Set the properties of the payoff object
      function set.assetParam(obj,val)
         if isfield(val,'pathType') %data for type of option
            assert(any(strcmp(val.pathType,obj.allowPathType)))
            obj.assetParam.pathType=val.pathType; %row
         end
         if isfield(val,'nAsset') %data for number of assets
            validateattributes(val.nAsset,{'numeric'}, ...
               {'nonnegative'})
            obj.assetParam.nAsset=val.nAsset; %row
         end
         if isfield(val,'initPrice') %data for type of option
            validateattributes(val.initPrice,{'numeric'}, ...
               {'nonnegative'})
           if numel(val.initPrice) == obj.assetParam.nAsset
                obj.assetParam.initPrice=val.initPrice(:);
           else
              obj.assetParam.initPrice ...
                    =repmat(val.initPrice(1),obj.assetParam.nAsset,1);
           end  
         end
         if isfield(val,'interest') %data for type of option
            validateattributes(val.interest,{'numeric'}, ...
               {'nonnegative'})
            obj.assetParam.interest=val.interest(:); 
         end
         if isfield(val,'volatility') %data for type of option
            validateattributes(val.volatility,{'numeric'}, ...
               {'nonnegative'})
           if numel(val.volatility) == obj.assetParam.nAsset
            obj.assetParam.volatility=val.volatility(:); %row
           else
                obj.assetParam.volatility ...
                    =repmat(val.volatility(1),obj.assetParam.nAsset,1);
           end
         end
         if isfield(val,'corrMat') %data for A
            validateattributes(val.corrMat,{'numeric'}, ...
               {'nonnegative'})
            obj.assetParam.corrMat=val.corrMat; %row
         end
         if isfield(val,'drift') %data for type of option
            validateattributes(val.drift,{'numeric'}, ...
               {'scalar'})
            obj.assetParam.drift=val.drift; %row
         end
      end
      
      % Generate square root of correlation matrix
       function val = get.sqCorr(obj)
          [U,S] = svd(obj.assetParam.corrMat);
          val = sqrt(S)*U';
       end
      % Generate coefficients of QE mathod
      function val = get.dT(obj)
          val = (obj.timeDim.endTime-obj.timeDim.startTime)/(obj.timeDim.nSteps-1);
      end
      function val = get.k1(obj)
          val = exp(-obj.assetParam.kappa*obj.dT);
      end
      function val = get.k2(obj)
          val = obj.assetParam.epsilon^2*obj.k1.*(1-obj.k1)/obj.assetParam.kappa;
      end
      function val = get.k3(obj)
          val = exp(obj.assetParam.kappa*obj.dT)*0.5.*obj.k2.*(1-obj.k1)...
              .*obj.assetParam.Vlong;
      end
      function val = get.K1(obj)
          val = obj.gamma1*obj.dT*(obj.assetParam.kappa*obj.assetParam.rho...
              /obj.assetParam.epsilon - .5)-obj.assetParam.rho/obj.assetParam.epsilon;
      end 
      function val = get.K2(obj)
          val = obj.gamma2*obj.dT*(obj.assetParam.kappa*obj.assetParam.rho...
              /obj.assetParam.epsilon - .5)+obj.assetParam.rho/obj.assetParam.epsilon;
      end 
      function val = get.K3(obj)
          val = obj.gamma1*obj.dT*(1-obj.assetParam.rho^2);
      end 
      function val = get.K4(obj)
          val = obj.gamma2*obj.dT*(1-obj.assetParam.rho^2); 
      end
      function val = get.c1(obj)
          val = (obj.assetParam.interest-obj.assetParam.dividend)*obj.dT;
      end
      function val = get.A(obj)
          val = obj.K2+0.5*obj.K4;
      end
                  
      % Generate asset paths
      function [paths]=genPaths(obj,val)
         bmpaths = genPaths@brownianMotion(obj,val);
         nPaths = size(bmpaths,1);             
         if strcmp(obj.assetParam.pathType,'GBM')
            tempc=zeros(nPaths,obj.timeDim.nSteps);
            paths=zeros(nPaths,obj.timeDim.nCols);
            for idx=1:obj.assetParam.nAsset
              colRange = ...
                 ((idx-1)*obj.timeDim.nSteps+1):idx*obj.timeDim.nSteps;
              for j=1:obj.timeDim.nSteps
                 tempc(:,j)=bmpaths(:,j:obj.timeDim.nSteps:obj.timeDim.nCols) ...
                    * obj.sqCorr(:,idx);
              end
              paths(:,colRange) = obj.assetParam.initPrice(idx) * ...
                 exp(bsxfun(@plus,(obj.assetParam.interest ...
                 - obj.assetParam.volatility(idx).^2/2) ...
                 .* obj.timeDim.timeVector, obj.assetParam.volatility(idx)...
                 .* tempc));
             end
         end
         if strcmp(obj.assetParam.pathType,'QE_m')
             paths = zeros(nPaths,obj.timeDim.nSteps+1); %output pathS
             
             lnS1 = zeros(nPaths,obj.timeDim.nSteps+1); %logspot price path
             lnS1(:,1)= log(obj.assetParam.initPrice...
                 *exp(-obj.assetParam.dividend*(obj.timeDim.endTime...
                 -obj.timeDim.startTime)));% set S(0) adjust with dividend     
             V2 = zeros(nPaths,obj.timeDim.nSteps+1); % Variance path
             V2(:,1) = obj.assetParam.Vinst; % set V0
             
             UV1 = rand(nPaths,obj.timeDim.nSteps);%uniforms
             dW2 = bmpaths;% 
             K0 = zeros(nPaths,1);        % K0 for martingale adjust  
             
             for i=2:obj.timeDim.nSteps+1             % time loop
                m = obj.assetParam.Vlong + (V2(:,i-1)-obj.assetParam.Vlong)...
                    *obj.k1;% mean (moment matching)
                s2 = V2(:,i-1)*obj.k2 + obj.k3;   % var (moment matching)
                psi = s2./m.^2;   % psi compared to psiC

                psihat = 1./psi;
                b2 = 2*psihat - 1 + sqrt(2*psihat.*(2*psihat-1));
                a = m ./ (1 + b2);

                % Non-Central Chi squared approximation for psi < psiC
                I1 = find(psi<=obj.psiC); 
                I2 = ~I1;
                V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
        %         if isempty(I1)
        %         else
        %             V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
        %         end
                p = (psi - 1)./(psi + 1);               % for switching rule
                V2((UV1(:,i-1)<=p) & (psi>obj.psiC),i) = 0; % case u<=p & psi>psiC
                I1b = find((UV1(:,i-1)>p) & (psi>obj.psiC));% find is faster here!

                beta = (1 - p)./m;                      % for switching rule
                if isempty(I1b)
                else    % Psi^(-1)
                    V2(I1b,i) = log((1-p(I1b))./(1-UV1(I1b,i-1)))./beta(I1b);
                end
                % K0 for martingale adjustment
                K0(I1)= obj.c1-obj.A*b2(I1).*a(I1)./(1-2*obj.A*a(I1)) + 0.5*log(1-2*obj.A*a(I1));
                K0(I2)= obj.c1-log(p(I2)+beta(I2).*(1-p(I2))./(beta(I2)-obj.A));

                % log Euler Predictor-Corrector step
                lnS1(:,i) = lnS1(:,i-1) + K0 - (obj.K1+0.5...
                    *obj.K3).*V2(:,i-1) + ...
                    obj.K1.*V2(:,i-1) + obj.K2.*V2(:,i)...
                    + sqrt(obj.K3.*V2(:,i-1) + ...
                    obj.K4.*V2(:,i)).*dW2(:,i-1);
            end
            paths(:,:) = exp(lnS1);          
         end
      end
                 
   end

   methods (Access = protected)

      function propList = getPropertyList(obj)
         propList = getPropertyList@brownianMotion(obj);
         propList.assetParam_pathType = obj.assetParam.pathType;
         propList.assetParam_initPrice = obj.assetParam.initPrice;
         propList.assetParam_interest = obj.assetParam.interest;
         propList.assetParam_volatility = obj.assetParam.volatility;
         if obj.assetParam.drift ~=0
            propList.assetParam_drift = obj.assetParam.drift;
         end
         propList.assetParam_nAsset = obj.assetParam.nAsset;
      end

   end

end
