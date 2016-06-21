classdef assetPath < brownianMotion

%% assetPath
% is a class of discretized stochastic processes that model the values of
% an asset with respect to time. Browniam motions are used to construct
% these asset paths.
% 
%5/29/2016

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
% This process inherits properties from the |brownianMotion| class.  Below
% are values assigned to that are abstractly defined in that class plus
% some properties particulary for this class.

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
         'nu',1,...%volatility of variance
         'rho',-0.9)%correlation of the Brownian motions
   end
   
   properties (Constant, Hidden) %do not change & not seen

      allowPathType = {'GBM','QE_m','QE'} 
         %kinds of asset paths that we can generate
   end
   
   properties (Constant)
       psiC=1.5,% Between 1 and two. Threshold for different estimation distributions
%        gamma1=0.5%For PredictorCorrector
%        gamma2=0.5%For PredictorCorrector
   end
   
   properties (Dependent = true)
       sqCorr
%      dT %time step
%        k1
%        k2
%        k3
%        K1
%        K2
%        K3
%        K4
%        c1 %adjustment due to drift
%        A %futher adjustment
       %%%%%%%%%%%%%%%%%%%%%
%        c2
%        K0
       %%%%%%%%%%%%%%%%%%%%%
   end
 


%% Methods
% The constructor for |assetPath| uses the |brownianMotion| constructor and
% then parses the other properties. The function |genStockPaths| generates
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
            'dim',obj.assetParam.nAsset* ...
            (1+strcmp(obj.assetParam.pathType,'QE'))); %add extra BM for variance process
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
         if isfield(val,'dividend')
             validateattributes(val.dividend,{'numeric'},...
                 {'nonnegative'})
             obj.assetParam.dividend=val.dividend;
         end
         if isfield(val,'Vinst')
             validateattributes(val.Vinst,{'numeric'},...
                 {'scalar'})
             obj.assetParam.Vinst=val.Vinst;
         end
         if isfield(val,'Vlong')
             validateattributes(val.Vlong,{'numeric'},...
                 {'scalar'})
             obj.assetParam.Vlong=val.Vlong;
         end
         if isfield(val,'kappa')
             validateattributes(val.kappa,{'numeric'},...
                 {'scalar'})
             obj.assetParam.kappa=val.kappa;
         end
         if isfield(val,'nu')
             validateattributes(val.nu,{'numeric'},...
                 {'scalar'})
             obj.assetParam.nu=val.nu;
         end
         if isfield(val,'rho')
             validateattributes(val.rho,{'numeric'},...
                 {'scalar'})
             obj.assetParam.rho=val.rho;
         end
                  
      end
      
      % Generate square root of correlation matrix
       function val = get.sqCorr(obj)
          [U,S] = svd(obj.assetParam.corrMat);
          val = sqrt(S)*U';
       end
       
      % Generate coefficients of QE mathod
%       function val = get.dT(obj)
%           val = obj.timeDim.timeIncrement(1);
%                 %(obj.timeDim.endTime-obj.timeDim.startTime)/(obj.timeDim.nSteps-1);
%       end
%       function val = get.k1(obj)
%           val = exp(-obj.assetParam.kappa*obj.dT);
%       end
%       function val = get.k2(obj)
%           val = obj.assetParam.nu^2*obj.k1.*(1-obj.k1)/obj.assetParam.kappa;
%       end
%       function val = get.k3(obj)
%           val = exp(obj.assetParam.kappa*obj.dT)*0.5.*obj.k2.*(1-obj.k1)...
%               .*obj.assetParam.Vlong;
%       end
%       function val = get.K1(obj)
%           val = obj.gamma1*obj.dT*(obj.assetParam.kappa*obj.assetParam.rho...
%               /obj.assetParam.nu - .5)-obj.assetParam.rho/obj.assetParam.nu;
%       end 
%       function val = get.K2(obj)
%           val = obj.gamma2*obj.dT*(obj.assetParam.kappa*obj.assetParam.rho...
%               /obj.assetParam.nu - .5)+obj.assetParam.rho/obj.assetParam.nu;
%       end 
%       function val = get.K3(obj)
%           val = obj.gamma1*obj.dT*(1-obj.assetParam.rho^2);
%       end 
%       function val = get.K4(obj)
%           val = obj.gamma2*obj.dT*(1-obj.assetParam.rho^2); 
%       end
%       function val = get.c1(obj)
%           val = (obj.assetParam.interest-obj.assetParam.dividend)*obj.dT;
%       end
%       function val = get.A(obj)
%           val = obj.K2+0.5*obj.K4;
%       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       function val = get.c2(obj)
%           val = -obj.assetParam.rho*obj.assetParam.kappa*obj.assetParam.Vlong*obj.dT/...
%               obj.assetParam.nu;%used to determine K0
%       end
%       function val = get.K0(obj)
%           val = obj.c1+obj.c2;
%       end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
      % Generate asset paths
      function [paths]=genPaths(obj,val)
         bmpaths = genPaths@brownianMotion(obj,val);
         nPaths = size(bmpaths,1);             
         if strcmp(obj.assetParam.pathType,'GBM')
            tempc=zeros(nPaths,obj.timeDim.nSteps);
            paths=zeros(nPaths,obj.timeDim.nCols);
            % size_GBM = size(paths)
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
             
          % vector = obj.timeDim.timeVector        
%          dT = obj.timeDim.timeIncrement(1);
%           dT = obj.timeDim.timeVector(2)-obj.timeDim.timeVector(1);
          dT = (obj.timeDim.endTime-obj.timeDim.startTime)/(obj.timeDim.nSteps-1);
          gamma1 = (1-exp(obj.assetParam.kappa*dT)+obj.assetParam.kappa*dT)...
              /(obj.assetParam.kappa*dT*(1-exp(obj.assetParam.kappa*dT)));
          gamma2 = -(1-exp(obj.assetParam.kappa*dT)+obj.assetParam.kappa*dT*exp(obj.assetParam.kappa*dT))...
              /(obj.assetParam.kappa*dT*(1-exp(obj.assetParam.kappa*dT)));
          k1 = exp(-obj.assetParam.kappa*dT);
     
          k2 = obj.assetParam.nu^2*k1.*(1-k1)/obj.assetParam.kappa;
    
          k3 = exp(obj.assetParam.kappa*dT)*0.5.*k2.*(1-k1)...
              .*obj.assetParam.Vlong;
          K1 = gamma1*dT*(obj.assetParam.kappa*obj.assetParam.rho...
              /obj.assetParam.nu - .5)-obj.assetParam.rho/obj.assetParam.nu;
          K2 = gamma2*dT*(obj.assetParam.kappa*obj.assetParam.rho...
              /obj.assetParam.nu - .5)+obj.assetParam.rho/obj.assetParam.nu;
          K3 = gamma1*dT*(1-obj.assetParam.rho^2);
          K4 = gamma2*dT*(1-obj.assetParam.rho^2); 
          c1 = (obj.assetParam.interest-obj.assetParam.dividend)*dT;
                    
  %change timeDim.nSteps to Ntime***********************************%    
%              bmpaths = genPaths@brownianMotion(obj,val);
%              % size_bmpaths = size(bmpaths)
%              nPaths = size(bmpaths,1);  
%              paths = zeros(nPaths,obj.timeDim.nSteps); %output pathS
%              % size_paths = size (paths)
%              lnS1 = zeros(nPaths,obj.timeDim.nSteps); %logspot price path
%              lnS1(:,1)= log(obj.assetParam.initPrice...
%                  *exp(-obj.assetParam.dividend*obj.timeDim.endTime));
%              % set S(0) adjust with dividend 
%              V2 = zeros(nPaths,obj.timeDim.nSteps); % Variance path
%              V2(:,1) = obj.assetParam.Vinst; % set V0
%              UV1 = rand(nPaths,obj.timeDim.nSteps-1);%uniforms
%              temp = repmat(bmpaths,1);
%              dW2 = (temp(:,2:obj.timeDim.nSteps)- bmpaths(:,1:obj.timeDim.nSteps-1))/sqrt(dT);
% %              dW2 = randn(nPaths,obj.timeDim.nSteps-1);
%              K0 = zeros(nPaths,1);        % K0 for martingale adjust
%              A = K2+0.5*K4; %further adjustment
%******************************************************************%
%   %change Ntime to timeDim.nSteps ***********************************%    
             Ntime = floor(obj.timeDim.endTime/dT);
%              Ntime = obj.timeDim.nSteps-1;
             t0 = floor(obj.timeDim.startTime/dT);
             paths = zeros(nPaths,Ntime+1-t0);
             lnS1 = zeros(nPaths,Ntime+1);
             lnS1(:,1)= log(obj.assetParam.initPrice...
                 *exp(-obj.assetParam.dividend*obj.timeDim.endTime));
             % set S(0) adjust with dividend 
             
             V2 = zeros(nPaths,Ntime+1);
             V2(:,1) = obj.assetParam.Vinst; % set V0
            
             UV1 = rand(nPaths,Ntime);
             
             
             dW2 = randn(nPaths,Ntime);
%              dW2 = randn(nPaths,obj.timeDim.nSteps-1);
             K0 = zeros(nPaths,1);        % K0 for martingale adjust
             A = K2+0.5*K4; %further adjustment
%******************************************************************%
             for i=2:Ntime+1;%obj.timeDim.nSteps             % time loop
                m = obj.assetParam.Vlong + (V2(:,i-1)-obj.assetParam.Vlong)*k1;% mean (moment matching)
                s2 = V2(:,i-1)*k2 + k3;   % var (moment matching)
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
                K0(I1)= c1-A*b2(I1).*a(I1)./(1-2*A*a(I1)) + 0.5*log(1-2*A*a(I1));
                K0(I2)= c1-log(p(I2)+beta(I2).*(1-p(I2))./(beta(I2)-A));

                % log Euler Predictor-Corrector step
%                 lnS1(:,i) = lnS1(:,i-1) + K0 - (obj.K1+0.5...
%                     *obj.K3).*V2(:,i-1) + ...
%                     obj.K1.*V2(:,i-1) + obj.K2.*V2(:,i)...
%                     + sqrt(obj.K3.*V2(:,i-1) + ...
%                     obj.K4.*V2(:,i)).*(bmpaths(:,i)-bmpaths(:,i-1))/sqrt(obj.dT);%bmpaths(:,i)/sqrt(obj.timeDim.timeVector(i));

                lnS1(:,i) = lnS1(:,i-1) + K0 - (K1+0.5*K3).*V2(:,i-1) + ...
                    K1.*V2(:,i-1) + K2.*V2(:,i) + sqrt(K3.*V2(:,i-1) + ...
                    K4.*V2(:,i)).*dW2(:,i-1);
             end
            
             paths(:,:) = exp(lnS1(:,t0+1:end));   
%               paths = exp(lnS1);
%             paths(1,:)
         end
      
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           bmpaths = genPaths@brownianMotion(obj,val);
%              % size_bmpaths = size(bmpaths)
%              nPaths = size(bmpaths,1);      
%              
%              paths = zeros(nPaths,obj.timeDim.nSteps); %output pathS
%              % size_paths = size (paths)
%              lnS1 = zeros(nPaths,obj.timeDim.nSteps); %logspot price path
%              lnS1(:,1)= log(obj.assetParam.initPrice...
%                  *exp(-obj.assetParam.dividend*obj.timeDim.endTime));
%              % set S(0) adjust with dividend     
%              V2 = zeros(nPaths,obj.timeDim.nSteps); % Variance path
%              V2(:,1) = obj.assetParam.Vinst; % set V0
%              
%              UV1 = rand(nPaths,obj.timeDim.nSteps-1);%uniforms  
%              temp = repmat(bmpaths,1);
%              dW2 = (temp(:,2:obj.timeDim.nSteps)- bmpaths(:,1:obj.timeDim.nSteps-1))/sqrt(obj.dT);
%              %dW2 = randn(nPaths,obj.timeDim.nSteps-1);
%              K0 = zeros(nPaths,1);        % K0 for martingale adjust
%              A = obj.K2+0.5*obj.K4; %further adjustment
%              
%              for i=2:obj.timeDim.nSteps             % time loop
%                 m = obj.assetParam.Vlong + (V2(:,i-1)-obj.assetParam.Vlong)...
%                     *obj.k1;% mean (moment matching)
%                 s2 = V2(:,i-1)*obj.k2 + obj.k3;   % var (moment matching)
%                 psi = s2./m.^2;   % psi compared to psiC
% 
%                 psihat = 1./psi;
%                 b2 = 2*psihat - 1 + sqrt(2*psihat.*(2*psihat-1));
%                 a = m ./ (1 + b2);
% 
%                 % Non-Central Chi squared approximation for psi < psiC
%                 I1 = find(psi<=obj.psiC); 
%                 I2 = ~I1;
%                 V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
%         %         if isempty(I1)
%         %         else
%         %             V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
%         %         end
%                 p = (psi - 1)./(psi + 1);               % for switching rule
%                 V2((UV1(:,i-1)<=p) & (psi>obj.psiC),i) = 0; % case u<=p & psi>psiC
%                 I1b = find((UV1(:,i-1)>p) & (psi>obj.psiC));% find is faster here!
% 
%                 beta = (1 - p)./m;                      % for switching rule
%                 if isempty(I1b)
%                 else    % Psi^(-1)
%                     V2(I1b,i) = log((1-p(I1b))./(1-UV1(I1b,i-1)))./beta(I1b);
%                 end
%                 % K0 for martingale adjustment
%                 K0(I1)= obj.c1-A*b2(I1).*a(I1)./(1-2*A*a(I1)) + 0.5*log(1-2*A*a(I1));
%                 K0(I2)= obj.c1-log(p(I2)+beta(I2).*(1-p(I2))./(beta(I2)-A));
% 
%                 % log Euler Predictor-Corrector step
% %                 lnS1(:,i) = lnS1(:,i-1) + K0 - (obj.K1+0.5...
% %                     *obj.K3).*V2(:,i-1) + ...
% %                     obj.K1.*V2(:,i-1) + obj.K2.*V2(:,i)...
% %                     + sqrt(obj.K3.*V2(:,i-1) + ...
% %                     obj.K4.*V2(:,i)).*(bmpaths(:,i)-bmpaths(:,i-1))/sqrt(obj.dT);%bmpaths(:,i)/sqrt(obj.timeDim.timeVector(i));
% 
%                 lnS1(:,i) = lnS1(:,i-1) + K0 - (obj.K1+0.5...
%                     *obj.K3).*V2(:,i-1) + ...
%                     obj.K1.*V2(:,i-1) + obj.K2.*V2(:,i)...
%                     + sqrt(obj.K3.*V2(:,i-1) + ...
%                     obj.K4.*V2(:,i)).*dW2(:,i-1);
%             end
%             paths(:,:) = exp(lnS1);      
% %             paths(1,:)
%          end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if strcmp(obj.assetParam.pathType,'QE')
             
          % vector = obj.timeDim.timeVector        
%          dT = obj.timeDim.timeIncrement(1);
%           dT = obj.timeDim.timeVector(2)-obj.timeDim.timeVector(1);
          dT = obj.timeDim.timeIncrement(1);
          gamma1 = (1-exp(obj.assetParam.kappa*dT)+obj.assetParam.kappa*dT)...
              /(obj.assetParam.kappa*dT*(1-exp(obj.assetParam.kappa*dT)));
          gamma2 = -(1-exp(obj.assetParam.kappa*dT)+obj.assetParam.kappa*dT*exp(obj.assetParam.kappa*dT))...
              /(obj.assetParam.kappa*dT*(1-exp(obj.assetParam.kappa*dT)));
%           K1 = gamma1*dT*(obj.assetParam.kappa*obj.assetParam.rho...
%               /obj.assetParam.nu - .5)-obj.assetParam.rho/obj.assetParam.nu;
%           K2 = gamma2*dT*(obj.assetParam.kappa*obj.assetParam.rho...
%               /obj.assetParam.nu - .5)+obj.assetParam.rho/obj.assetParam.nu;
%           K3 = gamma1*dT*(1-obj.assetParam.rho^2);
%           K4 = gamma2*dT*(1-obj.assetParam.rho^2); 
%           c1 = (obj.assetParam.interest-obj.assetParam.dividend)*dT;
%           c2 = -obj.assetParam.rho*obj.assetParam.kappa*obj.assetParam.Vlong*dT/obj.assetParam.nu;
%           K0= c1+c2;
%   %change Ntime to timeDim.nSteps ***********************************%    
             Ntime = obj.timeDim.nSteps;
%              Ntime = obj.timeDim.nSteps-1;
             %t0 = floor(obj.timeDim.startTime/dT);
             paths = zeros(nPaths,Ntime);
             lnS1 = zeros(nPaths,Ntime+1);
             lnS1(:,1)= log(obj.assetParam.initPrice...
                 *exp(-obj.assetParam.dividend*obj.timeDim.endTime));
             % set S(0) adjust with dividend 
             
             V2 = zeros(nPaths,Ntime+1);
             V2(:,1) = obj.assetParam.Vinst; % set V0
            
             UV1 = rand(nPaths,Ntime);
             dW2 = randn(nPaths,Ntime);
%              dW2 = randn(nPaths,obj.timeDim.nSteps-1);
%******************************************************************%
%              for i=2:Ntime+1;%obj.timeDim.nSteps             % time loop
%                 m = obj.assetParam.Vlong + (V2(:,i-1)-obj.assetParam.Vlong)*k1;% mean (moment matching)
%                 s2 = V2(:,i-1)*k2 + k3;   % var (moment matching)
%                 psi = s2./m.^2;   % psi compared to psiC
% 
%                 psihat = 1./psi;
%                 b2 = 2*psihat - 1 + sqrt(2*psihat.*(2*psihat-1));
%                 a = m ./ (1 + b2);
% 
%                 % Non-Central Chi squared approximation for psi < psiC
%                 I1 = find(psi<=obj.psiC); 
%                 if isempty(I1)
%                 else
%                     V2(I1,i) = a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
%                 end
%                 p = (psi - 1)./(psi + 1);               % for switchingg
%                 ruleg
%                 V2((UV1(:,i-1)<=p) & (psi>obj.psiC),i) = 0; % case u<=p & psi>psiC
%                 I1b = find((UV1(:,i-1)>p) & (psi>obj.psiC));% find is faster here!
% 
%                 beta = (1 - p)./m;                      % for switching rule
%                 if isempty(I1b)
%                 else    % Psi^(-1)
%                     V2(I1b,i) = log((1-p(I1b))./(1-UV1(I1b,i-1)))./beta(I1b);
%                 end
%                 % log Euler Predictor-Corrector step
%                     lnS1(:,i) = lnS1(:,i-1) + K0 + K1.*V2(:,i-1) + K2.*V2(:,i) +...
%                         sqrt(K3.*V2(:,i-1) + K4.*V2(:,i)).*dW2(:,i-1);
%              end
%             
%              paths(:,:) = exp(lnS1(:,t0+1:end));   
% %               paths = exp(lnS1);
% %             paths(1,:)
%          end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************************************************
%% set U=V-Vlong
             U = zeros(nPaths,Ntime+1);
             U(:,1)=V2(:,1) - obj.assetParam.Vlong; % set U0
             VRing = zeros(nPaths,Ntime+1);
             k1 = exp(-obj.assetParam.kappa*dT);
%              k2 = obj.assetParam.nu^2*k1.*(1-k1)/obj.assetParam.kappa;
% %             k3 = exp(obj.assetParam.kappa*dT)*0.5.*k2.*(1-k1).*obj.assetParam.Vlong;
             for i=2:Ntime+1;%obj.timeDim.nSteps             % time loop
                 m = obj.assetParam.Vlong + U(:,i-1)*k1;% mean (moment matching)
%                 s2 = k2*(U(:,i-1)*k1 + obj.assetParam.Vlong/2*(1+k1));   % var (moment matching)
%                 psi = s2./m.^2;   % psi compared to psiC
                 psi = ((U(:,i-1)+obj.assetParam.Vlong).*exp(-obj.assetParam.kappa*dT)./obj.assetParam.kappa*(1-exp(-obj.assetParam.kappa*dT))+obj.assetParam.Vlong...
                     /2/obj.assetParam.kappa*(1-exp(-2*obj.assetParam.kappa*dT)))./(obj.assetParam.Vlong+U(:,i-1).*exp(-obj.assetParam.kappa*dT)).^2;
%                psihat = 1./psi;
%                b2 = 2*psihat - 1 + sqrt(2*psihat*(2*psihat-1));
%                b2 = 2./psi.*(1-psi/2+sqrt(1-psi/2)); % rewrite b^2
%                 binv = sqrt(1./sqrt(1-psi./2)-1); %change b^2 to 1/b
                 binv2 = psi./(2*sqrt(1-psi.*obj.assetParam.nu^2/2).*(1+sqrt(1-psi*obj.assetParam.nu^2/2)));
%                a = m ./ (1 + b2);

                % Non-Central Chi squared approximation for psi < psiC
                I1 = find(psi<=obj.psiC); 
                if isempty(I1)
                else
                    %U(I1,i) = -obj.assetParam.Vlong + a(I1).*(sqrt(b2(I1)) + norminv(UV1(I1,i-1))).^2;
%                     U(I1,i) = -obj.assetParam.Vlong + m(I1)./(1+binv(I1).^2).*(1+norminv(UV1(I1,i-1)).*binv(I1)).^2;
                    U(I1,i) = -obj.assetParam.Vlong + m(I1)./(1+obj.assetParam.nu^2*binv2(I1)).*(1+obj.assetParam.nu*norminv(UV1(I1,i-1)).*sqrt(binv2(I1))).^2;
                    VRing(I1,i) = m(I1)./obj.assetParam.nu.*((1+obj.assetParam.nu.*sqrt(binv2(I1))...
                        .*norminv(UV1(I1,i-1))).^2./(1+obj.assetParam.nu^2.*binv2(I1))-1);
                end
                p = (psi - 1)./(psi + 1);               % for switching rule
                U((UV1(:,i-1)<=p) & (psi>obj.psiC),i) = -obj.assetParam.Vlong; % case u<=p & psi>psiC
                I1b = find((UV1(:,i-1)>p) & (psi>obj.psiC));% find is faster here!

                beta = (1 - p)./m;                      % for switching rule
                if isempty(I1b)
                else    % Psi^(-1)
                    U(I1b,i) = -obj.assetParam.Vlong + log((1-p(I1b))./(1-UV1(I1b,i-1)))./beta(I1b);
                end
                % log Euler Predictor-Corrector step
                Gammas = (1-exp(-obj.assetParam.kappa*dT))/obj.assetParam.kappa/dT*U(:,i-1)+gamma2*obj.assetParam.nu*VRing(:,i);
                lnS1(:,i) = lnS1(:,i-1) - obj.assetParam.Vlong*dT/2 - dT/2*Gammas + obj.assetParam.rho*(obj.assetParam.kappa*dT*exp(obj.assetParam.kappa*dT)/...
                    (exp(obj.assetParam.kappa*dT)-1)*VRing(:,i))+sqrt(dT*(1-obj.assetParam.rho^2)*(obj.assetParam.Vlong+Gammas)).*dW2(:,i-1);
%                 temp1 = -dT/2*(obj.assetParam.Vlong+1/2*(U(:,i-1)+U(:,i)));
%                 temp2 = obj.assetParam.kappa*dT/2*(U(:,i-1)+U(:,i))+(-U(:,i-1)+U(:,i));
%                 temp3 = 1/2*dT*(1-obj.assetParam.rho^2)*(U(:,i-1)+U(:,i));
%                 temp4 = dT*obj.assetParam.Vlong*(1-obj.assetParam.rho^2);
%                 if temp2 == 0
%                     lnS1(:,i) = lnS1(:,i-1) + temp1 + sqrt(temp4).*dW2(:,i-1);
%                 else
%                     lnS1(:,i) = lnS1(:,i-1) + temp1 + obj.assetParam.rho/obj.assetParam.nu*temp2 + sqrt(temp3+temp4).*dW2(:,i-1);
%                 end
%                     lnS1(:,i) = lnS1(:,i-1) + obj.assetParam.rho*obj.assetParam.kappa/obj.assetParam.nu*dT...
%                         *0.5*(U(:,i-1)+U(:,i)) + obj.assetParam.rho/obj.assetParam.nu*(-U(:,i-1) + U(:,i))...
%                          - dT/4*(U(:,i-1)+U(:,i)) - dT*obj.assetParam.Vlong/2 + sqrt(0.5*dT*(1-obj.assetParam.rho^2)...
%                          *(U(:,i-1)+U(:,i))+dT*obj.assetParam.Vlong*(1-obj.assetParam.rho^2)).*dW2(:,i-1);
             end
            
             paths(:,:) = exp(lnS1(:,2:end));   
%               paths = exp(lnS1);
%             paths(1,:)
         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
   end
  

   methods (Access = protected)

      function propList = getPropertyList(obj)
          if strcmp(obj.assetParam.pathType,'GBM')
            propList = getPropertyList@brownianMotion(obj);
            propList.assetParam_pathType = obj.assetParam.pathType;
            propList.assetParam_initPrice = obj.assetParam.initPrice;
            propList.assetParam_interest = obj.assetParam.interest;
            propList.assetParam_volatility = obj.assetParam.volatility;
            if obj.assetParam.drift ~=0
                propList.assetParam_drift = obj.assetParam.drift;
            end
            propList.assetParam_nAsset = obj.assetParam.nAsset;
         else
            propList = getPropertyList@brownianMotion(obj);
            propList.assetParam_pathType = obj.assetParam.pathType;
            propList.assetParam_initPrice = obj.assetParam.initPrice;
            propList.assetParam_interest = obj.assetParam.interest;
            propList.assetParam_dividend = obj.assetParam.dividend;
            propList.assetParam_Vinst = obj.assetParam.Vinst;
            propList.assetParam_Vlong = obj.assetParam.Vlong;
            propList.assetParam_kappa = obj.assetParam.kappa;
            propList.assetParam_nu = obj.assetParam.nu;
            propList.assetParam_rho = obj.assetParam.rho;
          end              
          
      end

   end

end
