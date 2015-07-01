classdef optPayoff < assetPath
 
%% optPayoff
% is a class of option payoffs based on asset paths.
% 
% Example 1
% >> obj = optPayoff
% obj = 
%   optPayoff with properties:
% 
%                   inputType: 'n'
%          timeDim_timeVector: [1 2 3]
%           timeDim_startTime: 1
%             timeDim_endTime: 3
%            timeDim_initTime: 0
%           timeDim_initValue: 10
%          wnParam_sampleKind: 'IID'
%         wnParam_distribName: 'Gaussian'
%            wnParam_xDistrib: 'Uniform'
%        bmParam_assembleType: 'diff'
%         assetParam_pathType: 'GBM'
%        assetParam_initPrice: 10
%         assetParam_interest: 0.0100
%       assetParam_volatility: 0.5000
%         payoffParam_optType: {'euro'}
%     payoffParam_putCallType: {'call'}
%          payoffParam_strike: 10
%                  exactPrice: 3.4501
 
%% Properties
% This process inherits properties from the |stochProcess| class.  Below are 
% values assigned to that are abstractly defined in that class plus some
% properties particulary for this class
 
   properties (SetAccess=public) %so they can only be set by the constructor
      payoffParam = struct('optType', {{'euro'}}, ... %type of option
         'putCallType', {{'call'}}, ... %put or call
         'strike', 10, ... %strike price
         'barrier', 12, ... %barrier
         'digitalPay', 100,...%digital payoff 
         'basketWeight',1)%basket constants
      
   end
 
   properties (Constant, Hidden) %do not change & not seen
      allowOptType = {'euro','upin', 'downin' 'upout', 'downout', 'look', ...
         'amean', 'gmean', 'digitalcash', 'digitalasset','basket'} 
         %kinds of payoffs that we can generate
      allowPutCallType = {'call','put'} 
         %kinds of payoffs that we can generate
         
   end
 
   properties (Dependent = true, SetAccess = private)
      exactPrice
   end
   
   properties (SetAccess=private, Hidden) %so they can only be set by the constructor
      defaultNPayoffs = 1e4;
   end
      
%% Methods
% The constructor for |assetPath| uses the |brownianMotion| constructor
% and then parses the other properties. The function |genStockPaths| generates
% the asset paths based on |whiteNoise| paths.
 
   methods
        
      % Creating an asset path process
      function obj = optPayoff(varargin)         
         obj@assetPath(varargin{:}) %parse basic input
         if isfield(obj.restInput,'payoffParam')
            val = obj.restInput.payoffParam;
            obj.payoffParam = val;
            obj.restInput = rmfield(obj.restInput,'payoffParam');
         end
         
      end
      
      % Set the properties of the payoff object
      function set.payoffParam(obj,val)
         if isfield(val,'optType') %data for type of option
            assert(any(strcmp(val.optType,obj.allowOptType)))
            obj.payoffParam.optType=val.optType; %row
         end
         if isfield(val,'putCallType') %data for type of option
            assert(all(any(strcmp( ...
               repmat(val.putCallType,numel(obj.allowPutCallType),1), ...
               repmat(obj.allowPutCallType',1,numel(val.putCallType))),1),2))
            obj.payoffParam.putCallType=val.putCallType; %row
         end
         assert(numel(obj.payoffParam.optType)  ...
            == numel(obj.payoffParam.putCallType))
         if isfield(val,'strike') %data for type of option
            validateattributes(val.strike,{'numeric'}, ...
               {'nonnegative'})
            obj.payoffParam.strike=val.strike; %row
         end
         if isfield(val,'barrier') %data for type of option
            validateattributes(val.barrier,{'numeric'}, ...
               {'nonnegative'})
            obj.payoffParam.barrier=val.barrier; %row
         end
         if isfield(val,'digitalPay') %data for type of option
            validateattributes(val.digitalPay,{'numeric'}, ...
               {'nonnegative'})
            obj.payoffParam.digitalPay=val.digitalPay; %row
         end
         if isfield(val,'basketWeight') %data for type of option
            validateattributes(val.basketWeight,{'numeric'}, ...
               {'nonnegative'})
            obj.payoffParam.basketWeight=val.basketWeight; %row
         end
      end
 
           
      % Generate payoffs of options
      function payoffs=genOptPayoffs(obj,val)
         paths = genPaths(obj,val);
         nOptType = numel(obj.payoffParam.optType);
         nPaths = size(paths,1);
         tempPay = zeros(nPaths, nOptType);
         
         wh=strcmp(obj.payoffParam.optType,'stockprice');
         if any(wh) %final stock price
            tempPay(:,wh)=paths(:,obj.timeDim.nSteps) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whdigitalcashcall = strcmp(obj.payoffParam.optType,'digitalcash') ...
             & strcmp(obj.payoffParam.putCallType,'call'); %digital cash call
             digitalPay=obj.payoffParam.digitalPay;

         if any(whdigitalcashcall) %digitalcash option
             tempPay(:,whdigitalcashcall) ...
               =  (paths(:,obj.timeDim.nSteps) > obj.payoffParam.strike) ...
               .* digitalPay...
               .*exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whdigitalcashput = strcmp(obj.payoffParam.optType,'digitalcash') ...
             & strcmp(obj.payoffParam.putCallType,'put'); %digital cash put
         if any(whdigitalcashput) %digital option
             tempPay(:,whdigitalcashput) ...
               =  (paths(:,obj.timeDim.nSteps) <= obj.payoffParam.strike) ...
               .* digitalPay...
               .*exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whdigitalassetcall = strcmp(obj.payoffParam.optType,'digitalasset') ...
             & strcmp(obj.payoffParam.putCallType,'call'); %digitalasset call

         if any(whdigitalassetcall) %digitalasset option
             tempPay(:,whdigitalassetcall) ...
               =  (paths(:,obj.timeDim.nSteps) > obj.payoffParam.strike) ...
               .* (paths(:,obj.timeDim.nSteps))...
               .*exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whdigitalassetput = strcmp(obj.payoffParam.optType,'digitalasset') ...
             & strcmp(obj.payoffParam.putCallType,'put'); %digital put
         if any(whdigitalassetput) %digital option
             tempPay(:,whdigitalassetput) ...
               =  (paths(:,obj.timeDim.nSteps) <= obj.payoffParam.strike) ...
               .* (paths(:,obj.timeDim.nSteps))...
               .*exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
      
          whbasketcall = strcmp(obj.payoffParam.optType,'basket') ...
             & strcmp(obj.payoffParam.putCallType,'call'); %basket call
             weight=obj.payoffParam.basketWeight;

         if any(whbasketcall) %basket option
             tempPay(:,whbasketcall) ...
               =  max(paths(:,obj.timeDim.nSteps*(1:obj.assetParam.nAsset)) ...
               * weight' - obj.payoffParam.strike, 0) ...
               * exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whbasketput= strcmp(obj.payoffParam.optType,'basket') ...
             & strcmp(obj.payoffParam.putCallType,'put'); %basket put
         if any(whbasketput) %basket option
             tempPay(:,whbasketput) ...
               =  max(-paths(:,obj.timeDim.nSteps*(1:obj.assetParam.nAsset)) ...
               * weight' + obj.payoffParam.strike, 0) ...
               * exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         wheurobarrier = any(strcmp(repmat(obj.payoffParam.optType,5,1), ...
            repmat({'euro','upin','upout','downin','downout'}',1,nOptType)),1);
         wheurobarriercall = wheurobarrier  ...
            & strcmp(obj.payoffParam.putCallType,'call');
         if any(wheurobarriercall) %call payoff
            tempPay(:,wheurobarriercall) ...
               =  max(paths(:,obj.timeDim.nSteps) ...
               - obj.payoffParam.strike, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
 
         wheurobarrierput = wheurobarrier  ...
            & strcmp(obj.payoffParam.putCallType,'put');
         if any(wheurobarrierput); %put payoff
            tempPay(:,wheurobarrierput) =  max(obj.payoffParam.strike ...
               - paths(:,obj.timeDim.nSteps), 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         wh=strcmp(obj.payoffParam.optType,'upin');
         if any(wh); %up and in barrier
            if obj.assetParam.initPrice < obj.payParam.barrier;
               tempPay(:,wh) = tempPay(:,wh) ...
                  .* any(paths >= obj.payParam.barrier,2);
            end
         end
 
         wh=strcmp(obj.payoffParam.optType,'downin');
         if any(wh); %down and in barrier
            if obj.assetParam.initPrice > obj.payParam.barrier;
               tempPay(:,wh) = tempPay(:,wh) ...
                  .* any(paths <= obj.payParam.barrier,2);
            end
         end
 
         wh=strcmp(obj.payoffParam.optType,'upout');
         if any(wh); %up and out barrier
            if obj.assetParam.initPrice < obj.payParam.barrier;
               tempPay(:,wh) = tempPay(:,wh) ...
                  .* all(paths < obj.payParam.barrier,2);
            else
               tempPay(:,wh) = zeros(nPaths,sum(wh));
            end
         end
 
         wh=strcmp(obj.payoffParam.optType,'downout');
         if any(wh); %down and out barrier
            if obj.assetParam.initPrice > obj.payParam.barrier;
               tempPay(:,wh) = tempPay(:,wh) ...
                  .* all(paths > obj.payParam.barrier,2);
            else
               tempPay(:,wh) = zeros(nPaths,sum(wh));
            end
         end
         
         whlook = strcmp(obj.payoffParam.optType,'look');
         wh = whlook & strcmp(obj.payoffParam.putCallType,'call');
         if any(wh)
            K = min([repmat(obj.assetParam.initPrice,nPaths,1) paths],[ ],2);
            tempPay(:,wh) ...
               =  max(paths(:,obj.timeDim.nSteps) - K, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
 
         wh = whlook & strcmp(obj.payoffParam.putCallType,'put');
         if any(wh)
            K = min([repmat(obj.assetParam.initPrice,nPaths,1) paths],[ ],2);
            tempPay(:,wh) ...
               =  max(K - paths(:,obj.timeDim.nSteps), 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
 
         whamean=strcmp(obj.payoffParam.optType,'amean');
         if any(whamean); %arithmetic mean
            meanstock=mean(paths,2);
         end  
         
         wh = whamean & strcmp(obj.payoffParam.putCallType,'call');
         if any(wh); %arithmetic mean call
            tempPay(:,wh)=max(meanstock - obj.payoffParam.strike, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         wh = whamean & strcmp(obj.payoffParam.putCallType,'put');
         if any(wh); %arithmetic mean put
            tempPay(:,wh)=max(obj.payoffParam.strike - meanstock, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         whgmean=strcmp(obj.payoffParam.optType,'gmean');
         if any(whgmean); %geometric mean
            meanstock=prod(paths,2).^(1./obj.timeDim.nSteps);
         end  
         
         wh = whgmean & strcmp(obj.payoffParam.putCallType,'call');
         if any(wh); %geometric mean call
            tempPay(:,wh)=max(meanstock - obj.payoffParam.strike, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         wh = whgmean & strcmp(obj.payoffParam.putCallType,'put');
         if any(wh); %geometric mean put
            tempPay(:,wh)=max(obj.payoffParam.strike - meanstock, 0) ...
               .* exp(- obj.assetParam.interest .* obj.timeDim.endTime);
         end
         
         payoffs = tempPay;
      end
      
      function val = get.exactPrice(obj)
         %Expected value of ending price of asset
         val = NaN(1,numel(obj.payoffParam.optType));
         wh = strcmp('stockprice',obj.payoffParam.optType);
         if any(wh); 
            val(wh)=obj.assetParam.initPrice * ...
               exp(obj.assetParam.interest * ...
               obj.timeDim.endTime);
         end
 
         %Pricing European geometric brownian motion
         wheuro = strcmp(obj.payoffParam.optType, 'euro');
         whcall = strcmp(obj.payoffParam.putCallType, 'call');
         whput = strcmp(obj.payoffParam.putCallType, 'put');
         wheurocall = wheuro & whcall;
         wheuroput = wheuro & whput;
         if any(wheuro); 
            [eurocall,europut] = eurogbmprice(obj.assetParam.initPrice, ...
               obj.assetParam.interest, obj.timeDim.endTime, ...
               obj.assetParam.volatility, obj.payoffParam.strike);
            val(wheurocall) = eurocall;
            val(wheuroput) = europut;
         end
         
         %Pricing digital option
         whdigitcash = strcmp(obj.payoffParam.optType, 'digitalcash');
         whdigitasset = strcmp(obj.payoffParam.optType, 'digitalasset');
         whcall = strcmp(obj.payoffParam.putCallType, 'call');
         whput = strcmp(obj.payoffParam.putCallType, 'put');
         whdigitcashcall = whdigitcash & whcall;
         whdigitassetcall= whdigitasset & whcall;
         whdigitcashput = whdigitcash & whput;
         whdigitassetput = whdigitasset & whput;
         if any(whdigitcash | whdigitasset); 
            [digitcashcall,digitassetcall,digitcashput,digitassetput] ...
             = digitgbmprice(obj.payoffParam.digitalPay,obj.assetParam.initPrice, ...
               obj.assetParam.interest, obj.timeDim.endTime, ...
               obj.assetParam.volatility, obj.payoffParam.strike,0);
            val(whdigitcashcall) = digitcashcall;
            val(whdigitassetcall) = digitassetcall;
            val(whdigitcashput) = digitcashput;
            val(whdigitassetput) = digitassetput;
         end
 
         %Pricing Asian geometric mean
         whgmean = strcmp(obj.payoffParam.optType, 'gmean');
         whgmeancall = whgmean & whcall;
         whgmeanput = whgmean & whput;
         if any(whgmean) 
            Tbar=(1+1/obj.timeDim.nSteps) * obj.timeDim.endTime / 2; 
            sigmabar=obj.assetParam.volatility * sqrt((2 + 1 ...
               / obj.timeDim.nSteps) / 3);
            rbar = obj.assetParam.interest + (sigmabar^2 ...
               - obj.assetParam.volatility^2) / 2;
            [gmeancall,gmeanput]=eurogbmprice(obj.assetParam.initPrice, ...
               rbar,Tbar,sigmabar,obj.payoffParam.strike);
            gmeancall=gmeancall * exp(rbar * Tbar ...
               - obj.assetParam.interest.*obj.timeDim.endTime);
            gmeanput=gmeanput * exp(rbar * Tbar ...
               - obj.assetParam.interest.*obj.timeDim.endTime);
            val(whgmeancall) = gmeancall;
            val(whgmeanput) = gmeanput;
         end
         
           
         function [callprice,putprice]=eurogbmprice(S0,r,T,sigma,K)
            priceratio = K * exp(-r * T) ./ S0;
            xbig = log(priceratio) ./ (sigma * sqrt(T)) + sigma * sqrt(T)/2;
            xsmall = log(priceratio) ./ (sigma * sqrt(T)) - sigma * sqrt(T)/2;
            putprice = S0 .* (priceratio.*normcdf(xbig) - normcdf(xsmall));
            callprice = putprice + S0 .* (1-priceratio);        
         end
         
         function [digitcashcall,digitassetcall,digitcashput, ...
                 digitassetput] = digitgbmprice(digitalPay,S0,r,T,sigma,K,q)
            digitpriceratio1 = (log(S0./K)+(r-q+(sigma^2)/2)*T)/(sigma*sqrt(T));
            digitpriceratio2 = digitpriceratio1-sigma*sqrt(T);
            digitcashcall = digitalPay*exp(-r*T)*normcdf(digitpriceratio2);
            digitcashput = digitalPay*digitalPay*exp(-r*T)*normcdf(-digitpriceratio2);
            digitassetcall = S0.*exp(-q*T).*normcdf(digitpriceratio1);
            digitassetput = S0.*exp(-q*T).*normcdf(-digitpriceratio1);
         end

      end
      
      function varargout = plot(obj,varargin)
         assert(strcmp(obj.inputType,'n'), ...
            'plot requires inputType to be ''n''')
         if numel(varargin)
            nPayoffs = varargin{1};
         else
            nPayoffs = obj.defaultNPayoffs; %default 
         end
         payoffs = genOptPayoffs(obj,nPayoffs);
         probs = (1/(2*nPayoffs)):(1/nPayoffs):(1 - 1/(2*nPayoffs));
         h = plot(sort(payoffs),probs,'-');
         if numel(varargin) > 1
            set(h,varargin{offset+2:end});
         else
            set(h,obj.defaultLineSpecs{:});
         end
         set(gca,'fontsize',20)
         if nargout
            varargout{1}=h;
         end
         xlabel('payoff')
         ylabel('probability')
      end
       
   end
    
   methods (Access = protected)
 
      function propList = getPropertyList(obj)
         propList = getPropertyList@assetPath(obj);
         propList.payoffParam_optType = obj.payoffParam.optType;
         propList.payoffParam_putCallType = obj.payoffParam.putCallType;
         propList.payoffParam_strike = obj.payoffParam.strike;
             if any(strncmp(obj.payoffParam.optType,'digital',7))
                propList.payoffParam_digitalPay = obj.payoffParam.digitalPay;
             end
         propList.exactPrice = obj.exactPrice;
      end
 
   end
end
 


