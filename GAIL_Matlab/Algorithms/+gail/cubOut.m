classdef cubOut < gail.cubParam & gail.outParam
   %GAIL.cubOUT is a class containing the parameters related to the
   %outputs from the algorithms that compute integrals.
   %   This class includes the approximation to the integral
   %
   % Example 1.
   % >> cubParamObj = gail.cubParam; %an input object
   % >> cubOutObj = gail.cubOut(cubParamObj); %copied to becom an output object
   % >> cubOutObj.mu = 1.467; %integral value is recorded
   % >> cubOutObj.nSample = 31415926; %sample size is recorded
   % >> cubOutObj.time = 0.0278 %time of computation is recorded
   % cubOutObj = 
   %   cubOut with properties:
   % 
   %             mu: 1.467000000000000
   %        measure: 'uniform'
   %       trueMuCV: [1×0 double]
   %            nMu: 1
   %             nf: 1
   %        inflate: @(m)5*2^-m
   %             ff: @(x)sum(x.^2,2)
   %            nCV: 0
   %         volume: 1
   %              f: @(x)sum(x.^2,2)
   %         domain: [2×1 double]
   %     domainType: 'box'
   %          nInit: 100
   %           nMax: 10000000
   %              d: 1
   %          nfOut: 1
   %         absTol: 0.010000000000000
   %         relTol: 0
   %         solFun: @(mu)mu
   %       solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %        nSample: 31415926
   %           time: 0.027800000000000
   %
   %
   % Author: Fred J. Hickernell

   properties
      mu %approximation to the mean
   end
   
   properties (Hidden, SetAccess = private)
   end
   
   methods
      
      % Creating a meanYParam process
      function obj = cubOut(val)
         %this constructor essentially parses inputs
         %the parser will look for a meanYParam object
         
         obj@gail.cubParam(val)
        
      end %of constructor
     
      function set.mu(obj,val)
         validateattributes(val, {'numeric'}, {'scalar'})
         obj.mu = val;
      end     
      
   end
   
end

