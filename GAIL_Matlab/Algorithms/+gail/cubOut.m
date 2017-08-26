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
   % cubOutObj = ***
   %
   %              f: @(x)sum(x.^2,2)
   %         domain: [2***1 double]
   %    measureType: 'uniform'
   %        measure: 'uniform'
   %         absTol: 0.0100
   %         relTol: 0
   %             mu: 1.4670
   %        nSample: 31415926
   %           time: 0.0278
   %
   %
   % Author: Fred J. Hickernell

   properties
      mu %approximation to the mean
   end

   properties (Hidden, SetAccess = private)
   end

   methods

      % Creating a cubOut process
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

    methods (Access = protected)

         function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubParam(obj);
         propList.mu = obj.mu;
         propList.nSample = obj.nSample;
         propList.time = obj.time;
         propList.errBd = obj.errBd;
         end
   end

end
