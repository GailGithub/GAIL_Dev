classdef cubBayesLatticeOut < gail.cubBayesLatticeParam & gail.outParam
   %GAIL.CUBBAYESLATTICEOUT is a class containing the parameters related to the
   %outputs from the cubBayesLattice algorithm that computes integrals.
   %   This class includes the approximation to the integral
   %
   % Example 1.
   % >> cubBayesLatticeParamObj = gail.cubBayesLatticeParam; %an input object
   % >> cubBayesLatticeOutObj = gail.cubBayesLatticeOut(cubBayesLatticeParamObj); %copied to becom an output object
   % >> cubBayesLatticeOutObj.mu = 1.467; %integral value is recorded
   % >> cubBayesLatticeOutObj.nSample = 31415926; %sample size is recorded
   % >> cubBayesLatticeOutObj.time = 0.0278 %time of computation is recorded
   % cubBayesLatticeOutObj = ***
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
   %          errBd: []
   %         tolVal: []
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
      function obj = cubBayesLatticeOut(val)
         %this constructor essentially parses inputs
         %the parser will look for a meanYParam object

         obj@gail.cubBayesLatticeParam(val)

      end %of constructor

      function set.mu(obj,val)
         validateattributes(val, {'numeric'}, {'scalar'})
         obj.mu = val;
      end

   end

    methods (Access = protected)

         function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubBayesLatticeParam(obj);
         propList.mu = obj.mu;
         propList.nSample = obj.nSample;
         propList.time = obj.time;
         propList.errBd = obj.errBd;
         propList.tolVal = obj.tolVal;
         end
   end

end
