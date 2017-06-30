classdef cubLatticeOut < gail.cubLatticeParam & gail.outParam
   %GAIL.CUBLATTICEOUT is a class containing the parameters related to the
   %outputs from the cubBayesLattice algorithm that computes integrals.
   %   This class includes the approximation to the integral
   %
   % Example 1.
   % >> cubLatticeParamObj = gail.cubLatticeParam; %an input object
   % >> cubLatticeParamObj = gail.cubLatticeOut(cubLatticeParamObj); %copied to becom an output object
   % >> cubLatticeParamObj.mu = 1.467; %integral value is recorded
   % >> cubLatticeParamObj.nSample = 31415926; %sample size is recorded
   % >> cubLatticeParamObj.time = 0.0278 %time of computation is recorded
   % cubLatticeParamObj = 
   %   cubLatticeOut with properties:
   % 
   %              f: @(x)sum(x.^2,2)
   %         domain: [2×1 double]
   %        measure: 'uniform'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %             mu: 1.467000000000000
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
      
      % Creating a cubOut process
      function obj = cubLatticeOut(val)
         %this constructor essentially parses inputs
         %the parser will look for a meanYParam object
         
         obj@gail.cubLatticeParam(val)
        
      end %of constructor
     
      function set.mu(obj,val)
         validateattributes(val, {'numeric'}, {'scalar'})
         obj.mu = val;
      end     
      
   end

    methods (Access = protected)
   
         function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubLatticeParam(obj);
         propList.mu = obj.mu;
         propList.nSample = obj.nSample;
         propList.time = obj.time;
         propList.errBd = obj.errBd;
         propList.tolVal = obj.tolVal;
         end
   end

end

