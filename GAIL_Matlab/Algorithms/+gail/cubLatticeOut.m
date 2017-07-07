classdef cubLatticeOut < gail.cubLatticeParam & gail.outParam
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
   % cubBayesLatticeOutObj = 
   %   cubBayesLatticeOut with properties:
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

