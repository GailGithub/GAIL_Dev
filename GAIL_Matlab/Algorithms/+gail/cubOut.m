classdef cubOut < gail.cubParam & gail.outParam
   %GAIL.cubOUT is a class containing the parameters related to the
   %outputs from the algorithms that compute integrals.
   %   This class includes the approximation to the integral
   
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

