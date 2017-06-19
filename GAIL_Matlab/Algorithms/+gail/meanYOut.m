classdef meanYOut < gail.meanYParam & gail.outParam
   %GAIL.MEANYOUT is a class containing the parameters related to the
   %outputs from the algorithms that find the mean of a random variable.
   %   This class includes the time and sample size required for the
   %   computation
   
   properties
      mu %approximation to the mean
      std %sample standard deviation of the random variable
   end
   
   properties (Hidden, SetAccess = private)
   end
   
   methods
      
      % Creating a meanYParam process
      function obj = meanYOut(val)
         %this constructor essentially parses inputs
         %the parser will look for a meanYParam object
         
         obj@gail.meanYParam(val)
        
      end %of constructor
     
      function set.mu(obj,val)
         validateattributes(val, {'numeric'}, {'scalar'})
         obj.mu = val;
      end                    
     
      function set.std(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative'})
         obj.std = val;
      end                    
     
      
   end
   
end

