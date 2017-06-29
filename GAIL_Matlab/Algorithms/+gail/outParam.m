classdef outParam < handle
   %GAIL.OUTPARAM is a class containing the outputs for GAIL algorithms
   %   This includes the time and sample size required for the
   %   computation
   
   properties
      nSample = 0 %total sample size
      time = 0 %time required for computation
      errBd %error bounds on the means of Y our the integrals
      tolVal %value of the tolerance function
   end
      
   methods
          
      function set.nSample(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','positive','integer'})
         obj.nSample = val;
      end
      
      function set.time(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative'})
         obj.time = val;
      end                     
     
      function set.errBd(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative'})
         obj.errBd = val;
      end                    
     
      function set.tolVal(obj,val)
         validateattributes(val, {'numeric'}, {'scalar','nonnegative'})
         obj.tolVal = val;
      end                    
     
   end
   
end

