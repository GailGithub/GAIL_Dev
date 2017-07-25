classdef cubMCOut < gail.cubMCParam & gail.outParam
   %GAIL.MEANYOUT is a class containing the parameters related to the
   %outputs from the algorithms that find the mean of a random variable.
   %   This class includes the time and sample size required for the
   %   computation
   %
   % Example 1.
   % >> meanParamObj = gail.meanYParam; %an input object
   % >> meanOutObj = gail.meanYOut(meanParamObj); %copied to becom an output object
   % >> meanOutObj.mu = 1.467; %integral value is recorded
   % >> meanOutObj.stddev = 1.23; %standard deviation is recorded
   % >> meanOutObj.nSample = 31415926; %sample size is recorded
   % >> meanOutObj.time = 0.0278; %time of computation is recorded
   % >> meanOutObj.errBd = 0.000456 %error bound is recorded
   % meanOutObj = 
   %   meanYOut with properties:
   % 
   %           Y: @(n)rand(n,1)
   %      absTol: 0.0100
   %      relTol: 0
   %       alpha: 0.0100
   %          mu: 1.4670
   %      stddev: 1.2300
   %     nSample: 31415926
   %        time: 0.0278
   %       errBd: 4.5600e-04
   %
   %
   % Author: Fred J. Hickernell

   
   properties
      stddev %sample standard deviation of the random variable
   end
   
   properties (Hidden, SetAccess = private)
   end
   
   methods
      
      % Creating a meanYParam process
      function obj = cubMCOut(val)
         %this constructor essentially parses inputs
         %the parser will look for a meanYParam object
         
         obj@gail.cubMCParam(val)
        
      end %of constructor
          
      function set.stddev(obj,val)
         validateattributes(val, {'numeric'}, {'nonnegative'})
         obj.stddev = val;
      end                    
          
      
   end
   
   methods (Access = protected)
   
         function propList = getPropertyList(obj)
         propList = getPropertyList@gail.cubParam(obj);
         propList.sol = obj.sol;
         propList.stddev = obj.stddev;
         propList.nSample = obj.nSample;
         propList.time = obj.time;
         propList.errBd = obj.errBd;
         
         propList.alpha = obj.alpha;
         propList.nSig = obj.nSig;
       
         end
   end
   
end

