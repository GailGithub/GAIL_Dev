classdef meanYOut < gail.meanYParam & gail.outParam
   %GAIL.MEANYOUT is a class containing the parameters related to the
   %outputs from the algorithms that find the mean of a random variable.
   %   This class includes the time and sample size required for the
   %   computation
   %
   % Example 1.
   % >> meanParamObj = gail.meanYParam; %an input object
   % >> meanOutObj = gail.meanYOut(meanParamObj); %copied to becom an output object
   % >> meanOutObj.mu = 1.467; %integral value is recorded
   % >> meanOutObj.nSample = 31415926; %sample size is recorded
   % >> meanOutObj.time = 0.0278 %time of computation is recorded
   % meanOutObj = 
   %   meanYOut with properties:
   % 
   %           mu: 1.467000000000000
   %          std: []
   %            Y: @(n)rand(n,1)
   %        alpha: 0.010000000000000
   %         nSig: 1000
   %      inflate: 1.200000000000000
   %         nMax: 100000000
   %          nMu: 1
   %           nY: 1
   %     trueMuCV: [1×0 double]
   %          nCV: 0
   %        nYOut: 1
   %       absTol: 0.010000000000000
   %       relTol: 0
   %       solFun: @(mu)mu
   %     solBdFun: @(muhat,errbd)[muhat-errbd,muhat+errbd]
   %      nSample: 31415926
   %         time: 0.027800000000000
   %
   %
   % Author: Fred J. Hickernell

   
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

