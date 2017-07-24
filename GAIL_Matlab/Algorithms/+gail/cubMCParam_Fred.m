classdef cubMCParam < gail.cubParam 
   %GAIL.CUBMCPARAM is a class containing the parameters related to
   %algorithms that find the mean of a random variable
   %   This class contains the number of integrands with the same integral,
   %   etc.
   %
   % Example 1. Construct a cubParam object with default parameters
   % >> cubParamObj = gail.cubParam
   % cubParamObj = 
   %   cubParam with properties:
   % 
   %              f: @(x)sum(x.^2,2)
   %         domain: [2×1 double]
   %        measure: 'uniform'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 2. Construct a cubParam object with properly ordered inputs
   % >> cubParamObj = gail.cubParam(@(x) sum(x.^3.2),[0 0; 2 2],'box','Lebesgue')
   % cubParamObj = 
   %   cubParam with properties:
   % 
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %        measure: 'Lebesgue'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 3. Using name/value pairs
   % >> cubParamObj = gail.cubParam('domain', [-Inf -Inf; Inf Inf], 'f', @(x) sum(x.^3.2), 'relTol', 0.1, 'measure', 'Gaussian')
   % cubParamObj = 
   %   cubParam with properties:
   % 
   %              f: @(x)sum(x.^3.2)
   %         domain: [2×2 double]
   %        measure: 'normal'
   %         absTol: 0.010000000000000
   %         relTol: 0.100000000000000
   %
   %
   % Example 4. Using a structure for input
   % >> inpStruct.f = @(x) sin(sum(x,2));
   % >> inpStruct.domain = [zeros(1,4); ones(1,4)];
   % >> inpStruct.nInit = 1000;
   % >> cubParamObj = gail.cubParam(inpStruct)
   % cubParamObj = 
   %   cubParam with properties:
   % 
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %        measure: 'uniform'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %
   %
   % Example 5. Copying a cubParam object and changing some properties
   % >> NewCubParamObj = gail.cubParam(cubParamObj,'measure','Lebesgue')
   % NewCubParamObj = 
   %   cubParam with properties:
   % 
   %              f: @(x)sin(sum(x,2))
   %         domain: [2×4 double]
   %        measure: 'Lebesgue'
   %         absTol: 0.010000000000000
   %         relTol: 0
   %          nInit: 1000
   %
   %

   
   % Author: Fred J. Hickernell
   
   properties
   
   end
   
   
   
   methods
      
      % Creating a cubParam process
      function obj = cubMCParam(varargin)
         %this constructor essentially parses inputs
         %the parser will look for the following in order
         %  # a copy of a cubParam object
         %  # a structure
         %  # a function
         %  # a domain
         %  # a measure
         %  # numbers: absTol, relTol, trueMuCV, nMu, nf, inflate
         %  # name-value pairs
         
         obj@gail.cubParam(varargin);

      end
      
   end
  
end 
 
